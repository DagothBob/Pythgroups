import cProfile
from io import StringIO
from typing import List, Dict, Optional, ValuesView, Iterator, TextIO
import sys, time, threading

import yaml
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from networkx import Graph

import InputPreprocessing
from BPGDistance import BPGDistance
from DCJOperation import OperationTypes
from DCJRearrangement import DCJRearrangement
from Genome import Genome, split_at_whitespace
from GenomeInString import GenomeInString
from GroupGraph import GroupGraph
from MedianIteration import MedianIteration
from PGMPathForAGenome import PGMPathForAGenome
from SmallPhylogeny import SmallPhylogeny
from TreeStructure import TreeStructure

"""
 Driver program for Pythgroups
 
 Command line usage: 
 $ python [directory of GenomeReconstruction.py] [algorithm]
 where [algorithm] is the desired genome reconstruction algorithm:
    - SmallPhylogeny
    - GenomeAliquoting
    - DCJRearrangements
    - GenomeHalving
 
 Based on the small phylogeny program developed by C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""

CONFIG_DIR = "config.yaml"
CONFIG_GENOME_FILE = "genome_file"
CONFIG_ALGORITHM = "algorithm"
CONFIG_TREE_STRUCTURE = "tree_structure"
CONFIG_OPERATIONS = "operations"
CONFIG_MINIMUM_CHROMOSOME = "minimum_chromosome"
CONFIG_MAXIMUM_CHROMOSOME = "maximum_chromosome"
CONFIG_WHICH_CHROMOSOME = "which_chromosome"
CONFIG_NUMBER_OPERATIONS = "number_of_operations"
CONFIG_GENOME_REPLACE = "genome_to_replace"

done = False  # used in the progress spinner


class NetworkxNode:
    """
    Represents a node in the networkx graph, containing additional information like a unique integer identifier

    Attributes
    ----------
    name : str
        Name of the genome
    neighbors : List[int]
        The genome_ids of this genome's neighbors
    genome_id : int
        Integer identifier for this genome
    """
    def __init__(self, name: str, neighbors: List[int], genome_id: int):
        """
        Constructor

        Attributes
        ----------
        name : str
            Name of the genome
        neighbors : List[int]
            The genome_ids of this genome's neighbors
        genome_id : int
            Integer identifier for this genome
        """
        self.name: str = name
        self.neighbors: List[int] = neighbors
        self.genome_id: int = genome_id


class SpinnerThread(threading.Thread):
    """
    Used for the progress spinner
    """

    def __init__(self):
        super().__init__(target=self._spin)
        self._stopevent = threading.Event()

    def stop(self):
        self._stopevent.set()

    def _spin(self):
        while not self._stopevent.isSet():
            for t in '|/-\\':
                print("\rPerforming initial reconstruction... {}".format(t), end="")
                time.sleep(0.4)
        print("\bComplete", end="")


def parse_genomes(config_dir: str) -> Dict[str, List[str]]:
    """
    Parses the genome data from a file into a dictionary containing each genome and their chromosomes

    Parameters
    ----------
    config_dir : str
        Directory of the config file from which to get the genome file directory

    Returns
    -------
    Dict[str, List[str]]
        A dictionary containing each genome's chromosomes,
        where keys are genome headers and values are lists of chromosomes
    """
    config_file = open(config_dir, "r")
    config_data = yaml.safe_load(config_file)
    genome_file = open(config_data.get(CONFIG_GENOME_FILE))
    genomes: Dict[str, List[str]] = dict()
    line: str = genome_file.readline()
    header: str = str()

    while len(line) != 0:

        # ">" indicates a header, sets that as the current genome to add chromosomes to
        if line.startswith(">"):
            header = line.replace(">", "").replace("\n", "")

        # Add chromosomes to the current genome until an empty line, another header, or the end of file is found
        else:
            chromosomes: List[str] = list()
            while len(line) != 0 and not line.startswith(">") and line != "\n":
                clean_input: str = line.replace("$", "").replace("\n", "")
                chromosome: str = ""
                for string in clean_input.strip().split(" "):
                    if string.strip() != "":
                        chromosome += string.strip() + " "

                chromosomes.append(chromosome)
                line = genome_file.readline()
            genomes[header] = chromosomes
        line = genome_file.readline()

    config_file.close()

    return genomes


def parse_tree(config_dir: str) -> Optional[Tree]:
    """
    Parses the newick tree structure given in config.yaml

    Parameters
    ----------
    config_dir : str
        Directory of config.yaml containing the tree structure

    Returns
    -------
    Optional[Tree]
        Tree object converted from the newick tree found in the yaml file
    """
    config_file = open(config_dir, "r")
    config_data = yaml.safe_load(config_file)
    tree_data = Phylo.read(StringIO(config_data.get(CONFIG_TREE_STRUCTURE)), "newick")
    config_file.close()
    return tree_data


def parse_medians(graph_nodes: List[NetworkxNode]) -> List[NetworkxNode]:
    """
    Get the median nodes from a list of GenomeNodes

    Parameters
    ----------
    graph_nodes : List[NetworkxNode]
        List of GraphNodes to parse the medians from

    Returns
    -------
    Dict[Clade, List[Clade]]
        Key: median to reconstruct, value: list of 3 neighboring clades
    """
    medians: List[NetworkxNode] = list()
    for node in graph_nodes:
        if (len(node.neighbors)) == 3:
            medians.append(NetworkxNode(node.name, node.neighbors, node.genome_id))
        # TODO: Exception handling for illegal tree structures (eg. node with > 3 neighbors)

    return medians


def genome_nodes_from_tree(tree: Tree) -> List[NetworkxNode]:
    """
    Converts a tree into a networkx graph consisting of GenomeNode representations of each node

    Parameters
    ----------
    tree : Tree
        Tree object to convert into a graph

    Returns
    -------
    List[NetworkxNode]
        A list of GenomeNodes representing each node, their neighbors, and their unique IDs
    """
    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()
    graph_indexes = dict()
    node_index = 0
    for t in terminals:
        graph_indexes[t.name] = node_index
        node_index += 1
    for nt in nonterminals:
        if nt.name is not None and len(nt.clades) > 0:
            graph_indexes[nt.name] = node_index
            node_index += 1

    graph: Graph = Graph()
    for nt in nonterminals:
        if nt.name is not None and len(nt.clades) > 0:
            for clade in nt.clades:
                graph.add_edge(nt.name, clade.name)

    graph_nodes = list()
    for node in graph:
        genome_id = graph_indexes[node]
        neighbor_ids = list()
        for neighbor in graph.neighbors(node):
            neighbor_ids.append(graph_indexes[neighbor])
        graph_nodes.append(NetworkxNode(node, neighbor_ids, genome_id))

    return graph_nodes


def count_genes(genomes: Dict[str, List[str]]) -> int:
    """
    Counts the number of genes in each chromosome.
    Genomes must all have the same number of genes as compared to the first genome in genomes

    Parameters
    ----------
    genomes : Dict[str, List[str]]
        Dictionary of all the genomes as produced from parse_genomes()

    Returns
    -------
    int
        The number of genes in each genome
    """
    final_count: int = 0
    for genome, chromosomes in genomes.items():
        count: int = 0
        for chromosome in chromosomes:
            count += len(split_at_whitespace(chromosome))

        # Each genome must have the same number of genes (as compared to the first genome in genomes)
        if genome == list(genomes.keys())[0]:
            final_count = count
        elif final_count != count:
            raise Exception("Different numbers of genes in the genomes (check [{}] or the first genome for anomalies). \
                            Exiting.\n".format(genome))

    return final_count


def small_phylogeny():
    """
    Reconstructs the ancestor(s) of known modern genomes given an unrooted binary phylogenetic tree
    using the Pathgroups algorithm.

    """
    tree: Tree = parse_tree(CONFIG_DIR)
    genomes: Dict[str, List[str]] = parse_genomes(CONFIG_DIR)
    genome_nodes: List[NetworkxNode] = genome_nodes_from_tree(tree)
    median_nodes: List[NetworkxNode] = parse_medians(genome_nodes)

    num_ancestor = len([node for node in tree.get_nonterminals() if node.name is not None])
    num_leaves = len([node for node in tree.get_terminals() if node.name is not None])  # don't count unnamed nodes
    num_genes = count_genes(genomes)
    all_genomes: List[GenomeInString] = list()
    for chromosomes in genomes.values():
        all_genomes.append(GenomeInString(Genome.from_strings(chromosomes)))

    ts: TreeStructure = TreeStructure(num_ancestor, num_leaves, num_genes, None, None, None, all_genomes)

    for median in median_nodes:
        ts.set_tree_structure(median.genome_id, median.neighbors[0],
                              median.neighbors[1], median.neighbors[2])

    sp: SmallPhylogeny = SmallPhylogeny(ts)

    task = threading.Thread(target=sp.get_result)
    task.start()

    # progress spinner for before optimization step
    spin_thread = SpinnerThread()
    spin_thread.start()

    task.join()
    spin_thread.stop()
    spin_thread.join()
    print("\n")

    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        print("before optimization, reconstructed ancestors")
        for j in range(0, len(ts.medians[i].median)):
            print("chr {}\n {}".format(j, ts.medians[i].median[j]))

    reconstructed_paths: List[PGMPathForAGenome] = []

    # Leaf paths added first
    if ts.number_of_leaves >= 0:
        reconstructed_paths = [ts.all_paths[i] for i in range(0, ts.number_of_leaves)]

    # Ancestor paths added second
    for i in range(ts.number_of_leaves, ts.number_of_leaves + ts.number_of_ancestors):
        chromosome_strings: List[str] = ts.medians[i - ts.number_of_leaves].median
        median_genome: Genome = Genome.from_strings(chromosome_strings)
        reconstructed_paths.append(PGMPathForAGenome(ts.get_pgm_path(median_genome, i)))

    relation: List[List[int]] = ts.get_relation()
    reconstructed_dist: int = int()
    before_optimization: str = str()

    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[Dict[str, int]] = reconstructed_paths[i].paths
                p2: List[Dict[str, int]] = reconstructed_paths[j].paths
                cur_dist: int = ts.medians[0].get_distance(p1, p2)
                reconstructed_dist += cur_dist
                before_optimization += "  d({r},{c})={d}".format(r=str(i), c=str(j), d=str(cur_dist))

    print("before optimization:\n" + before_optimization)
    print("total distance: " + str(reconstructed_dist))
    print("")

    # Section 4: optimize the result

    mi: MedianIteration = MedianIteration(ts.number_of_leaves, ts.number_of_ancestors, ts.gene_number,
                                          ts.leaves, reconstructed_paths, ts.medians, ts.node_int, ts.node_string)
    mi.optimize_result(1, 50)
    after_optimization: str = str()
    optimized_dist: int = int()

    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[Dict[str, int]] = reconstructed_paths[i].paths
                p2: List[Dict[str, int]] = reconstructed_paths[j].paths
                cur_dist: int = ts.medians[0].get_distance(p1, p2)
                optimized_dist += cur_dist
                after_optimization += "  d({r},{c})={d}".format(r=str(i), c=str(j), d=str(cur_dist))

    print("after optimization:\n" + after_optimization)
    print("total distance: " + str(optimized_dist))

    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        print("reconstructed ancestors:")
        for j in range(0, len(ts.medians[i].median)):
            print("chr {}\n {}".format(j, ts.medians[i].median[j]))


def genome_aliquoting():
    """
    See the 2010 paper, section 2.5

    """
    gene_data = InputPreprocessing.parse_gene_file(CONFIG_DIR)

    genome_data = InputPreprocessing.group_genomes(gene_data)

    pathgroups_data = InputPreprocessing.to_pathgroups_format(genome_data)

    print(pathgroups_data)


def dcj_rearrangements():
    """
    Performs double cut-and-join (DCJ) operations until the source genome matches the target genome,
    records the state of the intermediate genomes and their distances to the target genome at each step

    Performs DCJ operations on the first 2 genomes found in the genome file
    """
    # Create the dictionary of genomes from the input file
    genomes: Dict[str, List[str]] = parse_genomes(CONFIG_DIR)

    # Get the first 2 genomes from the input file
    values_view = genomes.values()
    value_iterator = iter(values_view)
    genome1: List[str] = next(value_iterator)
    genome2: List[str] = next(value_iterator)

    # Calculate the initial BPG distance
    bpg_dist: BPGDistance = BPGDistance(Genome.from_strings(genome1), Genome.from_strings(genome2))
    bpg_dist.calculate_distance()
    cur_dist: int = bpg_dist.distance

    # Get DCJ configuration options from config file
    config_file: TextIO = open(CONFIG_DIR)
    config_data = yaml.safe_load(config_file)
    config_file.close()
    operation_types: List[int] = list()

    for operation in config_data.get(CONFIG_OPERATIONS):
        if operation.lower() == "inversion" and OperationTypes.INVERSION not in operation_types:
            operation_types.append(OperationTypes.INVERSION)
        elif operation.lower() == "translocation" and OperationTypes.TRANSLOCATION not in operation_types:
            operation_types.append(OperationTypes.TRANSLOCATION)
        elif operation.lower() == "fission" and OperationTypes.FISSION not in operation_types:
            operation_types.append(OperationTypes.FISSION)
        elif operation.lower() == "fusion" and OperationTypes.FUSION not in operation_types:
            operation_types.append(OperationTypes.FUSION)
        else:
            raise Exception("Invalid operation type: {}\n"
                            "Valid types: [ inversion | translocation | fission | fusion ]\n".format(operation))

    minimum_chromosome: int = config_data.get(CONFIG_MINIMUM_CHROMOSOME)
    if type(minimum_chromosome) is not int:
        raise Exception("Config attribute \"minimum_chromosome\" must be a number.\n")

    maximum_chromosome: int = config_data.get(CONFIG_MAXIMUM_CHROMOSOME)
    if type(maximum_chromosome) is not int:
        raise Exception("Config attribute \"maximum_chromosome\" must be a number.\n")

    which_chromosome: int = config_data.get(CONFIG_WHICH_CHROMOSOME)
    if type(which_chromosome) is not int:
        raise Exception("Config attribute \"which_chromosome\" must be a number.\n")

    number_operations: int = config_data.get(CONFIG_NUMBER_OPERATIONS)
    if type(number_operations) is not int:
        raise Exception("Config attribute \"number_of_operations\" must be a number.\n")

    # Perform DCJ operations until distance == 0 or there are no more DCJ operations to perform
    dcj: DCJRearrangement = DCJRearrangement(Genome.from_strings(genome1), Genome.from_strings(genome2))
    dcj.initial_value()
    more: bool = True

    while cur_dist > 0 and more:
        more = False
        rearrange_state: List[GenomeInString] = dcj.get_result(minimum_chromosome,
                                                               maximum_chromosome,
                                                               which_chromosome,
                                                               operation_types,
                                                               number_operations)

        if len(rearrange_state) != 0:
            more = True
            for genome in rearrange_state:
                print("*******")
                for i in range(len(genome.chromosomes)):
                    print("Chromosome " + str(i) + "\n" + str(genome.chromosomes[i]))

            new_genome: GenomeInString = rearrange_state[len(rearrange_state) - 1]
            bpg_dist = BPGDistance(Genome(new_genome.chromosomes), Genome.from_strings(genome2))
            bpg_dist.calculate_distance()
            cur_dist = bpg_dist.distance
            print("a run, steps: " + str(len(rearrange_state)) + ", cur_dist: " + str(cur_dist))
            dcj = DCJRearrangement(Genome(new_genome.chromosomes), Genome.from_strings(genome2))
            dcj.initial_value()


def genome_halving():
    # Create the dictionary of genomes from the input file
    genomes: Dict[str, List[str]] = parse_genomes(CONFIG_DIR)

    # Get the first 2 genomes from the input file
    values_view: ValuesView[List[str]] = genomes.values()
    value_iterator: Iterator[List[str]] = iter(values_view)
    tetrad: List[str] = next(value_iterator)
    outgroup: List[str] = next(value_iterator)

    # Get GenomeHalving configuration options
    config_file: TextIO = open(CONFIG_DIR)
    config_data = yaml.safe_load(config_file)
    config_file.close()

    to_replace: int = config_data.get(CONFIG_GENOME_REPLACE)
    if type(to_replace) is not int:
        raise Exception("Config attribute \"genome_to_replace\" needs to be a number.\n")
    elif to_replace not in range(0, 3):
        raise Exception("Genome to replace must be 0, 1, or 2.\n")

    # Perform Guided Genome Halving on the given tetrad and outgroup genomes
    ggh: GroupGraph = GroupGraph(Genome.from_strings(tetrad), Genome.from_strings(outgroup), to_replace)
    ggh.get_result()

    bpg_distance: BPGDistance = BPGDistance(Genome.from_strings(tetrad), ggh.ancestor_AA)
    bpg_distance.calculate_distance()
    distance_1: int = bpg_distance.distance

    bpg_distance = BPGDistance(Genome.from_strings(outgroup), ggh.ancestor_A)
    bpg_distance.calculate_distance()
    distance_2: int = bpg_distance.distance

    total_distance: int = distance_1 + distance_2

    print("\nd(AA, tetra) = " + str(distance_1) + " | d(A,outgroup) = " + str(distance_2), " | total = " +
          str(total_distance) + "\n")

    print("\n-\nGenome ancestor_AA:\n")

    for i in range(len(ggh.ancestor_AA.chromosomes)):
        print("Chromosome " + str(i + 1))
        print(str(ggh.ancestor_AA.chromosomes[i]))

    print("\n-\nGenome ancestor_A:\n")

    for i in range(len(ggh.ancestor_A.chromosomes)):
        print("Chromosome " + str(i + 1))
        print(str(ggh.ancestor_A.chromosomes[i]))


def get_algorithm(alg: str):
    """
    Calls the appropriate function based on the user input.

    Parameters
    ----------
    alg : str
        User input, the name of the desired algorithm
    """
    if alg == "SmallPhylogeny":
        small_phylogeny()
    elif alg == "GenomeAliquoting":
        genome_aliquoting()
    elif alg == "DCJRearrangements":
        dcj_rearrangements()
    elif alg == "GenomeHalving":
        genome_halving()
    else:
        raise Exception("Algorithm not supported by this program. Supported algorithms:\n"
                        "- SmallPhylogeny\n"
                        "- GenomeAliquoting\n"
                        "- DCJRearrangements\n"
                        "- GenomeHalving\n\n")


def main():
    config_file: TextIO = open(CONFIG_DIR)
    config_data = yaml.safe_load(config_file)
    config_file.close()
    algorithm: str = config_data.get(CONFIG_ALGORITHM)

    get_algorithm(algorithm)


if __name__ == '__main__':
    cProfile.run('main()')
    # main()
