from io import StringIO
from typing import List, Dict, ValuesView, Iterator, TextIO, Any
import time
import threading

import yaml
from Bio import Phylo
from numpy.core.multiarray import ndarray

import InputPreprocessing
from BPGDistance import BPGDistance
from DCJOperation import OperationTypes
from DCJRearrangement import DCJRearrangement
from Genome import Genome, split_at_whitespace
from GroupGraph import GroupGraph
from MedianIteration import MedianIteration
from PGMPathForAGenome import PGMPathForAGenome
from SmallPhylogeny import SmallPhylogeny
from TreeStructure import TreeStructure
import NetworkxNode

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

# General
CONFIG_DIR = "config.yaml"
CONFIG_ALGORITHM = "algorithm"
CONFIG_GENOME_FILE = "genome_file"

# SmallPhylogeny
CONFIG_TREE_STRUCTURE = "tree_structure"
CONFIG_SHOW_DCJR = "show_DCJR"

# DCJRearrangements
CONFIG_OPERATIONS = "operations"
CONFIG_VERBOSE_OUTPUT = "verbose_output"
CONFIG_MINIMUM_CHROMOSOME = "minimum_chromosome"
CONFIG_MAXIMUM_CHROMOSOME = "maximum_chromosome"
CONFIG_WHICH_CHROMOSOME = "which_chromosome"
CONFIG_NUMBER_OPERATIONS = "number_of_operations"

# GGH
CONFIG_GENOME_REPLACE = "genome_to_replace"


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


def config_get(parameter: str) -> Any:
    """ Shorthand for getting a parameter from the config file

    Parameters
    ----------
    parameter :
        Parameter to get the data from in the config file

    Returns
    -------
    Any
        The data from the specified parameter
    """
    with open(CONFIG_DIR, "r") as config_file:
        config_contents = yaml.safe_load(config_file)
        parameter_data = config_contents.get(parameter)

    return parameter_data


def parse_genomes() -> Dict[str, List[str]]:
    """
    Parses the genome data from a file into a dictionary containing each genome and their chromosomes

    Returns
    -------
    Dict[str, List[str]]
        A dictionary containing each genome's chromosomes,
        where keys are genome headers and values are lists of chromosomes
    """
    genome_file = open(config_get(CONFIG_GENOME_FILE))
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
                chromosome: str = str().join([string.strip() + " "
                                              for string in clean_input.strip().split(" ")
                                              if string.strip() != ""])

                chromosomes.append(chromosome)
                line = genome_file.readline()
            genomes[header] = chromosomes
        line = genome_file.readline()

    return genomes


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
    # # # # # # # # # # # # # # #
    # Step 1: Parse input data  #
    # # # # # # # # # # # # # # #

    tree = Phylo.read(StringIO(config_get(CONFIG_TREE_STRUCTURE)), "newick")
    genomes: Dict[str, List[str]] = parse_genomes()
    genome_nodes: List[NetworkxNode] = NetworkxNode.genome_nodes_from_tree(tree, list(genomes.keys()))
    median_nodes: List[NetworkxNode] = NetworkxNode.parse_medians(genome_nodes)

    # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Step 2: Set up tree structure for the algorithm #
    # # # # # # # # # # # # # # # # # # # # # # # # # #

    num_ancestor = len([node for node in tree.get_nonterminals() if node.name is not None])
    num_leaves = len([node for node in tree.get_terminals() if node.name is not None])  # don't count unnamed nodes
    num_genes = count_genes(genomes)
    all_genomes: List[Genome] = list()

    for chromosomes in genomes.values():
        all_genomes.append(Genome.from_strings(chromosomes))

    ts: TreeStructure = TreeStructure(num_ancestor, num_leaves, num_genes, None, None, None, all_genomes)

    for median in median_nodes:
        ts.set_tree_structure(median.genome_id, median.neighbors[0],
                              median.neighbors[1], median.neighbors[2])

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Step 3: Ancestor reconstruction, before optimization  #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    sp: SmallPhylogeny = SmallPhylogeny(ts, True)

    # run sp.get_result in a separate thread from the progress spinner
    task = threading.Thread(target=sp.get_result)
    task.start()

    # progress spinner for before optimization step
    spin_thread = SpinnerThread()
    spin_thread.start()

    task.join()
    spin_thread.stop()
    spin_thread.join()

    print("\nReconstructed ancestors (pre-optimization):")
    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        median_node: NetworkxNode = NetworkxNode.get_node(genome_nodes, i + num_leaves)
        print(">{}".format(median_node.name))
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

    relation: ndarray = ts.get_relation()
    reconstructed_dist: int = int()
    before_optimization: str = str()

    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[Dict[str, int]] = reconstructed_paths[i].paths
                p2: List[Dict[str, int]] = reconstructed_paths[j].paths
                cur_dist: int = ts.medians[0].get_distance(p1, p2)
                reconstructed_dist += cur_dist

                node1: NetworkxNode = NetworkxNode.get_node(genome_nodes, i)
                node2: NetworkxNode = NetworkxNode.get_node(genome_nodes, j)
                before_optimization += "  d({r},{c})={d}".format(r=str(node1.name), c=str(node2.name), d=str(cur_dist))

    print("Calculated distances:\n" + before_optimization)
    print("Total distance: " + str(reconstructed_dist))
    print("")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Section 4: Ancestor reconstruction, optimization step #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    mi: MedianIteration = MedianIteration(ts.number_of_leaves, ts.number_of_ancestors, ts.gene_number,
                                          ts.leaves, reconstructed_paths, ts.medians, ts.node_int, ts.node_string)
    mi.optimize_result(1, 50)
    after_optimization: str = str()
    optimized_dist: int = int()

    print("Reconstructed ancestors (post-optimization):")
    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        median_node: NetworkxNode = NetworkxNode.get_node(genome_nodes, i + num_leaves)
        print(">{}".format(median_node.name))
        for j in range(0, len(ts.medians[i].median)):
            print("chr {}\n {}".format(j, ts.medians[i].median[j]))

    distances: Dict[(str, str), str] = {}

    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[Dict[str, int]] = mi.all_paths[i].paths
                p2: List[Dict[str, int]] = mi.all_paths[j].paths
                cur_dist: int = mi.medians[0].get_distance(p1, p2)
                optimized_dist += cur_dist

                node1: NetworkxNode = NetworkxNode.get_node(genome_nodes, i)
                node2: NetworkxNode = NetworkxNode.get_node(genome_nodes, j)
                after_optimization += "  d({r},{c})={d}".format(r=str(node1.name), c=str(node2.name), d=str(cur_dist))
                distances[(node1, node2)] = str(cur_dist)

    print("Calculated distances:\n" + after_optimization)
    print("Total distance: " + str(optimized_dist))

    NetworkxNode.print_tree_structure(genome_nodes, distances)


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
    verbose_output: bool = bool(config_get("verbose_output"))

    # Create the dictionary of genomes from the input file
    genomes: Dict[str, List[str]] = parse_genomes()

    # Get the first 2 genomes from the input file
    values_view = genomes.values()
    value_iterator = iter(values_view)
    genome1: List[str] = next(value_iterator)
    genome2: List[str] = next(value_iterator)

    keys_view = genomes.keys()
    key_iterator = iter(keys_view)
    genome1_name: str = next(key_iterator)
    genome2_name: str = next(key_iterator)

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
    operation_counts: Dict[int, int] = {OperationTypes.INVERSION: 0,
                                        OperationTypes.TRANSLOCATION: 0,
                                        OperationTypes.FISSION: 0,
                                        OperationTypes.FUSION: 0}
    dcj: DCJRearrangement = DCJRearrangement(Genome.from_strings(genome1), Genome.from_strings(genome2), verbose_output)
    dcj.initial_value()
    more: bool = True

    while cur_dist > 0 and more:
        more = False
        rearrange_state: List[Genome] = dcj.get_result(minimum_chromosome,
                                                       maximum_chromosome,
                                                       which_chromosome,
                                                       operation_types,
                                                       number_operations)

        if len(rearrange_state) != 0:
            more = True
            if verbose_output:
                for genome in rearrange_state:
                    print("*******")
                    for i in range(len(genome.chromosomes)):
                        print("Chromosome " + str(i) + "\n" + str(genome.chromosomes[i]))

            new_genome: Genome = rearrange_state[len(rearrange_state) - 1]
            bpg_dist = BPGDistance(Genome(new_genome.chromosomes), Genome.from_strings(genome2))
            bpg_dist.calculate_distance()
            cur_dist = bpg_dist.distance
            if verbose_output:
                print("a run, steps: " + str(len(rearrange_state)) + ", cur_dist: " + str(cur_dist))
            else:
                print("\rCalculating DCJRearrangements, ""current distance between {g1} and {g2}: {d}".format(
                        g1=genome1_name, g2=genome2_name, d=cur_dist), end="")
            
            for key, value in operation_counts.items():
                operation_counts[key] += dcj.operation_counts[key]
            dcj = DCJRearrangement(Genome(new_genome.chromosomes), Genome.from_strings(genome2), verbose_output)
            dcj.initial_value()

    print("\rTotal number of DCJ rearrangement events between {g1} and {g2} is {td}, including:\n"
          "x{inv} inversions, x{tra} translocations, x{fis} fissions, x{fus} fusions.\n"
          "Remaining BPG distance between genomes after operations: {d}".format(
            g1=genome1_name, g2=genome2_name, td=sum(operation_counts.values()),
            inv=operation_counts[OperationTypes.INVERSION], tra=operation_counts[OperationTypes.TRANSLOCATION],
            fis=operation_counts[OperationTypes.FISSION], fus=operation_counts[OperationTypes.FUSION], d=cur_dist))


def genome_halving():
    # Create the dictionary of genomes from the input file
    genomes: Dict[str, List[str]] = parse_genomes()

    # Get the first 2 genomes from the input file
    values_view: ValuesView[List[str]] = genomes.values()
    value_iterator: Iterator[List[str]] = iter(values_view)
    tetrad: List[str] = next(value_iterator)
    outgroup: List[str] = next(value_iterator)

    # Get GenomeHalving configuration options
    to_replace: int = config_get(CONFIG_GENOME_REPLACE)
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
    algorithm: str = config_get(CONFIG_ALGORITHM)
    get_algorithm(algorithm)


if __name__ == '__main__':
    # cProfile.run('main()')
    main()
