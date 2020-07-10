import sys
import yaml
import InputPreprocessing

from BPGDistance import BPGDistance
from DCJOperations import DCJOperations, OperationTypes
from GenomeInString import GenomeInString
from TreeStructure import TreeStructure
from SmallPhylogeny import SmallPhylogeny
from MedianIteration import MedianIteration
from PGMPathForAGenome import PGMPathForAGenome
from PGMPath import PGMPath

import networkx as nx
from networkx import Graph
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from io import StringIO
from typing import List, Dict, TextIO, Type, Optional
from enum import Enum
from copy import copy, deepcopy


"""
 Driver program for Pythgroups
 
 Command line usage: 
 $ python [directory of GenomeReconstruction.py] [algorithm]
 where [algorithm] is the desired genome reconstruction algorithm:
    - SmallPhylogeny
    - GenomeAliquoting
    - DCJRearrangements
 
 Based on the small phylogeny program developed by C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""

CONFIG_DIR = "config.yaml"
CONFIG_GENOME_FILE = "genome_file"
CONFIG_TREE_STRUCTURE = "tree_structure"


class GenomeNode:
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
    genomes: Dict[str, List[str]] = {}
    line: str = genome_file.readline()
    header: str = ""
    while len(line) != 0:

        # ">" indicates a header, sets that as the current genome to add chromosomes to
        if line.startswith(">"):
            header = line.replace(">", "").replace("\n", "")

        # Add chromosomes to the current genome until an empty line, another header, or the end of file is found
        else:
            chromosomes: List[str] = []
            while len(line) != 0 and not line.startswith(">") and line != "\n":
                chromosome: str = line.replace("$", "").replace("\n", "")
                chromosomes.append(chromosome)
                line = genome_file.readline()
            genomes[header] = chromosomes
        line = genome_file.readline()

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
    return Phylo.read(StringIO(config_data.get(CONFIG_TREE_STRUCTURE)), "newick")


def parse_medians(graph_nodes: List[GenomeNode]) -> List[GenomeNode]:
    """
    Get the median nodes from a list of GenomeNodes

    Parameters
    ----------
    graph_nodes : List[GenomeNode]
        List of GraphNodes to parse the medians from

    Returns
    -------
    Dict[Clade, List[Clade]]
        Key: median to reconstruct, value: list of 3 neighboring clades
    """
    medians: List[GenomeNode] = []
    for node in graph_nodes:
        if (len(node.neighbors)) == 3:
            medians.append(GenomeNode(node.name, node.neighbors, node.genome_id))
        # TODO: Exception handling for illegal tree structures (eg. node with > 3 neighbors)

    return medians


def genome_nodes_from_tree(tree: Tree) -> List[GenomeNode]:
    """
    Converts a tree into a networkx graph consisting of GenomeNode representations of each node

    Parameters
    ----------
    tree : Tree
        Tree object to convert into a graph

    Returns
    -------
    List[GenomeNode]
        A list of GenomeNodes representing each node, their neighbors, and their unique IDs
    """
    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()
    graph_indexes = {}
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

    graph_nodes = []
    for node in graph:
        genome_id = graph_indexes[node]
        neighbor_ids = []
        for neighbor in graph.neighbors(node):
            neighbor_ids.append(graph_indexes[neighbor])
        graph_nodes.append(GenomeNode(node, neighbor_ids, genome_id))

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
            count += len(chromosome.strip().split(" "))

        # Each genome must have the same number of genes (as compared to the first genome in genomes)
        if genome == list(genomes.keys())[0]:
            final_count = count
        elif final_count != count:
            print("Different numbers of genes in the genomes (check [{}] or the first genome for anomalies). \
                    Exiting.\n".format(genome))
            exit(1)

    return final_count


def small_phylogeny():
    """
    Reconstructs the ancestor(s) of known modern genomes given an unrooted binary phylogenetic tree
    using the Pathgroups algorithm.

    """
    tree: Tree = parse_tree(CONFIG_DIR)
    genomes: Dict[str, List[str]] = parse_genomes(CONFIG_DIR)
    genome_nodes: List[GenomeNode] = genome_nodes_from_tree(tree)
    median_nodes: List[GenomeNode] = parse_medians(genome_nodes)

    num_ancestor = len(tree.get_nonterminals())
    num_leaves = len(tree.get_terminals())
    num_genes = count_genes(genomes)
    all_genomes: List[GenomeInString] = []
    for chromosomes in genomes.values():
        all_genomes.insert(0, GenomeInString(chromosomes))

    ts: TreeStructure = TreeStructure(num_ancestor, num_leaves, num_genes, None, None, None, all_genomes)

    for median in median_nodes:
        print(str(median.genome_id) + str(median.neighbors))
        ts.set_tree_structure(median.genome_id, median.neighbors[0],
                              median.neighbors[1], median.neighbors[2])

    sp: SmallPhylogeny = SmallPhylogeny(ts)
    sp.get_result()

    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        print("before optimization, reconstructed ancestors")
        for j in range(0, len(ts.medians[i].medians)):
            print("chr {}\n {}".format(j, ts.medians[i].medians[j]))

    reconstructed_paths: List[PGMPathForAGenome] = []

    # Leaf paths added first
    if ts.number_of_leaves >= 0:
        reconstructed_paths = [ts.all_paths[i] for i in range(0, ts.number_of_leaves)]

    # Ancestor paths added second
    for i in range(ts.number_of_leaves, ts.number_of_leaves + ts.number_of_ancestors):
        reconstructed_paths.append(PGMPathForAGenome(ts.get_pgm_path(ts.medians[i - ts.number_of_leaves].medians, i)))

    relation: List[List[int]] = ts.get_relation()
    reconstructed_dist: int = 0
    before_optimization: str = ""
    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[PGMPath] = reconstructed_paths[i].paths
                p2: List[PGMPath] = reconstructed_paths[j].paths
                cur_dist: int = ts.medians[0].get_distance(p1, p2)
                reconstructed_dist += cur_dist
                before_optimization += "  d({r},{c})={d}".format(r=str(i), c=str(j), d=str(cur_dist))

    print("before optimization:\n" + before_optimization)
    print("total distance: " + str(reconstructed_dist))

    # Section 4: optimize the result

    mi: MedianIteration = MedianIteration(ts.number_of_leaves, ts.number_of_ancestors, ts.gene_number,
                                          ts.leaves, reconstructed_paths, ts.medians, ts.node_int, ts.node_string)
    mi.optimize_result(1, 50)
    after_optimization: str = ""
    optimized_dist: int = 0
    for i in range(0, len(relation) - 1):
        for j in range(i + 1, len(relation)):
            if relation[i][j] == 2 or relation[j][i] == 2:
                p1: List[PGMPath] = reconstructed_paths[i].paths
                p2: List[PGMPath] = reconstructed_paths[j].paths
                cur_dist: int = ts.medians[0].get_distance(p1, p2)
                optimized_dist += cur_dist
                after_optimization += "  d({r},{c})={d}".format(r=str(i), c=str(j), d=str(cur_dist))

    print("after optimization:\n" + after_optimization)
    print("total distance: " + str(optimized_dist))

    for i in range(0, len(ts.medians)):
        ts.medians[i].get_ancestors()
        print("reconstructed ancestors:")
        for j in range(0, len(ts.medians[i].medians)):
            print("chr {}\n {}".format(j, ts.medians[i].medians[j]))


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
    bpg_dist: BPGDistance = BPGDistance(genome1, genome2)
    bpg_dist.calculate_distance()
    cur_dist: int = bpg_dist.distance

    # Perform DCJ operations until distance == 0 or there are no more DCJ operations to perform (?)
    operation_types: List[int] = [OperationTypes.INVERSION, OperationTypes.TRANSLOCATION]
    dcj: DCJOperations = DCJOperations(genome1, genome2)
    dcj.initial_value()
    more: bool = True
    while cur_dist > 0 and more:
        more = False
        rearrange_state: List[GenomeInString] = dcj.get_result(1, 30, -2, operation_types, 10)
        if len(rearrange_state) > 0:
            more = True
            for genome in rearrange_state:
                print("*******")
                for i in range(0, len(genome.chromosomes)):
                    print("Chromosome " + str(i) + "\n" + genome.chromosomes[i])

            new_genome: GenomeInString = rearrange_state[len(rearrange_state) - 1]
            bpg_dist = BPGDistance(new_genome.chromosomes, genome2)
            bpg_dist.calculate_distance()
            cur_dist = bpg_dist.distance
            print("a run, steps: " + str(len(rearrange_state)) + ", cur_dist: " + str(cur_dist))
            dcj = DCJOperations(new_genome.chromosomes, genome2)
            dcj.initial_value()


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
    else:
        print("Algorithm not supported by this program. Supported algorithms:\n" +
              "- SmallPhylogeny\n" +
              "- GenomeAliquoting\n" +
              "- DCJRearrangements\n\n" +
              "Exiting.")
        exit(1)


def main():
    # algorithm = sys.argv[1]  # The first argument when calling the program
    algorithm = "GenomeAliquoting"
    get_algorithm(algorithm)


if __name__ == '__main__':
    main()
