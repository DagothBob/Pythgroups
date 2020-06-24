from io import StringIO
from typing import List, Dict, Optional

import yaml
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from networkx import Graph

from BPGDistance import BPGDistance
from DCJOperations import DCJOperations, OperationTypes
from GenomeInString import GenomeInString

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


def parse_medians(graph: Graph) -> Dict[Clade, List[Clade]]:
    """
    Given a tree converted into a graph, parse the medians (ie nodes with 3 neighbours)

    Parameters
    ----------
    graph : Graph
        Graph to parse medians from

    Returns
    -------
    Dict[Clade, List[Clade]]
        Key: median to reconstruct, value: list of 3 neighbouring clades
    """
    medians = {}
    for node in graph.nodes:
        neighbours = [n for n in graph.neighbors(node)]
        if (len(neighbours)) == 3:
            for neighbour in neighbours:
                if node not in medians.keys():
                    medians[node] = []
                medians[node].append(neighbour)
        # TODO: Exception handling for illegal tree structures (eg. node with > 3 neighbours)

    return medians


def graph_from_phylo_tree(tree: Tree) -> Graph:
    """
    Converts a tree into a networkx graph by getting all its nonterminals
    and adding edges between them and their subclades

    Parameters
    ----------
    tree : Tree
        Tree object to convert into a graph

    Returns
    -------
    Graph
        A graph object converted from the tree
    """
    graph: Graph = Graph()
    nonterminals: List[Clade] = tree.get_nonterminals()
    for nt in nonterminals:
        if nt.name is not None and len(nt.clades) > 0:
            for clade in nt.clades:
                graph.add_edge(nt, clade)

    return graph


def small_phylogeny() -> str:
    """
    Reconstructs the ancestor(s) of known modern genomes given an unrooted binary phylogenetic tree
    using the Pathgroups algorithm.

    Returns
    -------

    """
    tree: Tree = parse_tree(CONFIG_DIR)
    genomes: Dict[str, List[str]] = parse_genomes(CONFIG_DIR)
    graph: Graph = graph_from_phylo_tree(tree)
    medians: Dict[Clade, List[Clade]] = parse_medians(graph)
    for m, n in medians.items():
        print("median: " + str(m) + ", neighbours: " + str(n))

    Phylo.draw_ascii(tree)

    return "Tree: " + str(tree)


def genome_aliquoting() -> str:
    """
    See the 2010 paper, section 2.5

    Returns
    -------

    """

    return "Hello yes this is genome aliquoting"


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
                index: int = 0
                print("*******")
                for chromosome in genome.chromosomes:
                    print("Chromosome " + str(index) + "\n" + chromosome)
                    index += 1

            new_genome: GenomeInString = rearrange_state[len(rearrange_state) - 1]
            bpg_dist = BPGDistance(new_genome.chromosomes, genome2)
            bpg_dist.calculate_distance()
            cur_dist = bpg_dist.distance
            print("a run, steps: " + str(len(rearrange_state)) + ", cur_dist: " + str(cur_dist))
            dcj = DCJOperations(new_genome.chromosomes, genome2)
            dcj.initial_value()


def get_algorithm(alg: str) -> str:
    """
    Calls the appropriate function based on the user input.

    Parameters
    ----------
    alg : str
        User input, the name of the desired algorithm

    Returns
    -------
    function
        The appropriate function based on the user input, or an error if the algorithm is invalid
    """
    # Essentially a switch statement
    return {
        "SmallPhylogeny": small_phylogeny(),
        "GenomeAliquoting": genome_aliquoting(),
        "DCJRearrangements": dcj_rearrangements()
    }.get(alg, "Algorithm doesn't exist")


def main():
    # algorithm = sys.argv[1]  # The first argument when calling the program
    output = get_algorithm("DCJRearrangements")
    print(output)


if __name__ == '__main__':
    main()
