from copy import deepcopy
from sys import exit
from typing import Optional, List

from BPGPath import BPGPath
from Genome import Genome

"""                              
 Calculate the distance between two genomes using          
 a breakpoint graph and counting the number of cycles.     
                                                           
 Based on BPGDistance.java from C.Zheng & D.Sankoff (2011) 
                                                           
 Author: Holger Jensen, Oskar Jensen                                     
"""


def insert_character(s: str, i: int, c: str) -> str:
    """
    Utility function for inserting a character into a string.

    Parameters
    ----------
    s
        String to be modified
    i
        Index to change
    c
        Character to be inserted

    Returns
    -------
    str
        New string with the inserted character
    """
    return s[:i] + c + s[(i + 1):]


def is_valid_cycle(edge1: int, edge2: int) -> bool:
    """
    Checks if the given path is a valid cycle.

    Parameters
    ----------
    edge1
        Head or Tail of a path
    edge2
        Tail or Head of a path

    Returns
    -------
    bool
        True if edge1 and edge2 are negative AND edge1 or edge2 are even
    """
    if edge1 >= 0 or edge2 >= 0:
        return False

    if int(edge1 / 2) * 2 != edge1:
        return True

    return int(edge2 / 2) * 2 != edge2


def get_valid_path(start_from: int, paths: List[BPGPath]) -> Optional[BPGPath]:
    """
    Gets the first non-None value from a list of paths

    Parameters
    ----------
    start_from
        Which index to begin search
    paths
        List of paths

    Returns
    -------
    BPGPath
        Found BPGPath
    """
    for i in range(start_from, len(paths)):
        if paths[i] is not None:
            return BPGPath.from_path(paths[i])

    return None


class BPGDistance:
    """
    Attributes
    ----------
    genome_1 : [str]
        First genome for distance calculation
    genome_2 : [str]
        Genome to compare to the first
    node_ints : [int]
        Integer values of genes
    node_strings_1 : [Optional[str]]
        String representation for genes in genome 1
    node_strings_2 : [Optional[str]]
        String representation for genes in genome 2
    genome_paths_1 : [BPGPath]
        Paths for combining into cycles in genome 1
    genome_paths_2 : [BPGPath]
        Paths for combining into cycles in genome 2
    distance : int
        Integer distance between genome 1 and genome 2
    """

    def __init__(self, genome_1: Genome, genome_2: Genome):
        """
        Constructor

        Parameters
        ----------
        genome_1
            First genome for distance calculation
        genome_2
            Genome to test distance from the first one
        """
        gene_count_1: int = 0
        gene_count_2: int = 0
        self.genome_1: Genome = Genome(list())  # Genome 1 as String
        self.genome_2: Genome = Genome(list())  # Genome 2 as String

        for chromosome in genome_1.chromosomes:
            if len(chromosome.genes) != 0:
                gene_count_1 += len(chromosome.genes)
                self.genome_1.add_chromosome(chromosome)

        for chromosome in genome_2.chromosomes:
            if len(chromosome.genes) != 0:
                gene_count_2 += len(chromosome.genes)
                self.genome_2.add_chromosome(chromosome)

        if gene_count_1 == gene_count_2:  # Algorithm requires genomes are equal length
            self.gene_count: int = gene_count_1  # Number of genes
        else:
            print("Different numbers of genes in both genomes. Exiting.\n")
            exit(1)

        self.node_ints: List[int] = [int() for _ in range((self.gene_count * 2))]
        self.node_strings_1: List[Optional[str]] = [None for _ in range((self.gene_count * 2))]
        self.node_strings_2: List[Optional[str]] = [None for _ in range((self.gene_count * 2))]
        self.genome_paths_1: Optional[List[Optional[BPGPath]]] = None
        self.genome_paths_2: Optional[List[Optional[BPGPath]]] = None
        self.distance: int = int()

    def graph_init(self):
        """
        Initialize BPG graph to the state which contains the set of gene/telomere paths
        """
        index1: int = 0

        for chromosome in self.genome_1.chromosomes:
            for gene in chromosome.genes:
                first_character: str = deepcopy(gene.name[0])  # Sign indicating gene is head-tail or tail-head
                node1: str
                node2: str

                if first_character == "-":
                    node1 = insert_character(gene.name[1:], len(gene.name[1:]), "h")
                    node2 = insert_character(gene.name[1:], len(gene.name[1:]), "t")
                    self.node_ints[index1] = index1 + 1
                    self.node_strings_1[index1] = node2

                    index1 += 1
                    self.node_ints[index1] = index1 + 1
                    self.node_strings_1[index1] = node1
                else:
                    node1 = insert_character(gene.name, len(gene.name), "t")
                    node2 = insert_character(gene.name, len(gene.name), "h")
                    self.node_ints[index1] = index1 + 1
                    self.node_strings_1[index1] = node1

                    index1 += 1
                    self.node_ints[index1] = index1 + 1
                    self.node_strings_1[index1] = node2

                index1 += 1

        self.node_strings_2 = deepcopy(self.node_strings_1)

        self.genome_paths_1 = self.get_paths(self.genome_1, 1)
        self.genome_paths_2 = self.get_paths(self.genome_2, 2)

    def get_paths(self, genome: Genome, node: int) -> List[BPGPath]:
        """
        Gets the gene/telomere paths from the given genome and node

        Parameters
        ----------
        genome
            Genome which the path is constructed from
        node
            Current node to use for the path

        Returns
        -------
        [BPGPath]
            A list of BPGPaths
        """
        null_node: int = -node
        path1: List[Optional[BPGPath]] = [None for _ in range(((self.gene_count * 2) + 1))]

        for chromosome in genome.chromosomes:
            pre_node: int = 0

            for i in range(len(chromosome.genes)):
                first_character: str = chromosome.genes[i].name[0]  # Sign indicating gene is head-tail or tail-head
                node1: str
                node2: str

                if first_character == "-":
                    node1 = insert_character(chromosome.genes[i].name[1:], len(chromosome.genes[i].name[1:]), "h")
                    node2 = insert_character(chromosome.genes[i].name[1:], len(chromosome.genes[i].name[1:]), "t")
                else:
                    node1 = insert_character(chromosome.genes[i].name, len(chromosome.genes[i].name), "t")
                    node2 = insert_character(chromosome.genes[i].name, len(chromosome.genes[i].name), "h")

                node1_int: int = self.get_node_int(node1, node)
                node2_int: int = self.get_node_int(node2, node)

                if i == 0:
                    path1[node1_int] = BPGPath(node1_int, null_node, node, node)
                    pre_node: int = node2_int
                    null_node -= 2

                if i != 0 and i != len(chromosome.genes) - 1:
                    path1[node1_int] = BPGPath(node1_int, pre_node, node, node)
                    path1[pre_node] = BPGPath(pre_node, node1_int, node, node)
                    pre_node: int = node2_int

                if i == len(chromosome.genes) - 1:
                    if len(chromosome.genes) != 1:
                        path1[pre_node] = BPGPath(pre_node, node1_int, node, node)
                        path1[node1_int] = BPGPath(node1_int, pre_node, node, node)

                    path1[node2_int] = BPGPath(node2_int, null_node, node, node)
                    null_node -= 2

        return path1

    def get_node_int(self, gene_str: str, gene_int: int) -> int:
        """
        Gets the node's integer value from the full string representation

        Parameters
        ----------
        gene_str
            Gene as a string
        gene_int
            Gene as an integer

        Returns
        -------
        int
            Integer of node
        """
        if gene_int == 1:
            for i in range(len(self.node_strings_1)):
                if self.node_strings_1[i] is not None and self.node_strings_1[i] == gene_str:
                    self.node_strings_1[i] = None
                    return i + 1

        if gene_int == 2:
            for i in range(len(self.node_strings_2)):
                if self.node_strings_2[i] is not None and self.node_strings_2[i] == gene_str:
                    self.node_strings_2[i] = None
                    return i + 1

        return 0

    def calculate_distance(self):
        """
        Produces cycles from the gene/telomere paths and counts them to calculate genome distance
        """
        self.graph_init()

        cycle_number: int = 0
        good_path_number: int = 0

        # Get the first path out of a list of paths with None values filtered out
        ancestor_path1: Optional[BPGPath] = get_valid_path(0, self.genome_paths_1)

        while ancestor_path1 is not None:
            node1: int = ancestor_path1.head
            node2: int = ancestor_path1.tail
            start: int = node1

            if node1 > 0:
                self.genome_paths_1[node1] = None

            if node2 > 0:
                self.genome_paths_1[node2] = None

            node_big: int = node1
            node_small: int = node2

            if node_big < node_small:
                node_big = node2
                node_small = node1

            more: bool = True

            while more:
                if ancestor_path1.head < 0 and ancestor_path1.tail < 0:
                    if is_valid_cycle(ancestor_path1.head, ancestor_path1.tail):
                        good_path_number += 1

                    break

                l_path: BPGPath = self.genome_paths_2[node_big]
                l_node1: int = l_path.head
                l_node2: int = l_path.tail

                if l_node2 > 0 and l_node2 == node_small:
                    cycle_number += 1
                    more = False
                    self.genome_paths_2[l_node1] = None
                    self.genome_paths_2[l_node2] = None
                elif l_node2 > 0 and l_node2 != node_small:
                    a_path2: BPGPath = self.genome_paths_1[l_node2]
                    another_node: int = a_path2.tail

                    ancestor_path1 = BPGPath(another_node, node_small, 1, 1)
                    node1 = ancestor_path1.head
                    node2 = ancestor_path1.tail

                    node_big = node1
                    node_small = node2

                    if node_big < node_small:
                        node_big = node2
                        node_small = node1

                    other_path1: int = self.genome_paths_1[l_node2].tail

                    if other_path1 > 0:
                        self.genome_paths_1[other_path1] = None

                    self.genome_paths_2[l_node1] = None
                    self.genome_paths_2[l_node2] = None
                    self.genome_paths_1[l_node2] = None
                elif l_node2 < 0 and node_small < 0:
                    if is_valid_cycle(l_node2, node_small):
                        good_path_number += 1

                    more = False
                    self.genome_paths_2[l_node1] = None
                elif l_node2 < 0 < node_small:
                    ancestor_path1 = BPGPath(node_small, l_node2, 1, 2)
                    node1 = ancestor_path1.head
                    node2 = ancestor_path1.tail

                    node_big = node1
                    node_small = node2

                    if node_big < node_small:
                        node_big = node2
                        node_small = node1

                    self.genome_paths_2[l_node1] = None

            ancestor_path1 = get_valid_path(start, self.genome_paths_1)

        self.distance = self.gene_count + len(self.genome_1.chromosomes) - cycle_number - good_path_number
