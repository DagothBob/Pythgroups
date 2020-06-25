from copy import deepcopy
from sys import exit
from typing import Optional, List

from BPGPath import BPGPath

"""                              
 Calculate the distance between two genomes using          
 a breakpoint graph and counting the number of cycles.     
                                                           
 Based on BPGDistance.java from C.Zheng & D.Sankoff (2011) 
                                                           
 Author: Holger Jensen                                     
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


def split_at_whitespace(strings: str) -> List[str]:
    """
    Strips, then splits, a string, then strips the substrings again

    Parameters
    ----------
    strings
        List of strings to operate on

    Returns
    -------
    [str]
        Set of cleaned-up strings
    """
    result: List[str] = []

    for string in strings.strip().split(" "):
        if string.strip() != "":
            result.append(string.strip())

    return result


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
    genome1 : [str]
        First genome for distance calculation
    genome2 : [str]
        Genome to compare to the first
    node_ints : [int]
        Integer values of genes
    node_strs1 : [Optional[str]]
        String representation for genes in genome 1
    node_strs2 : [Optional[str]]
        String representation for genes in genome 2
    genome_paths1 : [BPGPath]
        Paths for combining into cycles in genome 1
    genome_paths2 : [BPGPath]
        Paths for combining into cycles in genome 2
    distance : int
        Integer distance between genome 1 and genome 2
    """

    def __init__(self, genome1: List[str], genome2: List[str]):
        """
        Constructor

        Parameters
        ----------
        genome1
            First genome for distance calculation
        genome2
            Genome to test distance from the first one
        """
        index1: int = 0
        index2: int = 0
        num_genes1: int = 0
        num_genes2: int = 0
        self.genome1: List[str] = [str() for _ in range(len(genome1))]  # Genome 1 as String
        self.genome2: List[str] = [str() for _ in range(len(genome2))]  # Genome 2 as String

        for chromosome in genome1:
            genes: List[str] = split_at_whitespace(chromosome)

            if len(genes) != 0:
                num_genes1 += len(genes)
                self.genome1[index1] = chromosome
                index1 += 1

        for chromosome in genome2:
            genes: List[str] = split_at_whitespace(chromosome)

            if len(genes) != 0:
                num_genes2 += len(genes)
                self.genome2[index2] = chromosome
                index2 += 1

        if num_genes1 == num_genes2:  # Algorithm requires genomes are equal length
            self.gene_number: int = num_genes1  # Number of genes
        else:
            print("Different numbers of genes in both genomes. Exiting.\n")
            exit(1)

        self.node_ints: List[int] = [int() for _ in range((self.gene_number * 2))]
        self.node_strs1: List[Optional[str]] = [None for _ in range((self.gene_number * 2))]
        self.node_strs2: List[Optional[str]] = [None for _ in range((self.gene_number * 2))]
        self.genome_paths1: List[Optional[BPGPath]] = list()
        self.genome_paths2: List[Optional[BPGPath]] = list()
        self.distance: int = int()

    def graph_init(self):
        """
        Initialize BPG graph to the state which contains the set of gene/telomere paths
        """
        index1: int = 0

        for chromosome in self.genome1:
            genes: List[str] = split_at_whitespace(chromosome)

            for gene in genes:
                first_character: str = deepcopy(gene[0])  # Sign indicating gene is head-tail or tail-head
                node1: str
                node2: str

                if first_character == "-":
                    node1 = insert_character(gene[1:], len(gene[1:]), "h")
                    node2 = insert_character(gene[1:], len(gene[1:]), "t")
                    self.node_ints[index1] = index1 + 1
                    self.node_strs1[index1] = node2

                    index1 += 1
                    self.node_ints[index1] = index1 + 1
                    self.node_strs1[index1] = node1
                else:
                    node1 = insert_character(gene, len(gene), "t")
                    node2 = insert_character(gene, len(gene), "h")
                    self.node_ints[index1] = index1 + 1
                    self.node_strs1[index1] = node1

                    index1 += 1
                    self.node_ints[index1] = index1 + 1
                    self.node_strs1[index1] = node2

                index1 += 1

        self.node_strs2 = deepcopy(self.node_strs1)

        self.genome_paths1 = self.get_paths(self.genome1, 1)
        self.genome_paths2 = self.get_paths(self.genome2, 2)

    def get_paths(self, genome: List[str], node: int) -> List[BPGPath]:
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
        path1: List[Optional[BPGPath]] = [None for _ in range(((self.gene_number * 2) + 1))]

        for chromosome in genome:
            genes: List[str] = split_at_whitespace(chromosome)
            pre_node: int = 0

            for i in range(len(genes)):
                first_character: str = genes[i][0]  # Sign indicating gene is head-tail or tail-head
                node1: str
                node2: str

                if first_character == "-":
                    node1 = insert_character(genes[i][1:], len(genes[i][1:]), "h")
                    node2 = insert_character(genes[i][1:], len(genes[i][1:]), "t")
                else:
                    node1 = insert_character(genes[i], len(genes[i]), "t")
                    node2 = insert_character(genes[i], len(genes[i]), "h")

                node1_int: int = self.get_node_int(node1, node)
                node2_int: int = self.get_node_int(node2, node)

                if i == 0:
                    path1[node1_int] = BPGPath(node1_int, null_node, node, node)
                    pre_node: int = node2_int
                    null_node -= 2

                if i != 0 and i != len(genes) - 1:
                    path1[node1_int] = BPGPath(node1_int, pre_node, node, node)
                    path1[pre_node] = BPGPath(pre_node, node1_int, node, node)
                    pre_node: int = node2_int

                if i == len(genes) - 1:
                    if len(genes) != 1:
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
            for i in range(len(self.node_strs1)):
                if self.node_strs1[i] is not None and self.node_strs1[i] == gene_str:
                    self.node_strs1[i] = None
                    return i + 1

        if gene_int == 2:
            for i in range(len(self.node_strs2)):
                if self.node_strs2[i] is not None and self.node_strs2[i] == gene_str:
                    self.node_strs2[i] = None
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
        ancestor_path1: Optional[BPGPath] = get_valid_path(0, self.genome_paths1)

        while ancestor_path1 is not None:
            node1: int = ancestor_path1.head
            node2: int = ancestor_path1.tail
            start: int = node1

            if node1 > 0:
                self.genome_paths1[node1] = None

            if node2 > 0:
                self.genome_paths1[node2] = None

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

                l_path: BPGPath = self.genome_paths2[node_big]
                l_node1: int = l_path.head
                l_node2: int = l_path.tail

                if l_node2 > 0 and l_node2 == node_small:
                    cycle_number += 1
                    more = False
                    self.genome_paths2[l_node1] = None
                    self.genome_paths2[l_node2] = None
                elif l_node2 > 0 and l_node2 != node_small:
                    a_path2: BPGPath = self.genome_paths1[l_node2]
                    another_node: int = a_path2.tail

                    ancestor_path1 = BPGPath(another_node, node_small, 1, 1)
                    node1 = ancestor_path1.head
                    node2 = ancestor_path1.tail

                    node_big = node1
                    node_small = node2

                    if node_big < node_small:
                        node_big = node2
                        node_small = node1

                    other_path1: int = self.genome_paths1[l_node2].tail

                    if other_path1 > 0:
                        self.genome_paths1[other_path1] = None

                    self.genome_paths2[l_node1] = None
                    self.genome_paths2[l_node2] = None
                    self.genome_paths1[l_node2] = None
                elif l_node2 < 0 and node_small < 0:
                    if is_valid_cycle(l_node2, node_small):
                        good_path_number += 1

                    more = False
                    self.genome_paths2[l_node1] = None
                elif l_node2 < 0 < node_small:
                    ancestor_path1 = BPGPath(node_small, l_node2, 1, 2)
                    node1 = ancestor_path1.head
                    node2 = ancestor_path1.tail

                    node_big = node1
                    node_small = node2

                    if node_big < node_small:
                        node_big = node2
                        node_small = node1

                    self.genome_paths2[l_node1] = None

            ancestor_path1 = get_valid_path(start, self.genome_paths1)

        self.distance = self.gene_number + len(self.genome1) - cycle_number - good_path_number
