from PGMPath import PGMPath
from PGMFragment import PGMFragment
from ChoiceStructure import ChoiceStructure
import math
from typing import List


def get_gene_next_node(end1: int) -> int:
    """
    If end1 is even then it's greater than end2, else it's less than end2 (I think)
    Recall the pattern described in __init__ for self.fragments[]: [(1, 2) (2, 1), (3, 4), (4, 3), ... ]

    Parameters
    ----------
    end1
        end1 from PGMFragment

    Returns
    -------
    int
        end1 - 1 if even, else end + 1
        Some possible hints on why it returns what it returns:
        (2, 1) is just (1, 2) flipped, so end1 - 1 from (2, 1) == end1 from (1, 2)
        In the case of (1, 2), end1 + 1 from (1, 2) == end2 from (1, 2).
    """

    if end1 % 2 == 0:
        return end1 - 1
    return end1 + 1


def find_gray_edge_node(end1: int, gray_edge: [PGMPath]) -> int or None:
    """


    Parameters
    ----------
    end1
        end1 from PGMFragment
    gray_edge
        The gray edge to search for

    Returns
    -------
    int or None
        If end1 found in gray_edge, return the opposite of what it matched with (head/tail), else return None
    """

    for path in gray_edge:
        if path is not None:
            if end1 == path.head:
                return path.tail
            if end1 == path.tail:
                return path.head

    return None


class MedianData:
    """
    Attributes
    ----------
    self.three_genome_paths : [[PGMPath]]
        All paths for genomes 1-3 in the tree structure
    self.three_genomes : [int]
        The identifiers for the three leaves (genomes) in the tree structure
    self.gene_num : int
        Total number of genes in each genome
    self.which_genome : int
        Which genome in the tree structure these medians exist in
    self.node_strings : [str]
        String representation of each node
    self.gray_edge : [PGMPath]
        List of PGMPaths representing a gray edge
    self.gray_edge_index : int
        Index of the gray edge
    self.fragments : [PGMFragment]
        List of each fragment
    self.choice_structures : [ChoiceStructure]
        List of each choice structure
    self.medians : [str]
        List of each median
    """

    def __init__(self, paths1: [PGMPath], paths2: [PGMPath], paths3: [PGMPath], gene_num: int, which_genome: int,
                 genome1: int, genome2: int, genome3: int, node_str: [str]):

        self.three_genome_paths = [paths1, paths2, paths3]
        self.three_genomes = [genome1, genome2, genome3]
        self.gene_num = gene_num
        self.which_genome = which_genome
        self.node_strings = node_str

        self.gray_edge = List[PGMPath]
        self.gray_edge_index = 0
        self.fragments = List[PGMFragment]
        self.choice_structures = [ChoiceStructure]
        self.medians = List[str]

        # Each gene in each genome has 1 fragment initially (see 2010 paper, section 3.1.1)
        # [(1, 2), (2, 1), (3, 4), (4, 3), ..., (1999, 2000), (2000, 1999)] for gene_num == 1000
        fragment_size = self.gene_num * 2 + 1
        for i in range(0, math.trunc(fragment_size / 2)):
            self.fragments[2 * i + 1] = PGMFragment(2 * i + 1, 2 * i + 2)
            self.fragments[2 * i + 2] = PGMFragment(2 * i + 2, 2 * i + 1)

        cs_index = 0
        for i in range(1, math.ceil(self.gene_num * 2 + 1)):
            cs = ChoiceStructure()
            cs.index_from = i
            cs.genome_1_path = self.three_genome_paths[0][i]
            cs.genome_2_path = self.three_genome_paths[1][i]
            cs.genome_3_path = self.three_genome_paths[2][i]
            cs.for_which_genome = self.which_genome
            cs.priority = 200
            cs.position = -1
            cs.gray_edge = None

            self.choice_structures[i] = cs
            cs_index += 1

    def get_ancestors(self):
        """
        Fills medians[] with all the genome's medians (ancestors),
        which are formatted for the purposes of printing to screen.
        """

        for i in range(0, len(self.fragments)):
            if self.fragments[i] is not None:
                other_frag = self.fragments[i].end2
                self.fragments[other_frag] = None

        gene_index = 0
        for frag in self.fragments:
            if frag is not None:
                start_index = frag.end1
                end_index = frag.end2

                node_str = self.node_strings[frag.end1 - 1]
                if node_str.endswith("h"):
                    self.medians[gene_index] = "-" + node_str[0: len(node_str) - 1]
                else:
                    self.medians[gene_index] = node_str[0: len(node_str) - 1]

                start_index = get_gene_next_node(start_index)
                node_gene_index = find_gray_edge_node(start_index, self.gray_edge)

                while not (start_index == end_index or node_gene_index is None):
                    node_gene = self.node_strings[node_gene_index - 1]

                    if node_gene.startswith("h", len(node_gene) - 1):
                        self.medians[gene_index] = self.medians[gene_index] + "  -" + \
                                                  node_gene[0: len(node_gene) - 1]
                    else:
                        self.medians[gene_index] = self.medians[gene_index] + "  " + \
                                                  node_gene[0: len(node_gene) - 1]

                    start_index = get_gene_next_node(node_gene_index)
                    node_gene_index = find_gray_edge_node(start_index, self.gray_edge)

                gene_index += 1
