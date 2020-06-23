from __future__ import annotations

import math
from typing import List
from typing import Optional

from ChoiceStructure import ChoiceStructure
from PGMFragment import PGMFragment
from PGMPath import PGMPath

"""                            
 Represents data related to medians (ancestors) in the Pathgroups algorithm, 
 with related static methods.

 Based on MedianData.java from C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""


def get_gene_next_node(end1: int) -> int:
    """
    If end1 is even then it's greater than end2, else it's less than end2 (I think).
    Recall the pattern described in __init__ for self.fragments[]: [(1, 2) (2, 1), (3, 4), (4, 3), ... ]::

        ...
        Some possible hints on why it returns what it returns:
        (2, 1) is just (1, 2) flipped, so end1 - 1 from (2, 1) == end1 from (1, 2)
        In the case of (1, 2), end1 + 1 from (1, 2) == end2 from (1, 2).

    Parameters
    ----------
    end1
        end1 from PGMFragment

    Returns
    -------
    int
        end1 - 1 if even, else end + 1
    """

    if end1 % 2 == 0:
        return end1 - 1
    return end1 + 1


def find_gray_edge_node(end1: int, gray_edge: List[PGMPath]) -> Optional[int]:
    """
    Find the gray edge end that matches the given fragment end

    Parameters
    ----------
    end1
        end1 from PGMFragment
    gray_edge
        The gray edge to search for

    Returns
    -------
    Optional[int]
        If end1 found in gray_edge, return the opposite of what it matched with (head/tail), else return None
    """

    for path in gray_edge:
        if path is not None:
            if end1 == path.head:
                return path.tail
            if end1 == path.tail:
                return path.head

    return None


def get_next_path(start_pos: int, paths: List[PGMPath]) -> Optional[PGMPath]:
    """
    Returns the first path it finds in paths. Renamed from getAPath() for accuracy

    Parameters
    ----------
    start_pos : int
        Where to start in the paths list
    paths : [PGMPath]
        List of paths to search

    Returns
    -------
        Optional[PGMPath]
            If a path is found, return that, else return None
    """
    for i in range(start_pos, len(paths)):
        if paths[i] is not None:
            return paths[i]
    return None


def good_cycle(end1: int, end2: int) -> bool:
    """
    A good cycle has either of its ends be -1,
    and if either are >= 0 then it's an incomplete path

    Parameters
    ----------
    end1 : int
        end1 of the path
    end2 : int
        end2 of the path

    Returns
    -------
    bool
        Good cycle or not
    """
    if end1 >= 0 or end2 >= 0:
        return False
    return end1 == -1 or end2 == -1


class MedianData:
    """ Represents data related to medians (ancestors) in the Pathgroups algorithm

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
    self.choice_structures : [Optional[ChoiceStructure]]
        List of each choice structure
    self.medians : [str]
        List of each median
    """

    def __init__(self, paths1: List[PGMPath],
                 paths2: List[PGMPath],
                 paths3: List[PGMPath],
                 gene_num: int,
                 which_genome: int,
                 genome1: int,
                 genome2: int,
                 genome3: int,
                 node_strings: List[str]):
        """
        Constructor

        Parameters
        ----------
        paths1 : [[PGMPath]]
            All paths for genome 1 in the tree structure
        paths2 : [[PGMPath]]
            All paths for genome 2 in the tree structure
        paths3 : [[PGMPath]]
            All paths for genome 3 in the tree structure
        genome1 : [int]
            The identifiers for leaf 1 (genome) in the tree structure
        genome2 : [int]
            The identifiers for leaf 2 (genome) in the tree structure
        genome3 : [int]
            The identifiers for leaf 3 (genome) in the tree structure
        gene_num : int
            Total number of genes in each genome
        which_genome : int
            Which genome in the tree structure these medians exist in
        node_strings : [str]
            String representation of each node
        """
        self.three_genome_paths: List[List[PGMPath]] = [paths1, paths2, paths3]
        self.three_genomes: List[int] = [genome1, genome2, genome3]
        self.gene_num: int = gene_num
        self.which_genome: int = which_genome
        self.node_strings: List[str] = node_strings

        self.gray_edge: List[PGMPath] = list()
        self.gray_edge_index: int = int()
        self.fragments: List[Optional[PGMFragment]] = list()
        self.choice_structures: List[Optional[ChoiceStructure]] = list()
        self.medians: List[str] = list()

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

            self.choice_structures[cs_index] = cs
            cs_index += 1

    def gray_edge_total_distance(self, median: MedianData) -> int:
        """
        Calculate the total distances between the three genome paths from the given median
        and this instance's gray edge. Renamed from countTotalDis

        Parameters
        ----------
        median : MedianData
            Median to calculate distances from

        Returns
        -------
        int
            Total distance between three genome paths from the given median and this instance's gray edge
        """
        d1 = self.get_distance(median.three_genome_paths[0], self.gray_edge)
        d2 = self.get_distance(median.three_genome_paths[1], self.gray_edge)
        d3 = self.get_distance(median.three_genome_paths[2], self.gray_edge)
        return d1 + d2 + d3

    def get_distance(self, p1: List[PGMPath], p2: List[PGMPath]) -> int:
        """
        Gets the distance between two sets of paths

        Parameters
        ----------
        p1 : [PGMPath]
            Path 1
        p2 : [PGMPath]
            Path 2

        Returns
        -------
        int
            The distance between two path sets, calculated as:
                (total number of genes) + (number of chromosomes) - (number of cycles) - (number of "good" paths)
        """

        cycle_count, good_path_count, chr_count = 0, 0, 0
        paths1 = list()

        # Populates paths1 based on the heads and tails of each path in p1
        # where negative values are set to -1
        for path in p1:
            if path is not None:
                end1 = path.head
                end2 = path.tail
                if end1 < 0:
                    end1 = -1
                if end2 < 0:
                    end2 = -1
                if end1 > 0:
                    paths1[end1] = PGMPath(end1, end2, 1, 1)
                if end2 > 0:
                    paths1[end2] = PGMPath(end2, end1, 1, 1)

        # Counts the number of chromosomes in paths1
        # Fills in None values with paths where head is i and tail is -1
        # Count gets divided by 2 in the end I think because there are two paths for each gene?
        for i in range(0, len(paths1)):
            if paths1[i] is not None and paths1[i].tail == -1:
                chr_count += 1
            elif paths1[i] is None and i != 0:
                paths1[i] = PGMPath(i, -1, 1, 1)
                chr_count += 1
        chr_count = math.trunc(chr_count / 2)

        # Same as paths1, except without chr_count
        paths2 = List[PGMPath]
        for path in p2:
            if path is not None:
                end1 = path.head
                end2 = path.tail
                if end1 < 0:
                    end1 = -1
                if end2 < 0:
                    end2 = -1
                if end1 > 0:
                    paths2[end1] = PGMPath(end1, end2, 2, 2)
                if end2 > 0:
                    paths2[end2] = PGMPath(end2, end1, 2, 2)

        for i in range(0, len(paths1)):
            if paths2[i] is None and i != 0:
                paths2[i] = PGMPath(i, -2, 2, 2)

        # Iterates through paths1 until all paths have been checked
        next_path = get_next_path(0, paths1)
        while next_path is not None:
            n1 = next_path.head
            n2 = next_path.tail
            start_point = n1

            if n1 > 0:
                paths1[n1] = None
            if n2 > 0:
                paths1[n2] = None

            if n1 > n2:
                n_big = n1
                n_small = n2
            else:
                n_big = n2
                n_small = n1

            more = True
            while more:
                if next_path.head < 0 and next_path.tail < 0:
                    if good_cycle(next_path.head, next_path.tail):
                        good_path_count += 1
                    break

                path_l = paths2[n_big]
                m1 = path_l.head
                m2 = path_l.tail
                if m2 > 0:
                    if m2 == n_small:
                        cycle_count += 1
                        more = False
                        paths2[m1] = None
                        paths2[m2] = None
                    else:
                        m2_path = paths1[m2]
                        next_path = PGMPath(m2_path.tail, n_small, 1, 1)
                        n1 = next_path.head
                        n2 = next_path.tail
                        if n1 > n2:
                            n_big = n1
                            n_small = n2
                        else:
                            n_big = n2
                            n_small = n1

                        other_path = paths1[m2].tail
                        if other_path > 0:
                            paths1[other_path] = None
                        paths2[m1] = None
                        paths2[m2] = None
                        paths1[m2] = None
                elif m2 < 0:
                    if n_small < 0:
                        if good_cycle(m2, n_small):
                            good_path_count += 1
                        more = False
                        paths2[m1] = None
                    elif n_small > 0:
                        next_path = PGMPath(n_small, m2, 1, 2)
                        n1 = next_path.head
                        n2 = next_path.tail
                        if n1 > n2:
                            n_big = n1
                            n_small = n2
                        else:
                            n_big = n2
                            n_small = n1
                        paths2[m1] = None

            next_path = get_next_path(start_point, paths1)

        return self.gene_num + chr_count - cycle_count - good_path_count

    def get_ancestors(self):
        """
        Fills medians[] with all the genome's medians (ancestors),
        which are formatted for the purpose of printing to screen.
        """

        for i in range(0, len(self.fragments)):
            if self.fragments[i] is not None:
                other_frag: int = self.fragments[i].end2
                self.fragments[other_frag] = None

        gene_index: int = 0
        for frag in self.fragments:
            if frag is not None:
                start_index: int = frag.end1
                end_index: int = frag.end2

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
