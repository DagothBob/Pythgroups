from random import random
from typing import List

from MedianData import MedianData
from PGMPath import PGMPath
from PGMPathForAGenome import PGMPathForAGenome
from SmallPhylogeny import SmallPhylogeny
from TreeStructure import TreeStructure

"""
 Used in result optimization in the Pathgroups algorithm

 Based on MedianIteration.java from C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""


class MedianIteration:
    """ Used in result optimization in the Pathgroups algorithm

    Attributes
    ----------
    leaf_num : int
        Number of leaves in the tree
    changed : List[int]
        List of changed ancestors. Unsure of what its purpose is, as it's not used anywhere meaningfully.
    gene_num : int
        Total number of genes in the tree
    leaves : List[List[int]]
        List of each leaf for each genome in the tree (in relation to the median)
    all_paths : List[PGMPathForAGenome]
        List of all paths for each genome in the tree
    medians : List[MedianData]
        List of all medians in the tree
    nodes_int : List[int]
        List of the integer representations of each node
    nodes_str : List[str]
        List of the string representations of each node
    """

    def __init__(self, leaf_num: int,
                 ancestor_num: int,
                 gene_num: int,
                 leaves: List[List[int]],
                 all_paths: List[PGMPathForAGenome],
                 medians: List[MedianData],
                 nodes_int: List[int],
                 nodes_str: List[str]):
        """
        Parameters
        ----------
        leaf_num : int
            Number of leaves in the tree
        ancestor_num : int
            Number of ancestors (medians) in the tree
        gene_num : int
            Total number of genes in the tree
        leaves : List[List[int]]
            List of each leaf for each genome in the tree (in relation to the median)
        all_paths : List[PGMPathForAGenome]
            List of all paths for each genome in the tree
        medians : List[MedianData]
            List of all medians in the tree
        nodes_int : List[int]
            List of the integer representations of each node
        nodes_str : List[str]
            List of the string representations of each node
        """
        self.leaf_num: int = leaf_num
        self.changed: List[int] = [1 for _ in range(ancestor_num)]
        self.gene_num: int = gene_num
        self.leaves: List[List[int]] = [row[:] for row in leaves]  # 2D array copy
        self.all_paths: List[PGMPathForAGenome] = all_paths.copy()
        self.medians: List[MedianData] = medians.copy()
        self.nodes_int: List[int] = nodes_int.copy()
        self.nodes_str: List[str] = nodes_str.copy()

    def median_total_distance(self, median_paths: List[PGMPath], paths1: List[PGMPath],
                              paths2: List[PGMPath], paths3: List[PGMPath]) -> int:
        """
        Returns the sum of all distances between the paths of the given median and three leaves
        Renamed from countTotalDistance

        Parameters
        ----------
        median_paths : [PGMPath]
            List of paths for the given median
        paths1 : [PGMPath]
            List of paths for the first leaf
        paths2 : [PGMPath]
            List of paths for the second leaf
        paths3 : [PGMPath]
            List of paths for the third leaf

        Returns
        -------
        int
            Total distance between median and 3 leaves
        """
        d1 = self.medians[0].get_distance(median_paths, paths1)
        d2 = self.medians[0].get_distance(median_paths, paths2)
        d3 = self.medians[0].get_distance(median_paths, paths3)
        return d1 + d2 + d3

    def optimize_result(self, prob_threshold_top: float, runs: int):
        """
        Optimizes the results of the ancestor reconstruction

        Renamed from optimizeResultTwoStepLA

        Parameters
        ----------
        prob_threshold_top : float
            Probability threshold for when total_dis == dis_before
        runs : int
            Total number of optimization runs
        """
        p1 = prob_threshold_top  # Probability threshold for when total_dis == dis_before
        p2 = 0  # Probability threshold for when total_dis > dis_before

        # Each optimization run
        for turn in range(0, runs):
            # Iterates through each changed ancestor
            for i in range(0, len(self.changed)):
                leaf1 = self.leaves[i + self.leaf_num][0]
                leaf2 = self.leaves[i + self.leaf_num][1]
                leaf3 = self.leaves[i + self.leaf_num][2]

                dis_before = self.median_total_distance(self.all_paths[i + self.leaf_num].paths,
                                                        self.all_paths[leaf1].paths,
                                                        self.all_paths[leaf2].paths,
                                                        self.all_paths[leaf3].paths)
                aps = List[PGMPathForAGenome]  # Local copy of all_paths containing the given leaves
                aps[0] = self.all_paths[leaf1]
                aps[1] = self.all_paths[leaf2]
                aps[2] = self.all_paths[leaf3]
                # Filters out all None values for each list of paths in each leaf
                for t in range(0, len(aps[0].paths)):
                    if aps[0].paths[t] is not None:
                        aps[0].paths[t] = PGMPath(aps[0].paths[t].head, aps[0].paths[t].tail, 0, 0)
                    if aps[1].paths[t] is not None:
                        aps[1].paths[t] = PGMPath(aps[1].paths[t].head, aps[1].paths[t].tail, 1, 1)
                    if aps[2].paths[t] is not None:
                        aps[2].paths[t] = PGMPath(aps[2].paths[t].head, aps[2].paths[t].tail, 2, 2)

                ts = TreeStructure(1, 3, self.gene_num, aps, self.nodes_str, self.nodes_int, None)
                sp = SmallPhylogeny(ts)
                sp.get_result()

                # Specifies thresholds to meet to perform optimizations
                total_dis = ts.medians[0].gray_edge_total_distance(ts.medians[0])
                rand1 = random()
                rand2 = random()
                total_less = total_dis < dis_before
                total_equal = total_dis == dis_before and rand1 < p1
                total_greater = total_dis > dis_before and rand2 < p2
                if total_less or total_equal or total_greater:
                    # Not sure what the point of these statements are, self.changed isn't used anywhere really
                    if leaf1 >= self.leaf_num:
                        self.changed[leaf1 - self.leaf_num] = 1
                    if leaf2 >= self.leaf_num:
                        self.changed[leaf2 - self.leaf_num] = 1
                    if leaf3 >= self.leaf_num:
                        self.changed[leaf3 - self.leaf_num] = 1

                    self.medians[i] = ts.medians[0]
                    self.all_paths[i + self.leaf_num].paths = List[PGMPath]
                    # Performs optimizations on each path in gray_edge
                    for g in range(0, len(self.medians[i].gray_edge)):
                        if self.medians[i].gray_edge[g] is not None:
                            end1 = self.medians[i].gray_edge[g].head
                            end2 = self.medians[i].gray_edge[g].tail
                            if end1 < 0:
                                end1 = -1
                            if end2 < 0:
                                end2 = -1
                            if end1 > 0:
                                self.all_paths[i + self.leaf_num].paths[end1] = PGMPath(end1, end2, 1, 1)
                            if end2 > 0:
                                self.all_paths[i + self.leaf_num].paths[end2] = PGMPath(end2, end1, 1, 1)
                    for g in range(0, len(self.all_paths[i + self.leaf_num].paths)):
                        if self.all_paths[i + self.leaf_num].paths[g] is None and g != 0:
                            self.all_paths[i + self.leaf_num].paths[g] = PGMPath(g, -1, 1, 1)

        # Get the ancestors and distances
        for median in self.medians:
            median.get_ancestors()
