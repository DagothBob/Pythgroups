from random import random
from typing import List, Dict

import PGMPath
from MedianData import MedianData
from PGMPathForAGenome import PGMPathForAGenome
from SmallPhylogeny import SmallPhylogeny
from TreeStructure import TreeStructure

"""
 Used in result optimization in the Pathgroups algorithm

 Based on MedianIteration.java from C.Zheng & D.Sankoff (2011)

 Author: Holger Jensen, Oskar Jensen
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
        self.leaves: List[List[int]] = leaves
        self.all_paths: List[PGMPathForAGenome] = all_paths
        self.medians: List[MedianData] = medians
        self.nodes_int: List[int] = nodes_int
        self.nodes_str: List[str] = nodes_str

    def median_total_distance(self,
                              median_paths: List[Dict[str, int]],
                              paths1: List[Dict[str, int]],
                              paths2: List[Dict[str, int]],
                              paths3: List[Dict[str, int]]) -> int:
        """
        Returns the sum of all distances between the paths of the given median and three leaves
        Renamed from countTotalDistance

        Parameters
        ----------
        median_paths : List[Dict[str, int]]
            List of paths for the given median
        paths1 : List[Dict[str, int]]
            List of paths for the first leaf
        paths2 : List[Dict[str, int]]
            List of paths for the second leaf
        paths3 : List[Dict[str, int]]
            List of paths for the third leaf

        Returns
        -------
        int
            Total distance between median and 3 leaves
        """
        return self.medians[0].get_distance(median_paths, paths1) + \
            self.medians[0].get_distance(median_paths, paths2) + \
            self.medians[0].get_distance(median_paths, paths3)

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
        p1: float = prob_threshold_top  # Probability threshold for when total_dis == dis_before
        p2: float = 0  # Probability threshold for when total_dis > dis_before

        # Each optimization run
        for turn in range(runs):
            # Iterates through each changed ancestor
            for i in range(len(self.changed)):
                print("\rOptimizating result: run {}/{}".format(turn + 1, runs), end="")

                leaf1: int = self.leaves[i + self.leaf_num][0]
                leaf2: int = self.leaves[i + self.leaf_num][1]
                leaf3: int = self.leaves[i + self.leaf_num][2]

                dis_before: int = self.median_total_distance(self.all_paths[i + self.leaf_num].paths,
                                                             self.all_paths[leaf1].paths,
                                                             self.all_paths[leaf2].paths,
                                                             self.all_paths[leaf3].paths)

                # Local copy of all_paths containing the given leaves
                aps: List[PGMPathForAGenome] = [self.all_paths[leaf1],
                                                self.all_paths[leaf2],
                                                self.all_paths[leaf3]]

                # Filters out all None values for each list of paths in each leaf
                for t in range(len(aps[0].paths)):
                    if aps[0].paths[t] is not None:
                        aps[0].paths[t] = PGMPath.create_pgm_path(aps[0].paths[t]["head"],
                                                                  aps[0].paths[t]["tail"],
                                                                  0,
                                                                  0)
                    if aps[1].paths[t] is not None:
                        aps[1].paths[t] = PGMPath.create_pgm_path(aps[1].paths[t]["head"],
                                                                  aps[1].paths[t]["tail"],
                                                                  1,
                                                                  1)
                    if aps[2].paths[t] is not None:
                        aps[2].paths[t] = PGMPath.create_pgm_path(aps[2].paths[t]["head"],
                                                                  aps[2].paths[t]["tail"],
                                                                  2,
                                                                  2)

                ts: TreeStructure = TreeStructure(1, 3, self.gene_num, aps, self.nodes_str, self.nodes_int)
                sp: SmallPhylogeny = SmallPhylogeny(ts)
                sp.get_result()

                # Specifies thresholds to meet to perform optimizations
                total_dis: int = ts.medians[0].gray_edge_total_distance(ts.medians[0])
                rand1: float = random()
                rand2: float = random()
                total_less: bool = total_dis < dis_before
                total_equal: bool = total_dis == dis_before and rand1 < p1
                total_greater: bool = total_dis > dis_before and rand2 < p2

                if total_less or total_equal or total_greater:
                    # Not sure what the point of these statements are, self.changed isn't used anywhere really
                    if leaf1 >= self.leaf_num:
                        self.changed[leaf1 - self.leaf_num] = 1
                    if leaf2 >= self.leaf_num:
                        self.changed[leaf2 - self.leaf_num] = 1
                    if leaf3 >= self.leaf_num:
                        self.changed[leaf3 - self.leaf_num] = 1

                    self.medians[i] = ts.medians[0]
                    self.all_paths[i + self.leaf_num].paths = [None for _ in range(2 * self.gene_num + 1)]

                    # Performs optimizations on each path in gray_edge
                    for g in range(len(self.medians[i].gray_edge)):
                        if self.medians[i].gray_edge[g] is not None:
                            end1: int = self.medians[i].gray_edge[g]["head"]
                            end2: int = self.medians[i].gray_edge[g]["tail"]

                            if end1 < 0:
                                end1 = -1
                            if end2 < 0:
                                end2 = -1
                            if end1 > 0:
                                self.all_paths[i + self.leaf_num].paths[end1] = PGMPath.create_pgm_path(end1,
                                                                                                        end2,
                                                                                                        1,
                                                                                                        1)
                            if end2 > 0:
                                self.all_paths[i + self.leaf_num].paths[end2] = PGMPath.create_pgm_path(end2,
                                                                                                        end1,
                                                                                                        1,
                                                                                                        1)

                    for g in range(len(self.all_paths[i + self.leaf_num].paths)):
                        if self.all_paths[i + self.leaf_num].paths[g] is None and g != 0:
                            self.all_paths[i + self.leaf_num].paths[g] = PGMPath.create_pgm_path(g, -1, 1, 1)

        # Get the ancestors and distances
        for median in self.medians:
            median.get_ancestors()
