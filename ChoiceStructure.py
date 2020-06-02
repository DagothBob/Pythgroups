from __future__ import annotations

import PGMPath


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Used in MedianData for solving the rearrangement median problem.#
#                                                                 #
# Based on ChoiceStructure.java from C.Zheng & D.Sankoff (2011)   #
#                                                                 #
# Author: Holger Jensen, Oskar Jensen                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class ChoiceStructure:
    def __init__(self):
        self.index_from = 0        # Object instance identifier
        self.for_which_genome = 0  # For which genome this structure refers to
        self.priority = 0          # Priority level for makeCycles algorithm (section 3.1.3)
        self.position = 0          # Position of the choice structure in Priority's cs_index[]
        self.genome_1_path = None  # Path for genome 1
        self.genome_2_path = None  # Path for genome 2
        self.genome_3_path = None  # Path for genome 3
        self.gray_edge = None      # Gets used when updating paths & priorities in the makeCycles algorithm

    @classmethod
    def from_cs(cls, cs: ChoiceStructure):
        cls.index_from = cs.index_from
        cls.for_which_genome = cs.for_which_genome
        cls.priority = cs.priority
        cls.position = cs.position
        cls.genome_1_path = cs.genome_1_path
        cls.genome_2_path = cs.genome_2_path
        cls.genome_3_path = cs.genome_3_path
        cls.gray_edge = cs.gray_edge

    def __setattr__(self, key, value):
        self.key = value

    def set_new_path(self, path: PGMPath):
        genome_here = path.genome_head

        if genome_here == self.genome_1_path.genome_head:
            self.genome_1_path = path

        if genome_here == self.genome_2_path.genome_head:
            self.genome_2_path = path

        if genome_here == self.genome_3_path.genome_head:
            self.genome_3_path = path
