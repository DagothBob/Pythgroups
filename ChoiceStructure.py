from __future__ import annotations

import PGMPath

"""                                 
 Used in MedianData for solving the rearrangement median problem.
                                                                 
 Based on ChoiceStructure.java from C.Zheng & D.Sankoff (2011)   
                                                                 
 Author: Holger Jensen, Oskar Jensen                             
"""


class ChoiceStructure:
    """
    Attributes:
        index_from: Object instance identifier
        for_which_genome: For which genome this structure refers to
        priority: Priority level for makeCycles algorithm
        position: Position of the choice structure in Priority's cs_index
        genome_1_path: Path for genome 1
        genome_2_path: Path for genome 2
        genome_3_path: Path for genome 3
        gray_edge: Used when updating paths and priorities in the makeCycles algorithm
    """
    def __init__(self):
        """
        Constructor
        """
        self.index_from = 0
        self.for_which_genome = 0
        self.priority = 0
        self.position = 0
        self.genome_1_path = None
        self.genome_2_path = None
        self.genome_3_path = None
        self.gray_edge = None

    @classmethod
    def from_cs(cls, cs: ChoiceStructure):
        """
        Constructs a new ChoiceStructure from an existing given one

        :param cs: ChoiceStructure to copy from
        :return: New ChoiceStructure
        """
        cls.index_from = cs.index_from
        cls.for_which_genome = cs.for_which_genome
        cls.priority = cs.priority
        cls.position = cs.position
        cls.genome_1_path = cs.genome_1_path
        cls.genome_2_path = cs.genome_2_path
        cls.genome_3_path = cs.genome_3_path
        cls.gray_edge = cs.gray_edge

    def set_new_path(self, path: PGMPath):
        """
        Sets instance paths to the given PGMPath if its genome matches

        :param path: PGMPath to copy from
        """
        genome_here = path.genome_head

        if genome_here == self.genome_1_path.genome_head:
            self.genome_1_path = path

        if genome_here == self.genome_2_path.genome_head:
            self.genome_2_path = path

        if genome_here == self.genome_3_path.genome_head:
            self.genome_3_path = path
