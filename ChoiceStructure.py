from __future__ import annotations

from typing import Optional

from PGMPath import PGMPath

"""                                 
 Used in MedianData for solving the rearrangement median problem.
                                                                 
 Based on ChoiceStructure.java from C.Zheng & D.Sankoff (2011)   
                                                                 
 Author: Holger Jensen, Oskar Jensen                             
"""


class ChoiceStructure:
    """
    Attributes
    ----------
    index_from : int
        Object instance identifier
    for_which_genome : int
        For which genome this structure refers to
    priority : int
        Priority level for makeCycles algorithm
    position : int
        Position of the choice structure in Priority's cs_index
    genome_1_path : Optional[PGMPath]
        Path for genome 1
    genome_2_path : Optional[PGMPath]
        Path for genome 2
    genome_3_path : Optional[PGMPath]
        Path for genome 3
    gray_edge : Optional[PGMPath]
        Used when updating paths and priorities in the makeCycles algorithm
    """

    def __init__(self):
        """
        Constructor
        """
        self.index_from: int = 0
        self.for_which_genome: int = 0
        self.priority: int = 0
        self.position: int = 0
        self.genome_1_path: Optional[PGMPath] = None
        self.genome_2_path: Optional[PGMPath] = None
        self.genome_3_path: Optional[PGMPath] = None
        self.gray_edge: Optional[PGMPath] = None

    @classmethod
    def from_cs(cls, cs: ChoiceStructure):
        """
        Constructs a new ChoiceStructure from an existing given one

        Parameters
        ----------
        cs
            ChoiceStructure to copy from
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

        Parameters
        ----------
        path
            PGMPath to copy from
        """
        genome_here: int = path.genome_head

        if genome_here == self.genome_1_path.genome_head:
            self.genome_1_path = path

        if genome_here == self.genome_2_path.genome_head:
            self.genome_2_path = path

        if genome_here == self.genome_3_path.genome_head:
            self.genome_3_path = path
