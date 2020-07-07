from __future__ import annotations

from copy import copy
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
        self.index_from: int = int()
        self.for_which_genome: int = int()
        self.priority: int = int()
        self.position: int = int()
        self.genome_1_path: Optional[PGMPath] = None
        self.genome_2_path: Optional[PGMPath] = None
        self.genome_3_path: Optional[PGMPath] = None
        self.gray_edge: Optional[PGMPath] = None

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
            self.genome_1_path = copy(path)

        if genome_here == self.genome_2_path.genome_head:
            self.genome_2_path = copy(path)

        if genome_here == self.genome_3_path.genome_head:
            self.genome_3_path = copy(path)
