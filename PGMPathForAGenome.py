from typing import List

from PGMPath import PGMPath

"""                         
 Used in TreeStructure, MedianIteration, PGMPath 
                                                 
 Based on PGMPathForAGenome.java from C.Zheng & D.Sankoff (2011)                 
                                                 
 Author: Holger Jensen, Oskar Jensen                        
"""


class PGMPathForAGenome:
    """
    Attributes
    ----------
    paths : List[PGMPath]
        PGMPaths for the genome
    """

    def __init__(self, paths: List[PGMPath]):
        """
        Constructor

        Parameters
        ----------
        paths
            PGMPaths for the genome
        """
        self.paths: List[PGMPath] = paths
