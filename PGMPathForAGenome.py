from PGMPath import PGMPath


"""                         
 Used in TreeStructure, MedianIteration, PGMPath 
                                                 
 Based on PGMPathForAGenome.java                 
 from C.Zheng & D.Sankoff (2011)                 
                                                 
 Author: Holger Jensen                           
"""


class PGMPathForAGenome:
    """
    Attributes
    ----------
    paths
        PGMPaths for the genome
    """
    def __init__(self, paths: [PGMPath]):
        """
        Constructor

        Parameters
        ----------
        paths
            PGMPaths for the genome
        """
        self.paths: [PGMPath] = paths
