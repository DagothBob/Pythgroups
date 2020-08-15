from numpy import ndarray

from Genome import Genome

"""
 GenomeInString class for use in various classes

 Based on GenomeInString.java from C.Zheng & D.Sankoff (2011)

 Author: Holger Jensen, Oskar Jensen
"""


class GenomeInString:
    """
    Attributes
    ----------
    genome : [str]
        List of chromosomes expressed as strings
    """

    def __init__(self, genome: Genome):
        """
        Constructor

        Parameters
        ----------
        genome
            Chromosomes to construct from
        """
        self.chromosomes: ndarray = genome.chromosomes
