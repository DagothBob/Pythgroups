from typing import List


"""
 GenomeInString class for use in various classes

 Based on GenomeInString.java from C.Zheng & D.Sankoff (2011)

 Author: Holger Jensen, Oskar Jensen
"""


class GenomeInString:
    """
    Attributes
    ----------
    chromosomes : [str]
        List of chromosomes expressed as strings
    """

    def __init__(self, chromosomes: List[str]):
        """
        Constructor

        Parameters
        ----------
        chromosomes
            Chromosomes to construct from
        """
        self.chromosomes: List[str] = chromosomes
