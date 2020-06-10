"""
 GenomeInString class for use in various classes

 Based on GenomeInString.java from C.Zheng & D.Sankoff (2011)

 Author: Holger Jensen
"""


class GenomeInString:
    """
    Attributes:
        chromosomes: List of chromosomes expressed as strings
    """
    def __init__(self, chromosomes: [str]):
        """
        Constructor

        :param chromosomes: Chromosomes to construct from
        """
        self.chromosomes: [str] = chromosomes
