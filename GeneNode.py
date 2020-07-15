"""
For use in DCJRearrangement representing a Gene Node

Author: Holger Jensen
"""


class GeneNode:
    """
    Attributes
    ----------
    adjacency: int
        GeneNode it is adjacent to
    chromosome_index: int
        Its chromosomes index in the genome
    chromosome_position: int
        Its location in the chromosome
    linked_adjacency: int
        Unused
    """
    def __init__(self):
        """
        Constructor
        """
        self.adjacency: int = 0
        self.chromosome_index: int = 0
        self.chromosome_position: int = 0
        self.linked_adjacency: int = 0

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        String representation of the object
        """
        return "Adjacency: " + str(self.adjacency) + "\n" + \
            "Chromosome index: " + str(self.chromosome_index) + "\n" + \
            "Chromosome position: " + str(self.chromosome_position) + "\n" + \
            "Linked adjacency: " + str(self.linked_adjacency) + "\n"
