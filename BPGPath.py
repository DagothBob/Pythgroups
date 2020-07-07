from __future__ import annotations


"""                            
 Path for use in a BreakPoint Graph for BPGDistance.   
                                                       
 Based on BPGPath.java from C.Zheng & D.Sankoff (2011) 
                                                       
 Author: Holger Jensen                                 
"""


class BPGPath:
    """
    Attributes
    ----------
    head : int
        Head gene of the path
    tail : int
        Tail gene of the path
    genome_head : int
        Genome identifier for the head gene
    genome_tail : int
        Genome identifier for the tail gene
    """

    def __init__(self, head: int, tail: int, ancestor_head: int, ancestor_tail: int):
        """
        Constructor

        Parameters
        ----------
        head
            Head gene
        tail
            Tail gene
        ancestor_head
            Genome which gene head belongs to
        ancestor_tail
            Genome which gene tail belongs to
        """
        self.head: int = head
        self.tail: int = tail
        self.genome_head: int = ancestor_head
        self.genome_tail: int = ancestor_tail

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        str
            String representation of the BPGPath
        """
        return "Head node is: " + str(self.head) + \
               "(" + str(self.genome_head) + ") | Tail node is: " + \
               str(self.tail) + "(" + str(self.genome_tail) + ")"
