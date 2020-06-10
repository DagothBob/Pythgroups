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
    head
        Head gene of the path
    tail
        Tail gene of the path
    genome_head
        Genome identifier for the head gene
    genome_tail
        Genome identifier for the tail gene
    """
    def __init__(self, head: int or None, tail: int or None, ghead: int or None, gtail: int or None):
        """
        Constructor

        Parameters
        ----------
        head
            Head gene
        tail
            Tail gene
        ghead
            Genome which gene head belongs to
        gtail
            Genome which gene tail belongs to
        """
        if head is None or tail is None or ghead is None or gtail is None:
            return
        else:
            self.head = head
            self.tail = tail
            self.genome_head = ghead
            self.genome_tail = gtail

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
