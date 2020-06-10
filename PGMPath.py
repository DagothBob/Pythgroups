from __future__ import annotations

from typing import Optional

"""                              
 Path for use in the PathGroups algorithm (Median problem) 
                                                           
 Based on PGMPath.java from C.Zheng & D.Sankoff (2011)     
                                                           
 Author: Holger Jensen & Oskar Jensen                                     
"""


class PGMPath:
    """
    Attributes
    ----------
    head
        Head gene
    tail
        Tail gene
    genome_head
        Genome for the head gene
    genome_tail
        Genome for the tail gene
    """
    def __init__(self, head: int, tail: int, ghead: int, gtail: int):
        """
        Constructor

        Parameters
        ----------
        head
            Head gene
        tail
            Tail gene
        ghead
            Genome for the head gene
        gtail
            Genome for the tail gene
        """
        self.head = head
        self.tail = tail
        self.genome_head = ghead
        self.genome_tail = gtail

    @staticmethod
    def connect(path1: PGMPath, path2: PGMPath, pathl: PGMPath, which_genome: int) -> Optional[PGMPath]:
        """
        Connect two PGMPaths

        Parameters
        ----------
        path1
            First path to connect
        path2
            Second path to connect
        pathl
            TODO: Figure out what this means
        which_genome
            Which genome this path belongs to

        Returns
        -------
        Optional[PGMPath]
            New PGMPath from the given paths or None if they are unable to be connected
        """
        h = 0
        t = 0
        gh = 0
        gt = 0

        # Step 1: Assigning the new head to path 1
        if path1.head == pathl.head or path1.head == pathl.tail:
            h = path1.tail
            gh = path1.genome_tail

            if h == 0 and gh == -1:
                h = path1.head
                gh = which_genome

        elif path1.tail == pathl.tail or path1.tail == pathl.head:
            h = path1.head
            gh = path1.genome_head

            if h == 0 and gh == -1:
                h = path1.tail
                gh = which_genome

        # Step 2: Assigning the new tail to path 2
        if path2.head == pathl.head or path2.head == pathl.tail:
            t = path2.tail
            gt = path2.genome_tail

            if t == 0 and gt == -1:
                t = path2.head
                gt = which_genome

        elif path2.tail == pathl.tail or path2.tail == pathl.head:
            t = path2.head
            gt = path2.genome_head

            if t == 0 and gt == -1:
                t = path2.tail
                gt = which_genome

        if h == 0 and t == 0 and gh == 0 and gt == 0:
            return None
        else:
            return PGMPath(h, t, gt, gh)
