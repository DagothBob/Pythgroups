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
    head : int
        Head gene
    tail : int
        Tail gene
    genome_head : int
        Genome for the head gene
    genome_tail : int
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
        self.head: int = head
        self.tail: int = tail
        self.genome_head: int = ghead
        self.genome_tail: int = gtail

    @staticmethod
    def connect(path1: PGMPath, path2: PGMPath, path_l: PGMPath, which_genome: int) -> Optional[PGMPath]:
        """
        Connect two PGMPaths

        Parameters
        ----------
        path1
            First path to connect
        path2
            Second path to connect
        path_l
            Path used for connecting two paths
        which_genome
            Which genome this path belongs to

        Returns
        -------
        Optional[PGMPath]
            New PGMPath from the given paths or None if they are unable to be connected
        """
        head: int
        tail: int
        genome_head: int
        genome_tail: int

        if path1.head == path_l.head and path1.genome_head == path_l.genome_head and \
                path2.head == path_l.tail and path2.genome_head == path_l.genome_tail:
            head = path1.tail
            genome_head = path1.genome_tail
            tail = path2.tail
            genome_tail = path2.genome_tail

            if head == 0 and genome_head == -1:
                head = path1.head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.head
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.head == path_l.tail and path1.genome_head == path_l.genome_tail and \
                path2.head == path_l.head and path2.genome_head == path_l.genome_head:
            head = path1.tail
            genome_head = path1.genome_tail
            tail = path2.tail
            genome_tail = path2.genome_tail

            if head == 0 and genome_head == -1:
                head = path1.head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.head
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.head == path_l.tail and path1.genome_head == path_l.genome_tail and \
                path2.tail == path_l.head and path2.genome_tail == path_l.genome_head:
            head = path1.tail
            genome_head = path1.genome_tail
            tail = path2.head
            genome_tail = path2.genome_head

            if head == 0 and genome_head == -1:
                head = path1.head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.tail
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.head == path_l.head and path1.genome_head == path_l.genome_head and \
                path2.tail == path_l.tail and path2.genome_tail == path_l.genome_tail:
            head = path1.tail
            genome_head = path1.genome_tail
            tail = path2.head
            genome_tail = path2.genome_head

            if head == 0 and genome_head == -1:
                head = path1.head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.tail
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.tail == path_l.head and path1.genome_tail == path_l.genome_head and \
                path2.head == path_l.tail and path2.genome_head == path_l.genome_tail:
            head = path1.head
            genome_head = path1.genome_head
            tail = path2.tail
            genome_tail = path2.genome_tail

            if head == 0 and genome_head == -1:
                head = path1.tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.head
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.tail == path_l.tail and path1.genome_tail == path_l.genome_tail and \
                path2.head == path_l.head and path2.genome_head == path_l.genome_head:
            head = path1.head
            genome_head = path1.genome_head
            tail = path2.tail
            genome_tail = path2.genome_tail

            if head == 0 and genome_head == -1:
                head = path1.tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.head
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.tail == path_l.head and path1.genome_tail == path_l.genome_head and \
                path2.tail == path_l.tail and path2.genome_tail == path_l.genome_tail:
            head = path1.head
            genome_head = path1.genome_head
            tail = path2.head
            genome_tail = path2.genome_head

            if head == 0 and genome_head == -1:
                head = path1.tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.tail
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        if path1.tail == path_l.tail and path1.genome_tail == path_l.genome_tail and \
                path2.tail == path_l.head and path2.genome_tail == path_l.genome_head:
            head = path1.head
            genome_head = path1.genome_head
            tail = path2.head
            genome_tail = path2.genome_head

            if head == 0 and genome_head == -1:
                head = path1.tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2.tail
                genome_tail = which_genome

            return PGMPath(head, tail, genome_head, genome_tail)

        return None
