from __future__ import annotations


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Path for use in the PathGroups algorithm (Median problem) #
#                                                           #
# Based on PGMPath.java from C.Zheng & D.Sankoff (2011)     #
#                                                           #
# Author: Holger Jensen                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class PGMPath:
    def __init__(self, head: int, tail: int, ghead: int, gtail: int):
        self.head = head          # Head gene
        self.tail = tail          # Tail gene
        self.genome_head = ghead  # Genome for head
        self.genome_tail = gtail  # Genome for tail

    # Connect two PGMPaths
    def connect(self, path1: PGMPath, path2: PGMPath, pathl: PGMPath, which_genome: int) -> PGMPath or None:
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
