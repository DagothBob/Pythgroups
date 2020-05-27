from __future__ import annotations


# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Path for use in a BreakPoint Graph for BPGDistance.   #
#                                                       #
# Based on BPGPath.java from C.Zheng & D.Sankoff (2011) #
#                                                       #
# Author: Holger Jensen                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class BPGPath:
    def __init__(self, head: int, tail: int, ghead: int, gtail: int):
        self.head = head          # Head gene
        self.tail = tail          # Tail gene
        self.genome_head = ghead  # Genome for head
        self.genome_tail = gtail  # Genome for tail

    def __str__(self):
        return "Head node is: " + str(self.head) + \
               "(" + str(self.genome_head) + ") | Tail node is: " + \
               str(self.tail) + "(" + str(self.genome_tail) + ")"
