import PGMPath


# # # # # # # # # # # # # # # # # # # # # # # # # #
# Used in TreeStructure, MedianIteration, PGMPath #
#                                                 #
# Based on PGMPathForAGenome.java                 #
# from C.Zheng & D.Sankoff (2011)                 #
#                                                 #
# Author: Holger Jensen                           #
# # # # # # # # # # # # # # # # # # # # # # # # # #
class PGMPathForAGenome:
    def __init__(self, paths: [PGMPath]):
        self.paths = paths
