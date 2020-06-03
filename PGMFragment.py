from __future__ import annotations
import PGMPath


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fragment for use in the PathGroups algorithm (Median problem) #
#                                                               #
# Based on PGMFragment.java from C. Zheng and D. Sankoff (2011) #
#                                                               #
# Author: Oskar Jensen                                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class PGMFragment:
    def __init__(self, end1: int, end2: int):
        self.end1 = end1
        self.end2 = end2

    @classmethod
    def from_fragment(cls, fragment: PGMFragment):
        end1 = fragment.end1
        end2 = fragment.end2
        return cls(end1, end2)

    def __str__(self):
        return "end1: " + str(self.end1) + " | end2: " + str(self.end2)


# Combines two fragments using the path between them to work out which ends to use
def combine(f_path: PGMPath, f1: PGMFragment, f2: PGMFragment) -> PGMFragment or None:
    h = f_path.head
    t = f_path.tail
    if (h == f1.end1 and t == f2.end1) or (t == f1.end1 and h == f2.end1):
        return PGMFragment(f1.end2, f2.end2)

    elif (h == f1.end1 and t == f2.end2) or (t == f1.end1 and h == f2.end2):
        return PGMFragment(f1.end2, f2.end1)

    elif (h == f1.end2 and t == f2.end1) or (t == f1.end2 and h == f2.end1):
        return PGMFragment(f1.end1, f2.end2)

    elif (h == f1.end2 and t == f2.end2) or (t == f1.end2 and h == f2.end2):
        return PGMFragment(f1.end1, f2.end1)

    return None
