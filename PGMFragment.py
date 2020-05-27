from __future__ import annotations
import PGMPath


class PGMFragment:
    def __init__(self, end1: int, end2: int):
        self.end1 = end1
        self.end2 = end2

    @classmethod
    def from_fragment(cls, fragment: PGMFragment):
        end1 = fragment.end1
        end2 = fragment.end2
        return cls(end1, end2)

    def str(self):
        return str(self.end1) + " " + str(self.end2)

    # TODO: once PGMPath is complete, work on this
    def combine(self, a1: PGMPath, f1: PGMFragment, f2: PGMFragment):
        return None
