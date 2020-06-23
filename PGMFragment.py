from __future__ import annotations

from typing import Optional

import PGMPath

"""                                
 Fragment for use in the PathGroups algorithm (Median problem) 
                                                               
 Based on PGMFragment.java from C. Zheng and D. Sankoff (2011) 
                                                               
 Author: Oskar Jensen                                          
"""


class PGMFragment:
    """
    Attributes
    ----------
    end1 : int
        One end of the fragment
    end2 : int
        Other end of the fragment
    """

    def __init__(self, end1: int, end2: int):
        """
        Constructor

        Parameters
        ----------
        end1
            One end of the fragment
        end2
            Other end of the fragment
        """
        self.end1: int = end1
        self.end2: int = end2

    @classmethod
    def from_fragment(cls, fragment: PGMFragment) -> PGMFragment:
        """
        Construct new PGMFragment from existing one

        Parameters
        ----------
        fragment
            To construct from
        """
        end1 = fragment.end1
        end2 = fragment.end2

        return cls(end1, end2)

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        str
            String representation of the object
        """
        return "end1: " + str(self.end1) + " | end2: " + str(self.end2)


def combine(f_path: PGMPath, f1: PGMFragment, f2: PGMFragment) -> Optional[PGMFragment]:
    """
    Combines two fragments using the path between them to work out which ends to use

    Parameters
    ----------
    f_path
        Path fragment belongs to
    f1
        First fragment to combine
    f2
        Second fragment to combine

    Returns
    -------
    PGMFragment
        New PGMFragment from the given fragments or None if they cannot be combined
    """
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
