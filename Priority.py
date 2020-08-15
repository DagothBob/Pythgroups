from typing import List, Optional

from numpy import full as npfull, empty as npempty, array as nparray, ndarray, arange, int32

"""
Priority for use in the PathGroups Algorithm

Based on Priority.java from C. Zheng and D. Sankoff (2011)

Author: Holger Jensen, Oskar Jensen
"""


class Priority:
    """
    Based on the priority system described in 'On the PATHGROUPS approach to rapid small phylogeny'

    Attributes
    ----------
    size : int
        Size of each array in this class
    cs_indexes : List[int]
        Position of choice structures
    median_indexes : List[int]
        Position of medians (ancestors)
    empty_start : int
        Starting empty position
    empty_end : int
        Ending empty position
    taken_start : int
        Starting taken position
    taken_end : int
        Ending taken position
    empty_next : List[int]
        Next empty values
    taken_previous : List[int]
        Previous taken values
    taken_next : List[int]
        Next taken values

    """

    def __init__(self, size: int):
        """
        Constructor

        Parameters
        ----------
        size
            Size of each array in this class
        """
        self.cs_indexes: ndarray = npfull(size, -1, int32)
        self.median_indexes: ndarray = npfull(size, -1, int32)
        self.taken_previous: ndarray = npempty(size, int32)
        self.taken_next: ndarray = npempty(size, int32)
        self.empty_next: ndarray = nparray([i + 1 for i in arange(size)], int32)

        self.empty_next[size - 1] = -1     # last entry set to -1

        self.empty_start: int = 0
        self.empty_end: int = size - 1
        self.taken_start: int = -1
        self.taken_end: int = -1

    def insert(self, cs_index: int, which_ancestor: Optional[int] = None) -> int:
        """
        Insert a choice structure position into the Priority object

        Parameters
        ----------
        cs_index : int
            Position of the choice structure in the median
        which_ancestor : int
            Position of the ancestor in the tree

        Returns
        -------
        int
            The position of empty_start, which shifts to the next empty position each time insert() is called
        """

        cur_pos: int = self.empty_start
        first_empty_pos: int = self.empty_next[cur_pos]

        # Update cs/median_indexes and empty_start/previous
        if which_ancestor is not None:  # Small phylogeny case
            self.median_indexes[self.empty_start] = which_ancestor

        self.cs_indexes[self.empty_start] = cs_index
        self.empty_start = first_empty_pos

        # Update various taken values based on whether or not
        # the taken_start/end values at the current empty start position are both empty.
        if self.taken_end == -1 and self.taken_start == -1:
            self.taken_end = cur_pos
            self.taken_start = cur_pos
            self.taken_previous[cur_pos] = -1
        else:
            last_taken_pos: int = self.taken_end
            self.taken_next[last_taken_pos] = cur_pos
            self.taken_end = cur_pos
            self.taken_previous[cur_pos] = last_taken_pos

        self.taken_next[cur_pos] = -1

        return cur_pos

    def remove(self, cs_pos: int):
        """
        Remove a choice structure position from the Priority object

        Parameters
        ----------
        cs_pos : int
            Position of the choice structure
        """
        t_prev: int = self.taken_previous[cs_pos]
        t_next: int = self.taken_next[cs_pos]
        last_empty_pos: int = self.empty_end

        # Update cs/median_indexes and empty_end/next
        self.cs_indexes[cs_pos] = -1
        self.median_indexes[cs_pos] = -1
        self.empty_next[last_empty_pos] = cs_pos
        self.empty_next[cs_pos] = -1
        self.empty_end = cs_pos

        # Remove (set to -1) various "taken" values or update them based on whether or not
        # the taken_previous/next values at the given choice structure position are empty.
        if t_prev == -1 and t_next == -1:
            self.taken_start = -1
            self.taken_end = -1
        elif t_prev == -1:
            self.taken_start = self.taken_next[cs_pos]
            self.taken_previous[self.taken_start] = -1
        elif t_next == -1:
            self.taken_end = self.taken_previous[cs_pos]
            self.taken_next[self.taken_end] = -1
        else:
            self.taken_previous[t_next] = t_prev
            self.taken_next[t_prev] = t_next
