"""
Priority for use in the PathGroups Algorithm

Based on Priority.java from C. Zheng and D. Sankoff (2011)

Author: Holger Jensen, Oskar Jensen
"""
from typing import List, Optional


class Priority:
    """
    Based on the priority system described in 'On the PATHGROUPS approach to rapid small phylogeny'

    Attributes
    ----------
    cn : int
        Current cycle
    bcla : int
        Best cycle look ahead
    bw : int
        Better or worse heuristic
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
    empty_previous : List[int]
        Previous empty values
    empty_next : List[int]
        Next empty values
    taken_previous : List[int]
        Previous taken values
    taken_next : List[int]
        Next taken values

    """

    def __init__(self,
                 cycle_now: int,
                 best_cycle_look_ahead: int,
                 better_or_worse: int,
                 best_cycle_look_ahead2: Optional[int],
                 size: int):
        """
        Constructor

        Parameters
        ----------
        cycle_now
            Current cycle
        best_cycle_look_ahead
            Best cycle look ahead
        better_or_worse
            Better or worse heuristic
        best_cycle_look_ahead2
            For certain operations
        size
            Size of each array in this class
        """
        self.bcla: int

        self.cn: int = cycle_now

        if best_cycle_look_ahead2 is not None:
            self.bcla = best_cycle_look_ahead2
        else:
            self.bcla = best_cycle_look_ahead

        self.bw: int = better_or_worse

        self.cs_indexes: List[int] = list()
        self.median_indexes: List[int] = list()
        self.empty_previous: List[int] = list()
        self.empty_next: List[int] = list()
        self.taken_previous: List[int] = list()
        self.taken_next: List[int] = list()

        # Initializes all the arrays to various default values
        for i in range(size):
            self.taken_previous.append(0)  # These two arrays are filled in so
            self.taken_next.append(0)      # that their sizes match the others
            self.cs_indexes.append(-1)
            self.median_indexes.append(-1)
            self.empty_previous.append(i - 1)
            self.empty_next.append(i + 1)

        self.empty_next[size - 1] = -1     # last entry set to -1

        self.empty_start: int = 0
        self.empty_end: int = size - 1
        self.taken_start: int = -1
        self.taken_end: int = -1

    def insert(self, cs_index: int, which_ancestor: Optional[int]) -> int:
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
        self.empty_previous[first_empty_pos] = -1

        # Update various taken values based on whether or not
        # the taken_start/end values at the current empty start position are both empty.
        if self.taken_end == -1 and self.taken_start == -1:
            self.taken_end = self.empty_start
            self.taken_start = self.empty_start
            self.taken_previous[self.empty_start] = -1
        else:
            last_taken_pos: int = self.taken_end
            self.taken_next[last_taken_pos] = cur_pos
            self.taken_end = cur_pos
            self.taken_previous[cur_pos] = last_taken_pos
        self.taken_next[cur_pos] = -1

        return self.empty_start

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
