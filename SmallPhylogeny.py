from copy import deepcopy, copy
from typing import Optional, List

from ChoiceStructure import ChoiceStructure
from PGMFragment import PGMFragment
from PGMFragment import combine
from PGMPath import PGMPath
from Priority import Priority
from TreeStructure import TreeStructure


def check_temp_list(temp: List[ChoiceStructure], ancestor_median: int, index_from: int) -> int:
    """
    Checks for valid genome and head value

    Parameters
    ----------
    temp
        Where to check
    ancestor_median
        The ancestor median
    index_from
        The head value

    Returns
    -------
    Index of the first valid ChoiceStructure in temp
    """
    for i in range(len(temp)):
        if temp[i] is not None and temp[i].for_which_genome == ancestor_median and temp[i].index_from == index_from:
            return i

    return -1


def count_cycle_look_ahead(ancestor_choice_structure: Optional[ChoiceStructure]) -> int:
    """
    Counts cycle look ahead in not-2-steps

    Parameters
    ----------
    ancestor_choice_structure
        The ancestor choice structure

    Returns
    -------
    Count of cycle look aheads
    """
    if ancestor_choice_structure.genome_1_path.tail < 0 and \
            ancestor_choice_structure.genome_2_path.tail < 0 and \
            ancestor_choice_structure.genome_3_path.tail < 0:
        return 3
    elif ancestor_choice_structure.genome_1_path.tail < 0 and ancestor_choice_structure.genome_2_path.tail < 0:
        return 2
    elif ancestor_choice_structure.genome_1_path.tail < 0 and ancestor_choice_structure.genome_3_path.tail < 0:
        return 2
    elif ancestor_choice_structure.genome_2_path.tail < 0 and ancestor_choice_structure.genome_3_path.tail < 0:
        return 2
    elif (
            ancestor_choice_structure.genome_1_path.genome_head ==
            ancestor_choice_structure.genome_1_path.genome_tail and
            ancestor_choice_structure.genome_2_path.genome_head ==
            ancestor_choice_structure.genome_2_path.genome_tail and
            ancestor_choice_structure.genome_3_path.genome_head ==
            ancestor_choice_structure.genome_3_path.genome_tail
    ) and (
            ancestor_choice_structure.genome_1_path.tail ==
            ancestor_choice_structure.genome_2_path.tail and
            ancestor_choice_structure.genome_1_path.tail ==
            ancestor_choice_structure.genome_3_path.tail):
        return 3
    elif (
            ancestor_choice_structure.genome_1_path.genome_head ==
            ancestor_choice_structure.genome_1_path.genome_tail and
            ancestor_choice_structure.genome_2_path.genome_head ==
            ancestor_choice_structure.genome_2_path.genome_tail
    ) and (
            ancestor_choice_structure.genome_1_path.tail ==
            ancestor_choice_structure.genome_2_path.tail):
        return 2
    elif (
            ancestor_choice_structure.genome_1_path.genome_head ==
            ancestor_choice_structure.genome_1_path.genome_tail and
            ancestor_choice_structure.genome_3_path.genome_head ==
            ancestor_choice_structure.genome_3_path.genome_tail
    ) and (
            ancestor_choice_structure.genome_1_path.tail ==
            ancestor_choice_structure.genome_3_path.tail):
        return 2
    elif (
            ancestor_choice_structure.genome_2_path.genome_head ==
            ancestor_choice_structure.genome_2_path.genome_tail and
            ancestor_choice_structure.genome_3_path.genome_head ==
            ancestor_choice_structure.genome_3_path.genome_tail
    ) and (
            ancestor_choice_structure.genome_2_path.tail ==
            ancestor_choice_structure.genome_3_path.tail):
        return 2

    return 1


def is_a_cycle(path1: PGMPath, path2: PGMPath) -> bool:
    return path1.head == path2.tail and \
           path1.tail == path2.head and \
           path1.genome_head == path2.genome_tail and \
           path1.genome_tail == path2.genome_head


class SmallPhylogeny:
    """
    Attributes
    ----------
    tree : TreeStructure
        TreeStructure
    priorities : List[Priority]
        List of priorities
    to_replace : int
        Which to replace
    """

    def __init__(self, tree_structure: TreeStructure):
        """
        Constructor

        Parameters
        ----------
        tree_structure
            TreeStructure
        """
        self.tree: TreeStructure = tree_structure
        self.priorities: List[Priority] = list()
        self.to_replace: int = 2

        priority_size: int = (self.tree.number_of_ancestors * (self.tree.gene_number + 500)) * 2
        self.priorities.append(Priority(3, 0, 0, None, priority_size))

        for cn in range(2, 0, -1):
            for bcla in range(3, 0, -1):
                for bw in range(4, -5, -1):
                    for bcla2 in range(3, 0, -1):
                        self.priorities.append(Priority(cn, bcla, bw, bcla2, priority_size))

    def group_pathgroup_into_priorities(self):
        """
        Groups path group into the instance priorities list
        """
        for j in range(len(self.tree.medians[0].choice_structures)):
            for i in range(len(self.tree.medians)):
                if self.tree.medians[i].choice_structures[j] is not None:
                    priority_count: int = self.get_priority_count(i, j)

                    if priority_count < len(self.priorities):
                        insert_position: int = self.priorities[priority_count].insert(j, i)
                        self.tree.medians[i].choice_structures[j].priority = priority_count
                        self.tree.medians[i].choice_structures[j].position = insert_position

    def get_priority_count(self, median_index: int, choice_structure_index: int) -> int:
        """
        Takes the result of calculate_case and checks it to get the priority count

        Parameters
        ----------
        median_index
            Which median to use
        choice_structure_index
            Which choice structure within the median, to use

        Returns
        -------
        int
            The priority count
        """
        result: int = 200
        index_from: int = self.tree.medians[median_index].choice_structures[choice_structure_index].index_from
        which_genome: int = self.tree.medians[median_index].choice_structures[choice_structure_index].for_which_genome
        tail_of1: int = self.tree.medians[median_index].choice_structures[choice_structure_index].genome_1_path.tail
        tail_of2: int = self.tree.medians[median_index].choice_structures[choice_structure_index].genome_2_path.tail
        tail_of3: int = self.tree.medians[median_index].choice_structures[choice_structure_index].genome_3_path.tail
        ancestor_priority: int
        add_tail1: bool = False
        add_tail2: bool = False

        if self.tree.medians[median_index].choice_structures[choice_structure_index].genome_1_path.genome_head == \
                self.tree.medians[median_index].choice_structures[choice_structure_index].genome_1_path.genome_tail:
            ancestor_priority = self.calculate_case(median_index, choice_structure_index, index_from, tail_of1)
            add_tail1 = True

            if 0 <= ancestor_priority <= 1:
                self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                    index_from, tail_of1, which_genome, which_genome)

                return ancestor_priority
            if ancestor_priority < result:
                result = ancestor_priority
                self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                    index_from, tail_of1, which_genome, which_genome)
            elif ancestor_priority == result and self.to_replace == 1:
                result = ancestor_priority
                self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                    index_from, tail_of1, which_genome, which_genome)

        self.to_replace = 3 - self.to_replace

        if self.tree.medians[median_index].choice_structures[choice_structure_index].genome_2_path.genome_head == \
                self.tree.medians[median_index].choice_structures[choice_structure_index].genome_2_path.genome_tail:
            if (not add_tail1) or tail_of2 != tail_of1:
                ancestor_priority = self.calculate_case(median_index, choice_structure_index, index_from, tail_of2)

                if ancestor_priority == 1:
                    self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                        index_from, tail_of2, which_genome, which_genome)

                    return ancestor_priority
                if ancestor_priority < result:
                    result = ancestor_priority
                    self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                        index_from, tail_of2, which_genome, which_genome)
                elif ancestor_priority == result and self.to_replace == 1:
                    result = ancestor_priority
                    self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                        index_from, tail_of2, which_genome, which_genome)

        self.to_replace = 3 - self.to_replace

        if self.tree.medians[median_index].choice_structures[choice_structure_index].genome_3_path.genome_head == \
                self.tree.medians[median_index].choice_structures[choice_structure_index].genome_3_path.genome_tail:
            if ((not add_tail1) or tail_of3 != tail_of1) and ((not add_tail2) or tail_of3 != tail_of2):
                ancestor_priority = self.calculate_case(median_index, choice_structure_index, index_from, tail_of3)

                if ancestor_priority < result:
                    result = ancestor_priority
                    self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                        index_from, tail_of3, which_genome, which_genome)
                elif ancestor_priority == result and self.to_replace == 1:
                    result = ancestor_priority
                    self.tree.medians[median_index].choice_structures[choice_structure_index].gray_edge = PGMPath(
                        index_from, tail_of3, which_genome, which_genome)

        self.to_replace = 3 - self.to_replace

        return result

    def calculate_case(self, median_index: int, choice_structure_index: int, index_from: int, tail: int) -> int:
        """
        Figures out which priority case to use (combination of bcla, bw, cn)

        Parameters
        ----------
        median_index
            Index of the current median
        choice_structure_index
            Index of the ChoiceStructure
        index_from
            Index of the starting node
        tail
            Index of the tail node

        Returns
        -------
        The index into self.priorities
        """
        if tail < 0:
            return 200

        cycle_now: int = self.count_cycle_for_edge(median_index, choice_structure_index, tail)

        if cycle_now == 3:
            return 0

        created_fragment: List[PGMFragment] = self.get_created_fragment(median_index, index_from, tail)
        self.get_new_fragment_list(created_fragment, median_index)
        created_choice_structure: List[ChoiceStructure] = self.get_new_choice_structure(median_index,
                                                                                        choice_structure_index,
                                                                                        index_from,
                                                                                        tail)

        look_ahead_cycles: List[int] = self.count_all_look_ahead_cycles(created_choice_structure)
        max_cycle_look_ahead: int = look_ahead_cycles[0]
        how_many_more_cycles: int = look_ahead_cycles[1]
        max_cycle_look_ahead2: int = look_ahead_cycles[2]

        if len(created_fragment) != 0:  # recoverOriginalFragment()
            self.tree.medians[median_index].fragments[created_fragment[2].end1] = None
            self.tree.medians[median_index].fragments[created_fragment[2].end2] = None
            self.tree.medians[median_index].fragments[created_fragment[0].end1] = created_fragment[0]
            self.tree.medians[median_index].fragments[created_fragment[0].end2] = PGMFragment(created_fragment[0].end2,
                                                                                              created_fragment[0].end1)
            self.tree.medians[median_index].fragments[created_fragment[1].end1] = created_fragment[1]
            self.tree.medians[median_index].fragments[created_fragment[1].end2] = PGMFragment(created_fragment[1].end2,
                                                                                              created_fragment[1].end1)

        index: int = 1

        for cn in range(2, 0, -1):
            for bcla in range(3, 0, -1):
                for bw in range(4, -5, -1):
                    for bcla2 in range(3, 0, -1):
                        if cycle_now == cn and max_cycle_look_ahead == bcla and \
                                how_many_more_cycles == bw and max_cycle_look_ahead2 == bcla2:
                            return index

                        index += 1

        return 200

    def count_cycle_for_edge(self, median_index: int, choice_structure_index: int, tail: int) -> int:
        """
        Counts the cycles for the current edge

        Parameters
        ----------
        median_index
            Index of the current median
        choice_structure_index
            Index of the ChoiceStructure
        tail
            Tail node

        Returns
        -------
        int
            Number of cycles
        """
        is_circular_chromosome: bool = False
        head: int = self.tree.medians[median_index].choice_structures[choice_structure_index].index_from

        if 0 < head == self.tree.medians[median_index].fragments[head].end1 and \
                self.tree.medians[median_index].fragments[head] is not None and \
                self.tree.medians[median_index].fragments[head].end2 == tail:
            is_circular_chromosome = True
        elif 0 < tail == self.tree.medians[median_index].fragments[tail].end1 and \
                self.tree.medians[median_index].fragments[tail] is not None and \
                self.tree.medians[median_index].fragments[tail].end2 == head:
            is_circular_chromosome = True

        if is_circular_chromosome:
            return -1

        result: int = 0

        if tail > 0:
            ancestor_choice_structure: ChoiceStructure = self.tree.medians[median_index].choice_structures[
                choice_structure_index]

            if ancestor_choice_structure.genome_1_path.genome_head == \
                    ancestor_choice_structure.genome_1_path.genome_tail and \
                    ancestor_choice_structure.genome_1_path.tail == tail:
                result += 1

            if ancestor_choice_structure.genome_2_path.genome_head == \
                    ancestor_choice_structure.genome_2_path.genome_tail and \
                    ancestor_choice_structure.genome_2_path.tail == tail:
                result += 1

            if ancestor_choice_structure.genome_3_path.genome_head == \
                    ancestor_choice_structure.genome_3_path.genome_tail and \
                    ancestor_choice_structure.genome_3_path.tail == tail:
                result += 1

        return result

    def get_created_fragment(self, median_index: int, node1: int, node2: int) -> List[PGMFragment]:
        """
        Returns PGMFragments based on given nodes, as well as the ancestor fragment of both of them

        Parameters
        ----------
        median_index
            Index of current median
        node1
            First node to connect
        node2
            Second node to connect

        Returns
        -------
        [PGMFragment]
            Both node fragments + their ancestor fragment
        """
        result: List[PGMFragment] = list()

        if node1 > 0 and node2 > 0:
            result.append(PGMFragment.from_fragment(self.tree.medians[median_index].fragments[node1]))
            result.append(PGMFragment.from_fragment(self.tree.medians[median_index].fragments[node2]))

            ancestor_line: PGMPath = PGMPath(node1, node2, 1, 1)

            if self.tree.medians[median_index].fragments[node1] is not None:
                ancestor_fragment: PGMFragment = combine(
                    ancestor_line,
                    self.tree.medians[median_index].fragments[node1],
                    self.tree.medians[median_index].fragments[node2])

                result.append(ancestor_fragment)

        return result

    def get_new_fragment_list(self, created_fragment: List[PGMFragment], median_index: int):
        """
        Sets the instance fragment list for performing small phylogeny operations

        Parameters
        ----------
        created_fragment
            New fragment to be operated on
        median_index
            Index of current median
        """
        if len(created_fragment) != 0:
            self.tree.medians[median_index].fragments[created_fragment[0].end1] = None
            self.tree.medians[median_index].fragments[created_fragment[0].end2] = None
            self.tree.medians[median_index].fragments[created_fragment[1].end1] = None
            self.tree.medians[median_index].fragments[created_fragment[1].end2] = None
            self.tree.medians[median_index].fragments[created_fragment[2].end1] = created_fragment[2]
            self.tree.medians[median_index].fragments[created_fragment[2].end2] = PGMFragment(created_fragment[2].end2,
                                                                                              created_fragment[2].end1)

    def count_all_look_ahead_cycles(self, created_choice_structures: List[ChoiceStructure]) -> List[int]:
        """
        Counts all look ahead cycles for the set of ChoiceStructures

        Parameters
        ----------
        created_choice_structures
            List of ChoiceStructures to count for

        Returns
        -------
        [int]
            Counts of look ahead cycles for each ChoiceStructure
        """
        result: List[int] = [int(), int(), int()]
        max_cycle: int = 1
        total: int = 0

        for choice_structure in created_choice_structures:
            if choice_structure is not None:
                index_from: int = choice_structure.index_from
                ancestor: int = choice_structure.for_which_genome - self.tree.number_of_leaves
                old_cycle: int = count_cycle_look_ahead(
                    self.tree.medians[ancestor].choice_structures[index_from - 1])
                new_cycle: int = count_cycle_look_ahead(choice_structure)

                if new_cycle > old_cycle and new_cycle > max_cycle:
                    max_cycle = new_cycle

                total += new_cycle - old_cycle

        result[0] = max_cycle
        result[1] = total
        result[2] = self.count_2_steps_look_ahead_cycle(created_choice_structures, max_cycle)

        return result

    def count_2_steps_look_ahead_cycle(self, created_choice_structures: List[ChoiceStructure], max_cycle: int) -> int:
        """
        Counts look ahead cycles in 2 steps

        Parameters
        ----------
        created_choice_structures
            The created choice structures
        max_cycle
            How high the cycles can go

        Returns
        -------
        Number of look ahead cycles
        """
        current_max: int = 1

        for ancestor_choice_structure in created_choice_structures:
            current_count: int = count_cycle_look_ahead(ancestor_choice_structure)
            new_count: int

            if current_count == max_cycle:
                index_from: int = ancestor_choice_structure.index_from
                current_genome: int = ancestor_choice_structure.for_which_genome

                tail1: int = ancestor_choice_structure.genome_1_path.tail
                genome_tail1: int = ancestor_choice_structure.genome_1_path.genome_tail
                genome_head1: int = ancestor_choice_structure.genome_1_path.genome_head

                tail2: int = ancestor_choice_structure.genome_2_path.tail
                genome_tail2: int = ancestor_choice_structure.genome_2_path.genome_tail
                genome_head2: int = ancestor_choice_structure.genome_2_path.genome_head

                tail3: int = ancestor_choice_structure.genome_3_path.tail
                genome_tail3: int = ancestor_choice_structure.genome_3_path.genome_tail
                genome_head3: int = ancestor_choice_structure.genome_3_path.genome_head

                if max_cycle > 1:
                    current_tail: int = -1000

                    if tail1 == tail2 or tail1 == tail3:
                        current_tail = tail1
                    if tail2 == tail3:
                        current_tail = tail2

                    current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                        index_from, current_tail, current_genome, ancestor_choice_structure, created_choice_structures)

                    for new_choice_structure in current_new_choice_structure:
                        new_count = count_cycle_look_ahead(new_choice_structure)

                        if new_count == 3:
                            return new_count
                        elif current_max < new_count:
                            current_max = new_count

                else:
                    if genome_tail1 == genome_head1:
                        current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                            index_from, tail1, current_genome, ancestor_choice_structure, created_choice_structures)

                        for new_choice_structure in current_new_choice_structure:
                            new_count = count_cycle_look_ahead(new_choice_structure)

                            if new_count == 3:
                                return new_count
                            elif current_max < new_count:
                                current_max = new_count
                    if genome_tail2 == genome_head2:
                        current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                            index_from, tail2, current_genome, ancestor_choice_structure, created_choice_structures)

                        for new_choice_structure in current_new_choice_structure:
                            new_count = count_cycle_look_ahead(new_choice_structure)

                            if new_count == 3:
                                return new_count
                            elif current_max < new_count:
                                current_max = new_count
                    if genome_tail3 == genome_head3:
                        current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                            index_from, tail3, current_genome, ancestor_choice_structure, created_choice_structures)

                        for new_choice_structure in current_new_choice_structure:
                            new_count = count_cycle_look_ahead(new_choice_structure)

                            if new_count == 3:
                                return new_count
                            elif current_max < new_count:
                                current_max = new_count

        return current_max

    def get_new_choice_structure(self, median_index: int, choice_structure_index: int, index_from: int, tail: int) -> \
            List[Optional[ChoiceStructure]]:
        """
        Gets new ChoiceStructure in some time that isn't 2 steps

        Parameters
        ----------
        median_index
            Index of the current median genome
        choice_structure_index
            Index of the current ChoiceStructure
        index_from
            Head index
        tail
            Tail index

        Returns
        -------
        New ChoiceStructure
        """
        path_1_1: PGMPath = \
            copy(self.tree.medians[median_index].choice_structures[choice_structure_index].genome_1_path)
        path_1_2: PGMPath = copy(self.tree.medians[median_index].choice_structures[tail - 1].genome_1_path)
        path_2_1: PGMPath = \
            copy(self.tree.medians[median_index].choice_structures[choice_structure_index].genome_2_path)
        path_2_2: PGMPath = copy(self.tree.medians[median_index].choice_structures[tail - 1].genome_2_path)
        path_3_1: PGMPath = \
            copy(self.tree.medians[median_index].choice_structures[choice_structure_index].genome_3_path)
        path_3_2: PGMPath = copy(self.tree.medians[median_index].choice_structures[tail - 1].genome_3_path)
        path_l1: PGMPath = PGMPath(index_from, tail, path_1_1.genome_head, path_1_1.genome_head)
        path_l2: PGMPath = PGMPath(index_from, tail, path_2_1.genome_head, path_2_1.genome_head)
        path_l3: PGMPath = PGMPath(index_from, tail, path_3_1.genome_head, path_3_1.genome_head)

        for_which_genome: int = self.tree.medians[median_index].choice_structures[
            choice_structure_index].for_which_genome

        new_path1: PGMPath = path_1_1.connect(path_1_1, path_1_2, path_l1, for_which_genome)
        new_path2: PGMPath = path_2_1.connect(path_2_1, path_2_2, path_l2, for_which_genome)
        new_path3: PGMPath = path_3_1.connect(path_3_1, path_3_2, path_l3, for_which_genome)

        path_1_genome: int = self.tree.leaves[median_index + self.tree.number_of_leaves][0]
        path_2_genome: int = self.tree.leaves[median_index + self.tree.number_of_leaves][1]
        path_3_genome: int = self.tree.leaves[median_index + self.tree.number_of_leaves][2]

        temp: List[Optional[ChoiceStructure]] = [None for _ in range(12)]

        if not is_a_cycle(path_1_1, path_1_2):
            temp = self.get_new_choice_structure_based_on_path(new_path1, temp, path_1_genome, median_index, None)
        if not is_a_cycle(path_2_1, path_2_2):
            temp = self.get_new_choice_structure_based_on_path(new_path2, temp, path_2_genome, median_index, None)
        if not is_a_cycle(path_3_1, path_3_2):
            temp = self.get_new_choice_structure_based_on_path(new_path3, temp, path_3_genome, median_index, None)

        new_choice_structure_number: int = 0

        for choice_structure in temp:
            if choice_structure is not None:
                new_choice_structure_number += 1

        result: List[Optional[ChoiceStructure]] = [None for _ in range(new_choice_structure_number)]

        if len(result) >= 0:
            result = deepcopy(temp)[:len(result)]

        return result

    def get_new_choice_structure_2_step(self, index_from: int, tail: int, current_genome: int,
                                        ancestor_choice_structure: ChoiceStructure,
                                        new_choice_structures: List[ChoiceStructure]) -> List[ChoiceStructure]:
        """
        Gets a new set of ChoiceStructures in 2 steps

        Parameters
        ----------
        index_from
            The head node
        tail
            Tail node
        current_genome
            Current genome operating on
        ancestor_choice_structure
            Choice structure for ancestor node
        new_choice_structures
            New ChoiceStructures to use

        Returns
        -------
        Set of new ChoiceStructures
        """
        if tail <= 0:
            return list()

        path_1_1: PGMPath = ancestor_choice_structure.genome_1_path
        path_2_1: PGMPath = ancestor_choice_structure.genome_2_path
        path_3_1: PGMPath = ancestor_choice_structure.genome_3_path

        new_choice_structure: Optional[ChoiceStructure] = None
        cont: bool = False

        for choice_structure in new_choice_structures:  # getTheCSWithGandStart(t, ghere, allncs)
            if choice_structure.index_from == tail and choice_structure.for_which_genome == current_genome:
                new_choice_structure = ChoiceStructure()
                new_choice_structure.from_cs(choice_structure)
                cont = True
                break

        if not cont:
            new_choice_structure = ChoiceStructure()
            new_choice_structure.from_cs(
                self.tree.medians[current_genome - self.tree.number_of_leaves].choice_structures[tail - 1])

        path_1_2: PGMPath = new_choice_structure.genome_1_path
        path_2_2: PGMPath = new_choice_structure.genome_2_path
        path_3_2: PGMPath = new_choice_structure.genome_3_path

        path_l1: PGMPath = PGMPath(index_from, tail, path_1_1.genome_head, path_1_1.genome_head)
        path_l2: PGMPath = PGMPath(index_from, tail, path_2_1.genome_head, path_2_1.genome_head)
        path_l3: PGMPath = PGMPath(index_from, tail, path_3_1.genome_head, path_3_1.genome_head)

        for_which_genome: int = ancestor_choice_structure.for_which_genome

        new_path1: PGMPath = path_1_1.connect(path_1_1, path_1_2, path_l1, for_which_genome)
        new_path2: PGMPath = path_2_1.connect(path_2_1, path_2_2, path_l2, for_which_genome)
        new_path3: PGMPath = path_3_1.connect(path_3_1, path_3_2, path_l3, for_which_genome)

        median_index: int = ancestor_choice_structure.for_which_genome - self.tree.number_of_leaves

        path_1_genome: int = self.tree.leaves[ancestor_choice_structure.for_which_genome][0]
        path_2_genome: int = self.tree.leaves[ancestor_choice_structure.for_which_genome][1]
        path_3_genome: int = self.tree.leaves[ancestor_choice_structure.for_which_genome][2]

        temp: List[Optional[ChoiceStructure]] = [None for _ in range(12)]

        if not is_a_cycle(path_1_1, path_1_2):
            temp = self.get_new_choice_structure_based_on_path(
                new_path1, temp, path_1_genome, median_index, new_choice_structures)
        if not is_a_cycle(path_2_1, path_2_2):
            temp = self.get_new_choice_structure_based_on_path(
                new_path2, temp, path_2_genome, median_index, new_choice_structures)
        if not is_a_cycle(path_3_1, path_3_2):
            temp = self.get_new_choice_structure_based_on_path(
                new_path3, temp, path_3_genome, median_index, new_choice_structures)

        new_choice_structure_number: int = 0

        for choice_structure in temp:
            if choice_structure is not None:
                new_choice_structure_number += 1

        result: List[Optional[ChoiceStructure]] = [None for _ in range(new_choice_structure_number)]

        if len(result) >= 0:
            result = deepcopy(temp)[:len(result)]

        return result

    def get_new_choice_structure_based_on_path(self, new_path1: PGMPath,
                                               temp: List[ChoiceStructure],
                                               for_which_genome: int,
                                               median_index: int,
                                               new_choice_structures: Optional[List[ChoiceStructure]]) \
            -> List[ChoiceStructure]:
        """
        Gets a new ChoiceStructure based on given PGMPath

        Parameters
        ----------
        new_path1
            Path to base ChoiceStructure on
        temp
            Modified to add the new ChoiceStructure to
        for_which_genome
            Which genome this is for
        median_index
            Index of the median
        new_choice_structures
            For certain circumstances

        Returns
        -------
        temp, modified with the new ChoiceStructure
        """
        index_from: int = new_path1.head
        genome_from: int = new_path1.genome_head
        from_tail: int = new_path1.tail
        genome_tail: int = new_path1.genome_tail
        temp_index: int = 0

        for choice_structure in temp:
            if choice_structure is not None:
                temp_index += 1

        if index_from > 0:
            ancestor_median: int = -1

            if genome_from == genome_tail and genome_from == for_which_genome:
                ancestor_median = median_index + self.tree.number_of_leaves
            elif genome_from == genome_tail and genome_from != for_which_genome:
                ancestor_median = for_which_genome
            elif genome_from != genome_tail:
                if genome_from == for_which_genome:
                    ancestor_median = median_index + self.tree.number_of_leaves
                elif genome_from != for_which_genome:
                    ancestor_median = for_which_genome

            if ancestor_median != -1:
                index_in_temp: int = check_temp_list(temp, ancestor_median, index_from)

                if index_in_temp == -1:
                    if new_choice_structures is None:
                        if self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                index_from - 1] is not None:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(self.tree.medians[ancestor_median - self.tree.number_of_leaves].
                                                     choice_structures[index_from - 1])
                            temp[temp_index].set_new_path(new_path1)
                            temp_index += 1
                    else:  # getNewCsBasedOnAPathLA2() case
                        find_now: bool = False

                        for choice_structure in new_choice_structures:
                            if choice_structure.index_from == index_from and \
                                    choice_structure.for_which_genome == ancestor_median:
                                temp[temp_index] = ChoiceStructure()
                                temp[temp_index].from_cs(choice_structure)
                                temp[temp_index].set_new_path(new_path1)
                                find_now = True
                                break

                        if (not find_now) and \
                                self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                    index_from - 1] is not None:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(
                                self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                    index_from - 1])
                            temp[temp_index].set_new_path(new_path1)
                            temp_index += 1
                else:
                    temp[index_in_temp].set_new_path(new_path1)

        if from_tail > 0:
            new_path2: PGMPath = PGMPath(from_tail, index_from, genome_tail, genome_from)
            ancestor_median: int = -1

            if genome_from == genome_tail and genome_from == for_which_genome:
                ancestor_median = median_index + self.tree.number_of_leaves
            elif genome_from == genome_tail and genome_from != for_which_genome:
                ancestor_median = for_which_genome
            elif genome_from != genome_tail:
                if genome_tail == for_which_genome:
                    ancestor_median = median_index + self.tree.number_of_leaves
                elif genome_tail != for_which_genome:
                    ancestor_median = for_which_genome

            if ancestor_median != -1:
                index_in_temp: int = check_temp_list(temp, ancestor_median, from_tail)

                if index_in_temp == -1:
                    if new_choice_structures is None:
                        if self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                from_tail - 1] is not None:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(
                                self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                    from_tail - 1])
                            temp[temp_index].set_new_path(new_path2)
                            temp_index += 1
                    else:  # getNewCsBasedOnAPathLA2() case
                        for choice_structure in new_choice_structures:
                            if choice_structure.index_from == from_tail and \
                                    choice_structure.for_which_genome == ancestor_median:
                                temp[temp_index] = ChoiceStructure()
                                temp[temp_index].from_cs(choice_structure)
                                temp[temp_index].set_new_path(new_path2)
                                break

                        if self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                from_tail - 1] is not None:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(
                                self.tree.medians[ancestor_median - self.tree.number_of_leaves].choice_structures[
                                    from_tail - 1])
                            temp[temp_index].set_new_path(new_path2)
                            temp_index += 1
                else:
                    temp[index_in_temp].set_new_path(new_path2)

        return temp

    def get_result(self):
        """
        Gets the result of the small phylogeny operation
        """
        self.group_pathgroup_into_priorities()

        best_choice_structure: List[int] = self.find_the_best_choice_structure()

        while best_choice_structure[0] != -1:
            priority: int = best_choice_structure[0]
            for_which_priority: int = best_choice_structure[1]
            current_ancestor: int = self.priorities[priority].median_indexes[for_which_priority]
            current_choice_structure_index: int = self.priorities[priority].cs_indexes[for_which_priority]

            self.add_gray_edge(current_ancestor, self.tree.medians[current_ancestor].choice_structures[
                current_choice_structure_index].gray_edge)
            best_choice_structure = self.find_the_best_choice_structure()

    def add_gray_edge(self, ancestor: int, gray_edge: PGMPath):
        """
        Adds a gray edge object and updates self

        Parameters
        ----------
        ancestor
            Which ancestor to use
        gray_edge
            Temporary variable used in calling function
        """
        if gray_edge.head >= 0 and gray_edge.tail >= 0:
            self.tree.medians[ancestor].gray_edge[self.tree.medians[ancestor].gray_edge_index] = gray_edge
            self.tree.medians[ancestor].gray_edge_index += 1

        self.update_all(ancestor, gray_edge.head - 1, gray_edge)

    def update_all(self, median_index: int, choice_structure_index: int, path_l: PGMPath):
        """
        Updates internal data

        Parameters
        ----------
        median_index
            Index of the current median candidate
        choice_structure_index
            Index of the current ChoiceStructure
        path_l
            Used for combining paths
        """
        new_fragment: List[PGMFragment] = self.get_created_fragment(median_index, path_l.head, path_l.tail)
        self.get_new_fragment_list(new_fragment, median_index)

        self.priorities[self.tree.medians[median_index].choice_structures[path_l.head - 1].priority].remove(
            self.tree.medians[median_index].choice_structures[path_l.head - 1].position)
        self.priorities[self.tree.medians[median_index].choice_structures[path_l.tail - 1].priority].remove(
            self.tree.medians[median_index].choice_structures[path_l.tail - 1].position)

        new_choice_structure: List[ChoiceStructure] = self.get_new_choice_structure(median_index,
                                                                                    choice_structure_index,
                                                                                    path_l.head, path_l.tail)

        for choice_structure in new_choice_structure:
            for_which_genome: int = choice_structure.for_which_genome
            index_from: int = choice_structure.index_from

            self.tree.medians[for_which_genome - self.tree.number_of_leaves].choice_structures[
                index_from - 1] = choice_structure

        self.tree.medians[median_index].choice_structures[path_l.head - 1] = None
        self.tree.medians[median_index].choice_structures[path_l.tail - 1] = None

        for choice_structure in new_choice_structure:
            current_ancestor: int = choice_structure.for_which_genome - self.tree.number_of_leaves
            current_choice_structure_index: int = choice_structure.index_from - 1
            priority: int = self.get_priority_count(current_ancestor, current_choice_structure_index)

            if priority != choice_structure.priority:
                old_priority: int = choice_structure.priority
                old_position: int = choice_structure.position

                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int = -1

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(choice_structure.index_from - 1, current_ancestor)

                self.tree.medians[current_ancestor].choice_structures[
                    choice_structure.index_from - 1].priority = priority
                self.tree.medians[current_ancestor].choice_structures[
                    choice_structure.index_from - 1].position = new_position

        if new_fragment[2].end1 > 0:
            self.update_priority(median_index, new_fragment[2].end1 - 1)

        if new_fragment[2].end2 > 0:
            self.update_priority(median_index, new_fragment[2].end2 - 1)

        for choice_structure in new_choice_structure:
            tail1: int = choice_structure.genome_1_path.tail
            tail2: int = choice_structure.genome_2_path.tail
            tail3: int = choice_structure.genome_3_path.tail
            which_ancestor: int = choice_structure.for_which_genome - self.tree.number_of_leaves

            if tail1 > 0:
                self.update_priority(which_ancestor, tail1 - 1)

            if tail2 > 0:
                self.update_priority(which_ancestor, tail2 - 1)

            if tail3 > 0:
                self.update_priority(which_ancestor, tail3 - 1)

    def update_priority(self, median_index: int, choice_structure_index: int):
        """
        Updates a priority for the given median

        Parameters
        ----------
        median_index
            Index of the current median candidate
        choice_structure_index
            Index of the current ChoiceStructure
        """
        if self.tree.medians[median_index].choice_structures[choice_structure_index] is not None:
            old_priority: int = self.tree.medians[median_index].choice_structures[choice_structure_index].priority
            old_position: int = self.tree.medians[median_index].choice_structures[choice_structure_index].position
            priority: int = self.get_priority_count(median_index, choice_structure_index)

            if priority != old_priority:
                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int = -1

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(choice_structure_index, median_index)

                self.tree.medians[median_index].choice_structures[choice_structure_index].priority = priority
                self.tree.medians[median_index].choice_structures[choice_structure_index].position = new_position

    def find_the_best_choice_structure(self) -> List[int]:
        """
        Finds the best choice structure

        Returns
        -------
        List[int]
            Result of the operation
        """
        result: List[int] = [int(), int()]
        result[0] = -1
        result[1] = -1

        for i in range(len(self.priorities)):
            if self.priorities[i].taken_start != -1:
                result[0] = i
                result[1] = self.priorities[i].taken_start

                return result

        return result
