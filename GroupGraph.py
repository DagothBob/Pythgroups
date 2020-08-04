from copy import copy, deepcopy
from random import Random
from typing import List, Optional

from ChoiceStructure import ChoiceStructure
from Genome import Genome
from PGMFragment import PGMFragment, combine
from PGMPath import PGMPath
from Priority import Priority


def insert_character(s: str, i: int, c: str) -> str:
    """
    Utility function for modifying strings in a similar way to lists.
    Character c is inserted at index i into string s and returned.

    Parameters
    ----------
    s
        String to replace character in
    i
        Index to replace in string
    c
        Character to insert

    Returns
    -------
    str
        New string with the character inserted
    """
    return s[:i] + c + s[(i + 1):]


def split_at_whitespace(strings: str) -> List[str]:
    """
    Strips, then splits, a string, then strips the substrings again

    Parameters
    ----------
    strings
        List of strings to operate on

    Returns
    -------
    [str]
        Set of cleaned-up strings
    """
    result: List[str] = list()

    for string in strings.strip().split(" "):
        if string.strip() != "":
            result.append(string.strip())

    return result


def count_gene_number(genome: Genome) -> int:
    result: int = 0

    for chromosome in genome.chromosomes:
        result += len(chromosome.genes)

    return result


def check_temp_list(temp_list: List[ChoiceStructure], f: int) -> int:
    for i in range(len(temp_list)):
        if temp_list[i] is not None:
            if temp_list[i].genome_1_path.head == f:
                return i

    return -1


def get_gene_next_node(index: int) -> int:
    if int((index / 2) * 2) == index:
        return index - 1

    return index + 1


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


class GroupGraph:
    def __init__(self, tetrad: Genome, outgroup: Genome, replace: int):
        self.node_int: List[int] = list()
        self.node_str: List[str] = list()
        self.ancestor_AA: List[str] = list()
        self.ancestor_A: List[str] = list()

        self.replace: int = replace
        self.gene_number: int = count_gene_number(outgroup)

        self.priorities: List[Priority] = self.initialize_priorities()
        self.set_nodes(outgroup)
        self.tetrad: List[PGMPath] = self.get_pgm_path(tetrad, 2)
        self.outgroup: List[PGMPath] = self.get_pgm_path(outgroup, 1)

        self.fragments: List[Optional[PGMFragment]] = [None for _ in range(self.gene_number * 2 + 1)]

        for i in range(int(len(self.fragments) / 2)):
            self.fragments[2 * i + 1] = PGMFragment(2 * i + 1, 2 * i + 2)
            self.fragments[2 * i + 2] = PGMFragment(2 * i + 2, 2 * i + 1)

        self.gray_edge: List[Optional[PGMPath]] = [None for _ in range(len(tetrad.chromosomes))]
        self.gray_edge_index: int = 0

        self.choice_structures: List[Optional[ChoiceStructure]] = [None for _ in range(self.gene_number * 2)]
        cs_index: int = 0

        for i in range(1, self.gene_number * 2 + 1):
            self.choice_structures[cs_index] = ChoiceStructure()
            self.choice_structures[cs_index].index_from = i
            self.choice_structures[cs_index].genome_1_path = copy(tetrad.chromosomes[i])
            self.choice_structures[cs_index].genome_2_path = copy(tetrad.chromosomes[i + self.gene_number * 2])
            self.choice_structures[cs_index].genome_3_path = copy(outgroup.chromosomes[i])
            self.choice_structures[cs_index].priority = 200
            self.choice_structures[cs_index].position = -1
            self.choice_structures[cs_index].gray_edge = None
            cs_index += 1

    def initialize_priorities(self) -> List[Priority]:
        priority_size: int = self.gene_number * 2 + 1000
        priorities: List[Priority] = list()

        priorities.append(Priority(3, 0, 0, None, priority_size))

        for cn in range(2, 0, -1):
            for bcla in range(3, 0, -1):
                for bw in range(4, -5, -1):
                    for bcla2 in range(3, 0, -1):
                        priorities.append(Priority(cn, bcla, bw, bcla2, priority_size))

        return priorities

    def set_nodes(self, genome: Genome):
        self.node_int = [int() for _ in range(self.gene_number * 2)]
        self.node_str = [str() for _ in range(self.gene_number * 2)]

        index: int = 0

        for chromosome in genome.chromosomes:
            for gene in chromosome.genes:
                node_1: str
                node_2: str

                if gene.name[0] == "-":
                    node_1 = gene.name[1:] + "h"
                    node_2 = gene.name[1:] + "t"

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_2
                    index += 1

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_1
                    index += 1
                else:
                    node_1 = gene.name + "t"
                    node_2 = gene.name + "h"

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_1
                    index += 1

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_2
                    index += 1

    def get_pgm_path(self, genome: Genome, ploidy: int) -> List[PGMPath]:
        """
        Gets list of PGMPaths for a genome

        Parameters
        ----------
        genome
            Genome to get the PGMPath for
        ploidy
            Monoploid or diploid

        Returns
        -------
        List[PGMPath]
            List of PGMPaths for the genome
        """
        path1: List[Optional[PGMPath]]

        if ploidy == 1:
            path1 = [None for _ in range((2 * self.gene_number) + 1)]
        else:
            path1 = [None for _ in range((4 * self.gene_number) + 1)]

        null_node: int = -1

        for chromosome in genome.chromosomes:
            pre_node: int = 0

            for j in range(len(chromosome.genes)):
                first_character: str = chromosome.genes[j].name[0]
                node1: str
                node2: str

                if first_character == '-':
                    node1 = insert_character(chromosome.genes[j].name[1:], len(chromosome.genes[j].name[1:]), "h")
                    node2 = insert_character(chromosome.genes[j].name[1:], len(chromosome.genes[j].name[1:]), "t")
                else:
                    node1 = insert_character(chromosome.genes[j].name, len(chromosome.genes[j].name), "t")
                    node2 = insert_character(chromosome.genes[j].name, len(chromosome.genes[j].name), "h")

                node1_int: int = self.find_node_int(node1, ploidy)
                node2_int: int = self.find_node_int(node2, ploidy)

                if node1_int == 0 or node2_int == 0:
                    print("Gene ", str(chromosome.genes[j]), " does not exist in the other genome.\n")

                if j == 0:
                    path1[node1_int] = PGMPath(node1_int, null_node, None, None)
                    pre_node = node2_int
                    null_node -= 1
                elif j != 0 and j != len(chromosome.genes) - 1:
                    path1[node1_int] = PGMPath(node1_int, pre_node, None, None)
                    path1[pre_node] = PGMPath(pre_node, node1_int, None, None)
                    pre_node = node2_int

                if j == len(chromosome.genes) - 1:
                    if len(chromosome.genes) != 1:
                        path1[pre_node] = PGMPath(pre_node, node1_int, None, None)
                        path1[node1_int] = PGMPath(node1_int, pre_node, None, None)

                    path1[node2_int] = PGMPath(node2_int, null_node, None, None)
                    null_node -= 1

        return path1

    def find_node_int(self, ancestor_string: str, ploidy: int) -> int:
        if ploidy == 1:
            for i in range(len(self.node_str)):
                if self.node_str[i] == ancestor_string:
                    return i + 1
        elif ploidy == 2:
            copy_: str = ancestor_string[-2]
            name: str = ancestor_string[:-2] + ancestor_string[-1]

            for i in range(len(self.node_str)):
                if self.node_str[i] == name:
                    if copy_ == "a":
                        return i + 1
                    elif copy_ == "b":
                        return (i + 1) + self.gene_number * 2

        return 0

    def get_result(self):
        self.group_pathgroup_into_priorities()

        best_cs: List[int] = self.find_the_best_choice_structure()

        while best_cs[0] != -1:
            self.add_gray_edge(self.choice_structures[self.priorities[best_cs[0]].cs_indexes[best_cs[1]]].gray_edge)
            best_cs = self.find_the_best_choice_structure()

        self.get_ancestors()

    def group_pathgroup_into_priorities(self):
        """
        Groups path group into the instance priorities list
        """
        for i in range(len(self.choice_structures)):
            priority = self.get_priority_count(i)

            if priority < len(self.priorities):
                insert_position: int = self.priorities[priority].insert(i, None)

                self.choice_structures[i].priority = priority
                self.choice_structures[i].position = insert_position

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
    
    def add_gray_edge(self, gray_edge: PGMPath):
        """
        Adds a gray edge object and updates self

        Parameters
        ----------
        gray_edge
            Temporary variable used in calling function
        """
        if gray_edge.head >= 0 and gray_edge.tail >= 0:
            self.gray_edge[self.gray_edge_index] = copy(gray_edge)
            self.gray_edge_index += 1

        self.update_all(gray_edge.head - 1, gray_edge)
    
    def update_all(self, choice_structure_index: int, path_l: PGMPath):
        """
        Updates internal data

        Parameters
        ----------
        choice_structure_index
            Index of the current ChoiceStructure
        path_l
            Used for combining paths
        """
        new_fragment: List[PGMFragment] = self.get_created_fragment(path_l.head, path_l.tail)
        self.get_new_fragment_list(new_fragment)

        head: int
        tail: int

        if path_l.head > self.gene_number * 2:
            head = path_l.head - (self.gene_number * 2)
        else:
            head = path_l.head

        if path_l.tail > self.gene_number * 2:
            tail = path_l.tail - (self.gene_number * 2)
        else:
            tail = path_l.tail

        self.priorities[self.choice_structures[head - 1].priority].remove(self.choice_structures[head - 1].position)
        self.priorities[self.choice_structures[tail - 1].priority].remove(self.choice_structures[tail - 1].position)

        new_choice_structure: List[ChoiceStructure] = self.get_new_choice_structure(choice_structure_index,
                                                                                    path_l.head,
                                                                                    path_l.tail)

        for choice_structure in new_choice_structure:
            index_from: int = choice_structure.index_from
            self.choice_structures[index_from - 1] = copy(choice_structure)

        self.choice_structures[path_l.head - 1] = None
        self.choice_structures[path_l.tail - 1] = None

        for choice_structure in new_choice_structure:
            current_choice_structure_index: int = choice_structure.index_from - 1
            priority: int = self.get_priority_count(current_choice_structure_index)

            if priority != choice_structure.priority:
                old_priority: int = choice_structure.priority
                old_position: int = choice_structure.position

                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int = -1

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(choice_structure.index_from - 1, None)

                self.choice_structures[choice_structure.index_from - 1].priority = priority
                self.choice_structures[choice_structure.index_from - 1].position = new_position

        if new_fragment[2].end1 > 0:
            self.update_priority(new_fragment[2].end1 - 1)

        if new_fragment[2].end2 > 0:
            self.update_priority(new_fragment[2].end2 - 1)

        for choice_structure in new_choice_structure:
            tail1: int
            tail2: int
            tail3: int

            if choice_structure.genome_1_path.tail > self.gene_number * 2:
                tail1 = choice_structure.genome_1_path.tail - (self.gene_number * 2)
            else:
                tail1 = choice_structure.genome_1_path.tail

            if choice_structure.genome_2_path.tail > self.gene_number * 2:
                tail2 = choice_structure.genome_2_path.tail - (self.gene_number * 2)
            else:
                tail2 = choice_structure.genome_2_path.tail

            if choice_structure.genome_3_path.tail > self.gene_number * 2:
                tail3 = choice_structure.genome_3_path.tail - (self.gene_number * 2)
            else:
                tail3 = choice_structure.genome_3_path.tail

            if tail1 > 0:
                self.update_priority(tail1 - 1)

            if tail2 > 0:
                self.update_priority(tail2 - 1)

            if tail3 > 0:
                self.update_priority(tail3 - 1)

    def get_ancestors(self):
        median_chromosome: int = 0

        for fragment in self.fragments:
            if fragment is not None:
                self.fragments[fragment.end2] = None

        for fragment in self.fragments:
            if fragment is not None:
                median_chromosome += 1

        self.ancestor_AA = list()

        for fragment in self.fragments:
            if fragment is not None:
                start_index: int = fragment.end1
                end_index: int = fragment.end2

                self.ancestor_AA.append(self.get_chromosome_using_start_gene(start_index, end_index, 2))
                self.ancestor_AA.append(self.get_chromosome_using_start_gene(start_index + self.gene_number * 2,
                                                                             end_index,
                                                                             2))

        self.ancestor_A = list()

        for fragment in self.fragments:
            if fragment is not None:
                start_index: int = fragment.end1
                end_index: int = fragment.end2

                self.ancestor_A.append(self.get_chromosome_using_start_gene(start_index, end_index, 1))

    def get_priority_count(self, cs_index: int) -> int:
        ancestor_priority: int
        result: int = 200
        to_replace: int = -1
        from_1: int = self.choice_structures[cs_index].genome_1_path.head
        from_2: int = self.choice_structures[cs_index].genome_2_path.head
        tail_1: int = self.choice_structures[cs_index].genome_1_path.tail
        tail_2: int = self.choice_structures[cs_index].genome_2_path.tail

        rng: Random = Random()
        rng.seed()

        if self.replace == 0:
            to_replace = rng.randint(0, 1)
        elif self.replace == 1:
            to_replace = 0
        elif self.replace == 2:
            to_replace = 1

        ancestor_priority = self.calculate_case(cs_index, from_1, tail_1)

        if ancestor_priority == 0 or ancestor_priority == 1:
            self.choice_structures[cs_index].gray_edge = PGMPath(from_1, tail_1, None, None)

            return ancestor_priority

        if ancestor_priority < result:
            result = ancestor_priority
            self.choice_structures[cs_index].gray_edge = PGMPath(from_1, tail_1, None, None)

        elif ancestor_priority == result and to_replace == 1:
            result = ancestor_priority
            self.choice_structures[cs_index].gray_edge = PGMPath(from_1, tail_1, None, None)
        
        to_replace = 1 - to_replace

        ancestor_priority = self.calculate_case(cs_index, from_2, tail_2)

        if ancestor_priority == 1:
            self.choice_structures[cs_index].gray_edge = PGMPath(from_2, tail_2, None, None)

            return ancestor_priority

        if ancestor_priority < result:
            result = ancestor_priority
            self.choice_structures[cs_index].gray_edge = PGMPath(from_2, tail_2, None, None)

        elif ancestor_priority == result and to_replace == 1:
            result = ancestor_priority
            self.choice_structures[cs_index].gray_edge = PGMPath(from_2, tail_2, None, None)

        return result

    def get_created_fragment(self, ancestor_1: int, ancestor_2: int) -> List[PGMFragment]:
        node_1: int = ancestor_1
        node_2: int = ancestor_2
        result: List[PGMFragment] = list()

        if node_1 > 0 and node_2 > 0:
            if node_1 > self.gene_number * 2:
                node_1 -= self.gene_number * 2

            if node_2 > self.gene_number * 2:
                node_2 -= self.gene_number * 2

            result.append(PGMFragment.from_fragment(self.fragments[node_1]))
            result.append(PGMFragment.from_fragment(self.fragments[node_2]))
            result.append(combine(PGMPath(node_1, node_2, None, None), self.fragments[node_1], self.fragments[node_2]))

        return result

    def get_new_fragment_list(self, created_frags: List[PGMFragment]):
        if len(created_frags) != 0:
            self.fragments[created_frags[0].end1] = None
            self.fragments[created_frags[0].end2] = None
            self.fragments[created_frags[1].end1] = None
            self.fragments[created_frags[1].end2] = None
            self.fragments[created_frags[2].end1] = created_frags[2]
            self.fragments[created_frags[2].end2] = PGMFragment(created_frags[2].end2, created_frags[2].end1)
    
    def get_new_choice_structure(self, choice_structure_index: int, index_from: int, tail: int) -> \
            List[Optional[ChoiceStructure]]:
        """
        Gets new ChoiceStructure in some time that isn't 2 steps

        Parameters
        ----------
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
        f: int = index_from
        t: int = tail

        if index_from > self.gene_number * 2:
            f = index_from - (self.gene_number * 2)

        if tail > self.gene_number * 2:
            t = tail - (self.gene_number * 2)

        if f != choice_structure_index + 1:
            print("Wrong choice_structure_index: " + str(choice_structure_index))

        path_1_1: PGMPath = \
            copy(self.choice_structures[choice_structure_index].genome_1_path)
        path_1_2: PGMPath = copy(self.choice_structures[tail - 1].genome_1_path)
        path_2_1: PGMPath = \
            copy(self.choice_structures[choice_structure_index].genome_2_path)
        path_2_2: PGMPath = copy(self.choice_structures[tail - 1].genome_2_path)
        path_3_1: PGMPath = \
            copy(self.choice_structures[choice_structure_index].genome_3_path)
        path_3_2: PGMPath = copy(self.choice_structures[tail - 1].genome_3_path)

        path_l1: PGMPath = PGMPath(index_from, tail, None, None)

        f2: int

        if index_from > self.gene_number * 2:
            f2 = f
        else:
            f2 = index_from + (self.gene_number * 2)

        t2: int

        if tail > self.gene_number * 2:
            t2 = t
        else:
            t2 = tail + (self.gene_number * 2)

        path_l2: PGMPath = PGMPath(f2, t2, None, None)

        new_path1: Optional[PGMPath] = None
        new_path2: Optional[PGMPath] = None

        if index_from <= self.gene_number * 2 and tail <= self.gene_number * 2:
            new_path1 = path_1_1.connect(path_1_1, path_1_2, path_l1, None)
            new_path2 = path_2_1.connect(path_2_1, path_2_2, path_l2, None)
        elif index_from > self.gene_number * 2 and tail > self.gene_number * 2:
            new_path1 = path_1_1.connect(path_1_1, path_1_2, path_l2, None)
            new_path2 = path_2_1.connect(path_2_1, path_2_2, path_l1, None)
        elif index_from <= self.gene_number * 2 < tail:
            new_path1 = path_1_1.connect(path_1_1, path_2_2, path_l1, None)
            new_path2 = path_2_1.connect(path_2_1, path_1_2, path_l2, None)
        elif index_from > self.gene_number * 2 >= tail:
            new_path1 = path_1_1.connect(path_2_1, path_1_2, path_l1, None)
            new_path2 = path_2_1.connect(path_1_1, path_2_2, path_l2, None)

        path_l3: PGMPath = PGMPath(f, t, None, None)
        new_path3: PGMPath = path_3_1.connect(path_3_1, path_3_2, path_l3, None)

        temp: List[Optional[ChoiceStructure]] = [None for _ in range(6)]

        if not self.is_a_cycle(new_path1, index_from, tail):
            temp = self.get_new_choice_structure_based_on_path(new_path1, temp, None, 2)
        if not self.is_a_cycle(new_path2, index_from, tail):
            temp = self.get_new_choice_structure_based_on_path(new_path2, temp, None, 2)
        if not self.is_a_cycle(new_path3, index_from, tail):
            temp = self.get_new_choice_structure_based_on_path(new_path3, temp, None, 1)

        new_choice_structure_number: int = 0

        for choice_structure in temp:
            if choice_structure is not None:
                new_choice_structure_number += 1

        result: List[Optional[ChoiceStructure]] = [None for _ in range(new_choice_structure_number)]

        if len(result) >= 0:
            result = deepcopy(temp)[:len(result)]

        return result

    def is_a_cycle(self, path1: PGMPath, f2: int, t2: int) -> bool:
        f1: int = path1.head
        t1: int = path1.tail

        if f1 > self.gene_number * 2:
            f1 -= self.gene_number * 2

        if t1 > self.gene_number * 2:
            t1 -= self.gene_number * 2

        f: int = f2
        t: int = t2

        if f > self.gene_number * 2:
            f -= self.gene_number * 2

        if t > self.gene_number * 2:
            t -= self.gene_number * 2

        return (f1 == f and t1 == t) or (f1 == t and t1 == f)

    def get_new_choice_structure_based_on_path(self, new_path1: PGMPath,
                                               temp: List[Optional[ChoiceStructure]],
                                               new_choice_structures: Optional[List[ChoiceStructure]],
                                               ploidy: int) \
            -> List[ChoiceStructure]:
        """
        Gets a new ChoiceStructure based on given PGMPath

        Parameters
        ----------
        new_path1
            Path to base ChoiceStructure on
        temp
            Modified to add the new ChoiceStructure to
        new_choice_structures
            New ChoiceStructures
        ploidy
            Monoploid or diploid

        Returns
        -------
        List[ChoiceStructure]
            temp, modified with the new ChoiceStructure
        """
        from_1: int = new_path1.head
        from_small: int

        if from_1 > self.gene_number * 2:
            from_small = from_1 - (self.gene_number * 2)
        else:
            from_small = from_1

        tail_1: int = new_path1.tail
        tail_small: int

        if tail_1 > self.gene_number * 2:
            tail_small = tail_1 - (self.gene_number * 2)
        else:
            tail_small = tail_1

        temp_index: int = len(temp) - temp.count(None)

        if from_1 > 0:
            small_index: int = check_temp_list(temp, from_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:
                    for choice_structure in new_choice_structures:
                        if choice_structure.index_from == from_small:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(choice_structure)
                            temp[temp_index].set_new_path(new_path1, ploidy, self.gene_number)
                            find_now = True
                            break

                if not find_now and self.choice_structures[from_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure()
                    temp[temp_index].from_cs(self.choice_structures[from_small - 1])
                    temp[temp_index].set_new_path(new_path1, ploidy, self.gene_number)
                    temp_index += 1
            else:
                temp[small_index].set_new_path(new_path1, ploidy, self.gene_number)

        if tail_1 > 0:
            np1: PGMPath = PGMPath(tail_1, from_1, None, None)
            small_index: int = check_temp_list(temp, tail_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:
                    for choice_structure in new_choice_structures:
                        if choice_structure.index_from == tail_small:
                            temp[temp_index] = ChoiceStructure()
                            temp[temp_index].from_cs(choice_structure)
                            temp[temp_index].set_new_path(new_path1, ploidy, self.gene_number)
                            find_now = True
                            break

                if not find_now and self.choice_structures[tail_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure()
                    temp[temp_index].from_cs(self.choice_structures[tail_small - 1])
                    temp[temp_index].set_new_path(np1, ploidy, self.gene_number)
                    temp_index += 1
            else:
                temp[small_index].set_new_path(np1, ploidy, self.gene_number)

        return temp

    def update_priority(self, cs_index: int):
        if self.choice_structures[cs_index] is not None:
            old_priority: int = self.choice_structures[cs_index].priority
            old_position: int = self.choice_structures[cs_index].position
            priority: int = self.get_priority_count(cs_index)

            if priority != old_priority:
                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(cs_index, None)
                else:
                    new_position = -1

                self.choice_structures[cs_index].priority = priority
                self.choice_structures[cs_index].position = new_position

    def get_chromosome_using_start_gene(self, start_index: int, end_index: int, ploidy: int) -> str:
        ancestor_chromosome: str = str()

        if ploidy == 1:
            start: str = self.node_str[start_index - 1]

            if start.endswith("h"):
                ancestor_chromosome = "-" + start[:-1]
            else:
                ancestor_chromosome = start[:-1]

            start_index = get_gene_next_node(start_index)

            while True:
                next_gene_index: int = self.find_gray_edge_node(start_index, self.gray_edge, 1)

                if next_gene_index == -1000:
                    break

                next_gene: str = self.node_str[next_gene_index - 1]
                current_start_index: int

                if start_index > self.gene_number * 2:
                    current_start_index = start_index - (self.gene_number * 2)
                else:
                    current_start_index = start_index

                if current_start_index == end_index:
                    break
                else:
                    if next_gene.endswith("h"):
                        ancestor_chromosome += " -" + next_gene[:-1]
                    else:
                        ancestor_chromosome += " " + next_gene[:-1]

                    start_index = get_gene_next_node(next_gene_index)
        elif ploidy == 2:
            if start_index > self.gene_number * 2:
                start: str = self.node_str[start_index - self.gene_number * 2 - 1]

                if start.endswith("h"):
                    ancestor_chromosome = "-" + start[:-1] + "b"
                else:
                    ancestor_chromosome = start[:-1] + "b"
            else:
                start: str = self.node_str[start_index - 1]

                if start.endswith("h"):
                    ancestor_chromosome = "-" + start[:-1] + "a"
                else:
                    ancestor_chromosome = start[:-1] + "a"

            start_index = get_gene_next_node(start_index)

            while True:
                next_gene_index: int = self.find_gray_edge_node(start_index, self.gray_edge, 2)

                if next_gene_index == -1000:
                    break

                next_gene: str

                if next_gene_index > self.gene_number * 2:
                    next_gene = self.node_str[next_gene_index - self.gene_number * 2 - 1]
                else:
                    next_gene = self.node_str[next_gene_index - 1]

                current_start_index: int

                if start_index > self.gene_number * 2:
                    current_start_index = start_index - self.gene_number * 2
                else:
                    current_start_index = start_index

                if current_start_index == end_index:
                    break
                else:
                    if next_gene.endswith("h"):
                        if next_gene_index > self.gene_number * 2:
                            ancestor_chromosome += " -" + next_gene[:-1] + "b"
                        else:
                            ancestor_chromosome += " -" + next_gene[:-1] + "a"
                    else:
                        if next_gene_index > self.gene_number * 2:
                            ancestor_chromosome += " " + next_gene[:-1] + "b"
                        else:
                            ancestor_chromosome += " " + next_gene[:-1] + "a"

                    start_index = get_gene_next_node(next_gene_index)

        return ancestor_chromosome
    
    def calculate_case(self, choice_structure_index: int, index_from: int, tail: int) -> int:
        """
        Figures out which priority case to use (combination of bcla, bw, cn)

        Parameters
        ----------
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

        cycle_now: int = self.count_cycle_for_edge(choice_structure_index, tail)

        if cycle_now == 3:
            return 0

        created_fragment: List[PGMFragment] = self.get_created_fragment(index_from, tail)
        self.get_new_fragment_list(created_fragment)
        created_choice_structure: List[ChoiceStructure] = self.get_new_choice_structure(choice_structure_index,
                                                                                        index_from,
                                                                                        tail)

        look_ahead_cycles: List[int] = self.count_all_look_ahead_cycles(created_choice_structure)
        max_cycle_look_ahead: int = look_ahead_cycles[0]
        how_many_more_cycles: int = look_ahead_cycles[1]
        max_cycle_look_ahead2: int = look_ahead_cycles[2]

        if len(created_fragment) != 0:  # recoverOriginalFragment()
            self.fragments[created_fragment[2].end1] = None
            self.fragments[created_fragment[2].end2] = None
            self.fragments[created_fragment[0].end1] = created_fragment[0]
            self.fragments[created_fragment[0].end2] = PGMFragment(created_fragment[0].end2, created_fragment[0].end1)
            self.fragments[created_fragment[1].end1] = created_fragment[1]
            self.fragments[created_fragment[1].end2] = PGMFragment(created_fragment[1].end2, created_fragment[1].end1)

        index: int = 0

        for cn in range(2, 0, -1):
            for bcla in range(3, 0, -1):
                for bw in range(4, -5, -1):
                    for bcla2 in range(3, 0, -1):
                        if cycle_now == cn and max_cycle_look_ahead == bcla and \
                                how_many_more_cycles == bw and max_cycle_look_ahead2 == bcla2:
                            return index

                        index += 1

        return 200
    
    def find_gray_edge_node(self, node: int, gray_edge: List[Optional[PGMPath]], ploidy: int) -> int:
        if ploidy == 1:
            for gray_node in gray_edge:
                if gray_node is not None:
                    head: int = gray_node.head
                    tail: int = gray_node.tail
                    
                    if node == head:
                        if tail > self.gene_number * 2:
                            return tail - self.gene_number * 2
                        else:
                            return tail
                    
                    if node == tail:
                        if head > self.gene_number * 2:
                            return head - self.gene_number * 2
                        else:
                            return head
        elif ploidy == 2:
            for gray_node in gray_edge:
                if gray_node is not None:
                    head: int = gray_node.head
                    tail: int = gray_node.tail
                    
                    if node == head:
                        return tail
                    
                    if node == tail:
                        return head
                    
        return -1000
    
    def count_cycle_for_edge(self, choice_structure_index: int, tail: int) -> int:
        """
        Counts the cycles for the current edge

        Parameters
        ----------
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
        head: int = self.choice_structures[choice_structure_index].index_from

        if 0 < head == self.fragments[head].end1 and \
                self.fragments[head] is not None and \
                self.fragments[head].end2 == tail:
            is_circular_chromosome = True
        elif 0 < tail == self.fragments[tail].end1 and \
                self.fragments[tail] is not None and \
                self.fragments[tail].end2 == head:
            is_circular_chromosome = True

        if is_circular_chromosome:
            return -1

        if head == tail + self.gene_number * 2 or tail == head + self.gene_number * 2:
            return -1

        result: int = 0

        if tail > 0:
            ancestor_choice_structure: ChoiceStructure = self.choice_structures[
                choice_structure_index]

            f1: int = ancestor_choice_structure.genome_1_path.head
            t1: int = ancestor_choice_structure.genome_1_path.tail
            f2: int = ancestor_choice_structure.genome_2_path.head
            t2: int = ancestor_choice_structure.genome_2_path.tail
            f3: int = ancestor_choice_structure.genome_3_path.head
            t3: int = ancestor_choice_structure.genome_3_path.tail

            if f1 > self.gene_number * 2:
                f1 -= self.gene_number * 2

            if t1 > self.gene_number * 2:
                t1 -= self.gene_number * 2

            if f2 > self.gene_number * 2:
                f2 -= self.gene_number * 2

            if t2 > self.gene_number * 2:
                t2 -= self.gene_number * 2

            if f1 == head and t1 == tail:
                result += 1

            if f2 == head and t2 == tail:
                result += 1

            if f3 == head and t3 == tail:
                result += 1

        return result

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
                old_cycle: int = count_cycle_look_ahead(self.choice_structures[index_from - 1])
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
                    elif tail2 == tail3:
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
                    elif genome_tail2 == genome_head2:
                        current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                            index_from, tail2, current_genome, ancestor_choice_structure, created_choice_structures)

                        for new_choice_structure in current_new_choice_structure:
                            new_count = count_cycle_look_ahead(new_choice_structure)

                            if new_count == 3:
                                return new_count
                            elif current_max < new_count:
                                current_max = new_count
                    elif genome_tail3 == genome_head3:
                        current_new_choice_structure: [ChoiceStructure] = self.get_new_choice_structure_2_step(
                            index_from, tail3, current_genome, ancestor_choice_structure, created_choice_structures)

                        for new_choice_structure in current_new_choice_structure:
                            new_count = count_cycle_look_ahead(new_choice_structure)

                            if new_count == 3:
                                return new_count
                            elif current_max < new_count:
                                current_max = new_count

        return current_max

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
            new_choice_structure.from_cs(self.choice_structures[tail - 1])

        path_1_2: PGMPath = new_choice_structure.genome_1_path
        path_2_2: PGMPath = new_choice_structure.genome_2_path
        path_3_2: PGMPath = new_choice_structure.genome_3_path

        from_2: int

        if index_from > self.gene_number * 2:
            from_2 = index_from - self.gene_number * 2
        else:
            from_2 = index_from + self.gene_number * 2

        to_2: int

        if tail > self.gene_number * 2:
            to_2 = tail - self.gene_number * 2
        else:
            to_2 = tail + self.gene_number * 2

        from_3: int = index_from

        if index_from > self.gene_number * 2:
            from_3 = from_2

        to_3: int = tail

        if tail > self.gene_number * 2:
            to_3 = to_2

        path_l1: PGMPath = PGMPath(index_from, tail, None, None)
        path_l2: PGMPath = PGMPath(from_2, to_2, None, None)
        path_l3: PGMPath = PGMPath(from_3, to_3, None, None)

        new_path1: Optional[PGMPath] = None
        new_path2: Optional[PGMPath] = None

        if index_from <= self.gene_number * 2 and tail <= self.gene_number * 2:
            new_path1 = path_1_1.connect(path_1_1, path_1_2, path_l1, None)
            new_path2 = path_2_1.connect(path_2_1, path_2_2, path_l2, None)

        if index_from > self.gene_number * 2 and tail > self.gene_number * 2:
            new_path1 = path_1_1.connect(path_1_1, path_1_2, path_l2, None)
            new_path2 = path_2_1.connect(path_2_1, path_2_2, path_l1, None)

        if index_from <= self.gene_number * 2 < tail:
            new_path1 = path_1_1.connect(path_1_1, path_2_2, path_l1, None)
            new_path2 = path_2_1.connect(path_2_1, path_1_2, path_l2, None)

        if index_from > self.gene_number * 2 >= tail:
            new_path1 = path_1_1.connect(path_2_1, path_1_2, path_l1, None)
            new_path2 = path_2_1.connect(path_1_2, path_2_2, path_l2, None)

        new_path3: Optional[PGMPath] = path_3_1.connect(path_3_1, path_3_2, path_l3, None)

        temp: List[Optional[ChoiceStructure]] = [None for _ in range(6)]

        if not is_a_cycle(new_path1, PGMPath(index_from, tail, None, None)):
            temp = self.get_new_choice_structure_based_on_path(new_path1, temp, new_choice_structures, 2)
        if not is_a_cycle(new_path2, PGMPath(index_from, tail, None, None)):
            temp = self.get_new_choice_structure_based_on_path(new_path2, temp, new_choice_structures, 2)
        if not is_a_cycle(new_path3, PGMPath(index_from, tail, None, None)):
            temp = self.get_new_choice_structure_based_on_path(new_path3, temp, new_choice_structures, 1)

        new_choice_structure_number: int = 0

        for choice_structure in temp:
            if choice_structure is not None:
                new_choice_structure_number += 1

        result: List[Optional[ChoiceStructure]] = [None for _ in range(new_choice_structure_number)]

        if len(result) >= 0:
            result = deepcopy(temp)[:len(result)]

        return result
