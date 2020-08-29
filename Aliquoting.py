from random import Random
from typing import List, Optional, Dict, Any

import ChoiceStructure
import PGMPath
from Chromosome import Chromosome
from Genome import Genome, split_at_whitespace
from PGMFragment import PGMFragment, combine
from Priority import Priority


def count_gene_number(genome: Genome) -> int:
    """
    Count number of genes in a genome

    Parameters
    ----------
    genome
        Genome to count genes of

    Returns
    -------
    int
        Number of genes in the genome
    """
    result: int = 0

    for chromosome in genome.chromosomes:
        result += len(chromosome.genes)

    return result


def check_temp_list(temp_list: List[Dict[str, Any]], f: int) -> int:
    """
    Finds the index of the first ChoiceStructure in the list with the given head

    Parameters
    ----------
    temp_list
        List to check
    f
        Head value to check for

    Returns
    -------
    int
        Index of the ChoiceStructure with the given head
    """
    for i in range(len(temp_list)):
        if temp_list[i] is not None:
            if temp_list[i]["genome_1_path"]["head"] == f:
                return i

    return -1


def get_gene_next_node(index: int) -> int:
    """
    Find paired gene_node

    Parameters
    ----------
    index
        Index to search

    Returns
    -------
    int
        Index of the paired gene_node
    """
    if (index // 2) * 2 == index:
        return index - 1

    return index + 1


ploidys: Dict[int, str] = dict([(c - 97, chr(c)) for c in range(ord("a"), ord("z") + 1)])


def find_gray_edge_node(node: int, gray_edge: List[Optional[Dict[str, int]]]) -> int:
    """
    Find gray edge node for the given node index

    Parameters
    ----------
    node
        Node index
    gray_edge
        GrayEdge object

    Returns
    -------
    int
        Index of the gray edge node
    """
    for gray_node in gray_edge:
        if gray_node is not None:
            head: int = gray_node["head"]
            tail: int = gray_node["tail"]

            if node == head:
                return tail

            if node == tail:
                return head

    return -1000


class Aliquoting:
    """
    Attributes
    ----------
    node_int
        Nodes in integer form
    node_str
        Nodes in string form
    ancestors
        Ancestor genomes
    replace
        Which genome to replace
    gene_number
        Number of genes in outgroup genome
    priorities
        Priority list
    polyd
        Polyploid genome
    ploidy
        Polyploid ploidy / 2 (1 = diploid, 2 = tetraploid, 3 = hexaploid, etc)
    fragments
        List of PGMFragments for PathGroups
    gray_edge
        Gray edge
    gray_edge_index
        Index of gray edge
    choice_structures
        List of ChoiceStructures
    """
    def __init__(self, polyd: Genome, reference: Genome, replace: int, ploidy: int):
        """
        Constructor

        Parameters
        ----------
        polyd
            Polyploid genome to use
        reference
            Reference genome
        replace
            Which genome to replace (< ploidy, or 0 for random)
        ploidy
            Number of gene copy sets in the polyploid genome (n)
        """
        self.node_int: List[int] = list()
        self.node_str: List[str] = list()
        self.ancestors: List[Genome] = list()

        self.replace: int = replace
        self.gene_number: int = count_gene_number(reference)
        self.ploidy: int = ploidy

        self.priorities: List[Priority] = self.initialize_priorities()
        self.set_nodes(reference)
        self.polyd: List[Dict[str, int]] = self.get_pgm_path(polyd)

        self.fragments: List[Optional[PGMFragment]] = [
            None for _ in range(self.gene_number * (self.ploidy - 1) * 2 + 1)]

        for i in range(int(len(self.fragments) / 2)):
            self.fragments[2 * i + 1] = PGMFragment(2 * i + 1, 2 * i + 2)
            self.fragments[2 * i + 2] = PGMFragment(2 * i + 2, 2 * i + 1)

        self.gray_edge: List[Optional[Dict[str, int]]] = [None for _ in range(len(self.polyd))]
        self.gray_edge_index: int = 0

        self.choice_structures: List[Optional[Dict[str, Any]]] = [
            None for _ in range(self.gene_number * (self.ploidy - 1) * 2)]

        cs_index: int = 0

        for i in range(1, self.gene_number * (self.ploidy - 1) * 2 + 1):
            self.choice_structures[cs_index] = ChoiceStructure.create_cs(ploidy=self.ploidy)
            self.choice_structures[cs_index]["index_from"] = i

            for j in range(self.ploidy):
                self.choice_structures[cs_index]["genome_paths"][j] = self.polyd[i + (self.gene_number * j * 2)]

            # self.choice_structures[cs_index]["genome_1_path"] = self.polyd[i]
            # self.choice_structures[cs_index]["genome_2_path"] = self.polyd[i + self.gene_number * 2]
            # self.choice_structures[cs_index]["genome_3_path"] = self.outgroup[i]
            self.choice_structures[cs_index]["priority"] = 200
            self.choice_structures[cs_index]["position"] = -1
            self.choice_structures[cs_index]["gray_edge"] = None
            cs_index += 1

    def initialize_priorities(self) -> List[Priority]:
        """
        Initialize priority list

        Returns
        -------
        List[Priority]
            Priority list
        """
        priority_size: int = self.gene_number * (self.ploidy - 1) * 2 + 1000
        priorities: List[Priority] = list()

        priorities.append(Priority(priority_size))

        for cn in range(2, 0, -1):
            for bcla in range(3, 0, -1):
                for bw in range(4, -5, -1):
                    for bcla2 in range(3, 0, -1):
                        priorities.append(Priority(priority_size))

        return priorities

    def set_nodes(self, genome: Genome):
        """
        Sets node_int and node_str from the given Genome

        Parameters
        ----------
        genome
            Genome to generate from
        """
        self.node_int = [int() for _ in range(self.gene_number * (self.ploidy - 1) * 2)]
        self.node_str = [str() for _ in range(self.gene_number * (self.ploidy - 1) * 2)]

        index: int = 0

        for chromosome in genome.chromosomes:
            for gene in chromosome.genes:
                node_1: str
                node_2: str

                if gene.name[0] == "-":
                    node_1 = str().join([gene.name[1:], "h"])
                    node_2 = str().join([gene.name[1:], "t"])

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_2
                    index += 1

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_1
                    index += 1
                else:
                    node_1 = str().join([gene.name, "t"])
                    node_2 = str().join([gene.name, "h"])

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_1
                    index += 1

                    self.node_int[index] = index + 1
                    self.node_str[index] = node_2
                    index += 1

    def get_pgm_path(self, genome: Genome) -> List[Dict[str, int]]:
        """
        Gets list of PGMPaths for a genome

        Parameters
        ----------
        genome
            Genome to get the PGMPath for

        Returns
        -------
        List[Dict[str, int]]
            List of PGMPaths for the genome
        """
        path1: List[Optional[Dict[str, int]]] = [None for _ in range((4 * self.gene_number) + 1)]

        null_node: int = -1

        for chromosome in genome.chromosomes:
            pre_node: int = 0

            for j in range(len(chromosome.genes)):
                first_character: str = chromosome.genes[j].name[0]
                node1: str
                node2: str

                if first_character == '-':
                    node1 = str().join([chromosome.genes[j].name[1:], "h"])
                    node2 = str().join([chromosome.genes[j].name[1:], "t"])
                else:
                    node1 = str().join([chromosome.genes[j].name, "t"])
                    node2 = str().join([chromosome.genes[j].name, "h"])

                node1_int: int = self.find_node_int(node1)
                node2_int: int = self.find_node_int(node2)

                if node1_int == 0 or node2_int == 0:
                    print("Gene ", str(chromosome.genes[j]), " does not exist in the other genome.\n")

                if j == 0:
                    path1[node1_int] = PGMPath.create_pgm_path(node1_int, null_node)
                    pre_node = node2_int
                    null_node -= 1
                elif j != 0 and j != len(chromosome.genes) - 1:
                    path1[node1_int] = PGMPath.create_pgm_path(node1_int, pre_node)
                    path1[pre_node] = PGMPath.create_pgm_path(pre_node, node1_int)
                    pre_node = node2_int

                if j == len(chromosome.genes) - 1:
                    if len(chromosome.genes) != 1:
                        path1[pre_node] = PGMPath.create_pgm_path(pre_node, node1_int)
                        path1[node1_int] = PGMPath.create_pgm_path(node1_int, pre_node)

                    path1[node2_int] = PGMPath.create_pgm_path(node2_int, null_node)
                    null_node -= 1

        return path1

    def find_node_int(self, ancestor_string: str) -> int:
        """
        Finds the node_int from node_str for the given ancestor_string

        Parameters
        ----------
        ancestor_string
            Ancestor string

        Returns
        -------
        int
            Index in node_int of the ancestor_string
        """
        copy_: str = ancestor_string[-2]
        name: str = ancestor_string[:-2] + ancestor_string[-1]

        for i in range(len(self.node_str)):
            if self.node_str[i] == name:
                if copy_ == "a":
                    return i + 1
                elif copy_ in ploidys.values():
                    return (i + 1) + self.gene_number * (ord(copy_) - 97) * 2

        return 0

    def get_result(self):
        """
        Generates GGH result from the list of ChoiceStructures
        """
        self.group_pathgroup_into_priorities()

        best_cs: List[int] = self.find_the_best_choice_structure()

        while best_cs[0] != -1:
            self.add_gray_edge(self.choice_structures[self.priorities[best_cs[0]].cs_indexes[best_cs[1]]]["gray_edge"])
            best_cs = self.find_the_best_choice_structure()

        self.get_ancestors()

    def group_pathgroup_into_priorities(self):
        """
        Groups path group into the instance priorities list
        """
        for i in range(len(self.choice_structures)):
            priority = self.get_priority_count(i)

            if priority < len(self.priorities):
                insert_position: int = self.priorities[priority].insert(i)

                self.choice_structures[i]["priority"] = priority
                self.choice_structures[i]["position"] = insert_position

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

    def add_gray_edge(self, gray_edge: Dict[str, int]):
        """
        Adds a gray edge object and updates self

        Parameters
        ----------
        gray_edge
            Temporary variable used in calling function
        """
        if gray_edge["head"] >= 0 and gray_edge["tail"] >= 0:
            from_1: int = 0
            tail_1: int = 0

            for i in range(self.ploidy - 1, 0, -1):
                if from_1 == 0 and gray_edge["head"] > self.gene_number * i * 2:
                    from_1 = gray_edge["head"] - self.gene_number * i * 2

                if tail_1 == 0 and gray_edge["tail"] > self.gene_number * i * 2:
                    tail_1 = gray_edge["tail"] - self.gene_number * i * 2

            if from_1 == 0:
                from_1 = gray_edge["head"] + self.gene_number * (self.ploidy - 1) * 2

            if tail_1 == 0:
                tail_1 = gray_edge["tail"] + self.gene_number * (self.ploidy - 1) * 2

            # if gray_edge["head"] > self.gene_number * 2:
            #     from_1 = gray_edge["head"] - self.gene_number * 2
            # else:
            #     from_1 = gray_edge["head"] + self.gene_number * 2
            #
            # if gray_edge["tail"] > self.gene_number * 2:
            #     tail_1 = gray_edge["tail"] - self.gene_number * 2
            # else:
            #     tail_1 = gray_edge["tail"] + self.gene_number * 2

            self.gray_edge[self.gray_edge_index] = gray_edge
            self.gray_edge_index += 1
            self.gray_edge[self.gray_edge_index] = PGMPath.create_pgm_path(from_1, tail_1)
            self.gray_edge_index += 1

        cs_index: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if gray_edge["head"] > self.gene_number * i * 2:
                cs_index = (gray_edge["head"] - 1) - self.gene_number * i * 2

        if cs_index == -10000:
            cs_index = gray_edge["head"] - 1

        # if gray_edge["head"] > self.gene_number * 2:
        #     cs_index = (gray_edge["head"] - 1) - self.gene_number * 2
        # else:
        #     cs_index = gray_edge["head"] - 1

        self.update_all(cs_index, gray_edge)

    def update_all(self, choice_structure_index: int, path_l: Dict[str, int]):
        """
        Updates internal data

        Parameters
        ----------
        choice_structure_index
            Index of the current ChoiceStructure
        path_l
            Used for combining paths
        """
        new_fragment: List[PGMFragment] = self.get_created_fragment(path_l["head"], path_l["tail"])
        self.get_new_fragment_list(new_fragment)

        start_index: int = -10000
        end_index: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if start_index == -10000 and path_l["head"] > self.gene_number * i * 2:
                start_index = path_l["head"] - self.gene_number * i * 2

            if end_index == - 10000 and path_l["tail"] > self.gene_number * i * 2:
                end_index = path_l["tail"] - self.gene_number * i * 2

        if start_index == -10000:
            start_index = path_l["head"]

        if end_index == -10000:
            end_index = path_l["tail"]

        # if path_l["head"] > self.gene_number * 2:
        #     start_index = path_l["head"] - (self.gene_number * 2)
        # else:
        #     start_index = path_l["head"]
        #
        # if path_l["tail"] > self.gene_number * 2:
        #     end_index = path_l["tail"] - (self.gene_number * 2)
        # else:
        #     end_index = path_l["tail"]

        self.priorities[self.choice_structures[start_index - 1]["priority"]].\
            remove(self.choice_structures[start_index - 1]["position"])
        self.priorities[self.choice_structures[end_index - 1]["priority"]].\
            remove(self.choice_structures[end_index - 1]["position"])

        new_choice_structure: List[Dict[str, Any]] = self.get_new_choice_structure(choice_structure_index,
                                                                                   path_l["head"],
                                                                                   path_l["tail"])

        for choice_structure in new_choice_structure:
            index_from: int = choice_structure["index_from"]
            self.choice_structures[index_from - 1] = choice_structure

        self.choice_structures[start_index - 1] = None
        self.choice_structures[end_index - 1] = None

        for choice_structure in new_choice_structure:
            current_choice_structure_index: int = choice_structure["index_from"] - 1
            priority: int = self.get_priority_count(current_choice_structure_index)

            if priority != choice_structure["priority"]:
                old_priority: int = choice_structure["priority"]
                old_position: int = choice_structure["position"]

                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int = -1

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(choice_structure["index_from"] - 1)

                self.choice_structures[choice_structure["index_from"] - 1]["priority"] = priority
                self.choice_structures[choice_structure["index_from"] - 1]["position"] = new_position

        if new_fragment[2].end1 > 0:
            self.update_priority(new_fragment[2].end1 - 1)

        if new_fragment[2].end2 > 0:
            self.update_priority(new_fragment[2].end2 - 1)

        for choice_structure in new_choice_structure:
            tails: List[int] = list()

            for i in range(self.ploidy):
                tails.append(-10000)

                for j in range(self.ploidy - 1, 0, -1):
                    if choice_structure["genome_paths"][i]["tail"] > self.gene_number * j * 2:
                        tails[i] = choice_structure["genome_paths"][i]["tail"] - self.gene_number * j * 2

                if tails[i] == -10000:
                    tails[i] = choice_structure["genome_paths"][i]["tail"]

            for i in range(self.ploidy):
                if tails[i] > 0:
                    self.update_priority(tails[i] - 1)

            # tail1: int
            # tail2: int
            # tail3: int
            #
            # if choice_structure["genome_1_path"]["tail"] > self.gene_number * 2:
            #     tail1 = choice_structure["genome_1_path"]["tail"] - (self.gene_number * 2)
            # else:
            #     tail1 = choice_structure["genome_1_path"]["tail"]
            #
            # if choice_structure["genome_2_path"]["tail"] > self.gene_number * 2:
            #     tail2 = choice_structure["genome_2_path"]["tail"] - (self.gene_number * 2)
            # else:
            #     tail2 = choice_structure["genome_2_path"]["tail"]
            #
            # if choice_structure["genome_3_path"]["tail"] > self.gene_number * 2:
            #     tail3 = choice_structure["genome_3_path"]["tail"] - (self.gene_number * 2)
            # else:
            #     tail3 = choice_structure["genome_3_path"]["tail"]
            #
            # if tail1 > 0:
            #     self.update_priority(tail1 - 1)
            #
            # if tail2 > 0:
            #     self.update_priority(tail2 - 1)
            #
            # if tail3 > 0:
            #     self.update_priority(tail3 - 1)

    def get_ancestors(self):
        """
        Finds ancestor genomes using PathGroups algorithm
        """
        median_chromosome: int = 0

        for fragment in self.fragments:
            if fragment is not None:
                self.fragments[fragment.end2] = None

        for fragment in self.fragments:
            if fragment is not None:
                median_chromosome += 1

        for i in range(self.ploidy):
            self.ancestors.append(Genome(list()))

            for fragment in self.fragments:
                if fragment is not None:
                    start_index: int = fragment.end1
                    end_index: int = fragment.end2

                    for j in range(i + 1):
                        self.ancestors[i].add_chromosome(
                            self.get_chromosome_using_start_gene(start_index + self.gene_number * j * 2, end_index))

        # self.ancestor_AA = Genome(list())
        #
        # for fragment in self.fragments:
        #     if fragment is not None:
        #         start_index: int = fragment.end1
        #         end_index: int = fragment.end2
        #
        #         self.ancestor_AA.add_chromosome(self.get_chromosome_using_start_gene(start_index, end_index, 2))
        #         self.ancestor_AA.add_chromosome(self.get_chromosome_using_start_gene(
        #             start_index + self.gene_number * 2, end_index, 2))
        #
        # self.ancestor_A = Genome(list())
        #
        # for fragment in self.fragments:
        #     if fragment is not None:
        #         start_index: int = fragment.end1
        #         end_index: int = fragment.end2
        #
        #         self.ancestor_A.add_chromosome(self.get_chromosome_using_start_gene(start_index, end_index, 1))

    def get_priority_count(self, cs_index: int) -> int:
        """
        Gets the priority count for the ChoiceStructure

        Parameters
        ----------
        cs_index
            Index of the ChoiceStructure

        Returns
        -------
        int
            Priority count for the ChoiceStructure index
        """
        ancestor_priority: int
        result: int = 200
        to_replace: int = -1

        froms: List[int] = [self.choice_structures[cs_index]["genome_paths"][i]["head"] for i in range(self.ploidy)]
        tails: List[int] = [self.choice_structures[cs_index]["genome_paths"][i]["tail"] for i in range(self.ploidy)]

        # from_1: int = self.choice_structures[cs_index]["genome_1_path"]["head"]
        # from_2: int = self.choice_structures[cs_index]["genome_2_path"]["head"]
        # tail_1: int = self.choice_structures[cs_index]["genome_1_path"]["tail"]
        # tail_2: int = self.choice_structures[cs_index]["genome_2_path"]["tail"]

        rng: Random = Random()

        if self.replace == 0:
            rng.seed()
            to_replace = rng.randint(0, 1)
        elif self.replace <= self.ploidy:
            to_replace = self.replace - 1

        # elif self.replace == 1:
        #     to_replace = 0
        # elif self.replace == 2:
        #     to_replace = 1

        for i in range(self.ploidy):
            ancestor_priority = self.calculate_case(cs_index, froms[i], tails[i])

            if ancestor_priority == 0 or ancestor_priority == 1:
                self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(froms[i], tails[i])

                return ancestor_priority

            if ancestor_priority < result:
                result = ancestor_priority
                self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(froms[i], tails[i])

            if ancestor_priority == result and to_replace == 1:
                result = ancestor_priority
                self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(froms[i], tails[i])

            to_replace = 1 - to_replace

        # ancestor_priority = self.calculate_case(cs_index, from_1, tail_1)
        #
        # if ancestor_priority == 0 or ancestor_priority == 1:
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_1, tail_1)
        #
        #     return ancestor_priority
        #
        # if ancestor_priority < result:
        #     result = ancestor_priority
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_1, tail_1)
        #
        # if ancestor_priority == result and to_replace == 1:
        #     result = ancestor_priority
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_1, tail_1)
        #
        # to_replace = 1 - to_replace
        #
        # ancestor_priority = self.calculate_case(cs_index, from_2, tail_2)
        #
        # if ancestor_priority == 1:
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_2, tail_2)
        #
        #     return ancestor_priority
        #
        # if ancestor_priority < result:
        #     result = ancestor_priority
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_2, tail_2)
        #
        # elif ancestor_priority == result and to_replace == 1:
        #     result = ancestor_priority
        #     self.choice_structures[cs_index]["gray_edge"] = PGMPath.create_pgm_path(from_2, tail_2)

        return result

    def get_created_fragment(self, ancestor_1: int, ancestor_2: int) -> List[PGMFragment]:
        """
        Creates a fragment from the given ancestor node indices

        Parameters
        ----------
        ancestor_1
            First ancestor node index
        ancestor_2
            Second ancestor node index

        Returns
        -------
        List[PGMFragment]
            New fragment combination of the given ancestor nodes
        """
        node_1: int = ancestor_1
        node_2: int = ancestor_2
        result: List[PGMFragment] = list()

        if node_1 > 0 and node_2 > 0:
            for i in range(self.ploidy - 1, 0, -1):
                if node_1 == ancestor_1 and node_1 > self.gene_number * i * 2:
                    node_1 -= self.gene_number * i * 2

                if node_2 == ancestor_2 and node_2 > self.gene_number * i * 2:
                    node_2 -= self.gene_number * i * 2

            # if node_1 > self.gene_number * 2:
            #     node_1 -= self.gene_number * 2
            #
            # if node_2 > self.gene_number * 2:
            #     node_2 -= self.gene_number * 2

            result.append(PGMFragment.from_fragment(self.fragments[node_1]))
            result.append(PGMFragment.from_fragment(self.fragments[node_2]))
            result.append(combine(PGMPath.create_pgm_path(node_1, node_2), self.fragments[node_1],
                                  self.fragments[node_2]))

        return result

    def get_new_fragment_list(self, created_frags: List[PGMFragment]):
        """
        Get new fragment list

        Parameters
        ----------
        created_frags
            Created fragments
        """
        if len(created_frags) != 0:
            self.fragments[created_frags[0].end1] = None
            self.fragments[created_frags[0].end2] = None
            self.fragments[created_frags[1].end1] = None
            self.fragments[created_frags[1].end2] = None
            self.fragments[created_frags[2].end1] = created_frags[2]
            self.fragments[created_frags[2].end2] = PGMFragment(created_frags[2].end2, created_frags[2].end1)

    def get_new_choice_structure(self, choice_structure_index: int, index_from: int, tail: int) -> \
            List[Optional[Dict[str, Any]]]:
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
        List[Optional[Dict[str, Any]]]
            New ChoiceStructure
        """
        froms: List[int] = [index_from]
        tails: List[int] = [tail]

        for i in range(self.ploidy - 1, 0, -1):
            if froms[0] == index_from and froms[0] > self.gene_number * i * 2:
                froms[0] -= self.gene_number * i * 2

            if tails[0] == tail and tails[0] > self.gene_number * i * 2:
                tails[0] -= self.gene_number * i * 2

        if froms[0] != choice_structure_index + 1:
            raise Exception("Wrong choice_structure_index: " + str(choice_structure_index))

        # f1: int = index_from
        # t1: int = tail
        #
        # if index_from > self.gene_number * 2:
        #     f1 = index_from - self.gene_number * 2
        #
        # if tail > self.gene_number * 2:
        #     t1 = tail - self.gene_number * 2
        #
        # if f1 != choice_structure_index + 1:
        #     raise Exception("Wrong choice_structure_index: " + str(choice_structure_index))

        paths_x_1: List[Dict[str, int]] = list()
        paths_x_2: List[Dict[str, int]] = list()

        for i in range(self.ploidy):
            paths_x_1.append(self.choice_structures[choice_structure_index]["genome_paths"][i])
            paths_x_2.append(self.choice_structures[tails[0] - 1]["genome_paths"][i])

        # path_1_1: Dict[str, int] = self.choice_structures[choice_structure_index]["genome_1_path"]
        # path_1_2: Dict[str, int] = self.choice_structures[t1 - 1]["genome_1_path"]
        # path_2_1: Dict[str, int] = self.choice_structures[choice_structure_index]["genome_2_path"]
        # path_2_2: Dict[str, int] = self.choice_structures[t1 - 1]["genome_2_path"]
        # path_3_1: Dict[str, int] = self.choice_structures[choice_structure_index]["genome_3_path"]
        # path_3_2: Dict[str, int] = self.choice_structures[t1 - 1]["genome_3_path"]

        l_paths: List[Dict[str, int]] = [PGMPath.create_pgm_path(froms[0], tails[0])]

        for i in range(1, self.ploidy):
            if froms[i - 1] > self.gene_number * i * 2:
                froms.append(froms[i - 1])
            else:
                froms.append(index_from + self.gene_number * i * 2)

            if tails[i - 1] > self.gene_number * i * 2:
                tails.append(tails[i - 1])
            else:
                tails.append(tail + self.gene_number * i * 2)

            l_paths.append(PGMPath.create_pgm_path(froms[i], tails[i]))

        # path_l1: Dict[str, int] = PGMPath.create_pgm_path(index_from, tail)
        #
        # f2: int
        #
        # if index_from > self.gene_number * 2:
        #     f2 = f1
        # else:
        #     f2 = index_from + (self.gene_number * 2)
        #
        # t2: int
        #
        # if tail > self.gene_number * 2:
        #     t2 = t1
        # else:
        #     t2 = tail + (self.gene_number * 2)
        #
        # path_l2: Dict[str, int] = PGMPath.create_pgm_path(f2, t2)

        new_paths: List[Optional[Dict[str, int]]] = [None for _ in range(len(paths_x_1))]

        for i in range(self.ploidy - 1, 0, -1):
            for j in range(self.ploidy):
                if index_from <= self.gene_number * i * 2 and tail <= self.gene_number * i * 2:
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[j], l_paths[j])

                if index_from > self.gene_number * i * 2 and tail > self.gene_number * i * 2:
                    circ_j: int = 0 if j == len(l_paths) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[j], l_paths[circ_j])

                if index_from <= self.gene_number * i * 2 < tail:
                    circ_j: int = 0 if j == len(paths_x_2) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[circ_j], l_paths[j])

                if index_from > self.gene_number * i * 2 >= tail:
                    circ_j: int = 0 if j == len(paths_x_1) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[circ_j], paths_x_2[j], l_paths[j])

        # new_path1: Optional[Dict[str, int]] = None
        # new_path2: Optional[Dict[str, int]] = None
        #
        # if index_from <= self.gene_number * 2 and tail <= self.gene_number * 2:
        #     new_path1 = PGMPath.connect(path_1_1, path_1_2, path_l1)
        #     new_path2 = PGMPath.connect(path_2_1, path_2_2, path_l2)
        #
        # if index_from > self.gene_number * 2 and tail > self.gene_number * 2:
        #     new_path1 = PGMPath.connect(path_1_1, path_1_2, path_l2)
        #     new_path2 = PGMPath.connect(path_2_1, path_2_2, path_l1)
        #
        # if index_from <= self.gene_number * 2 < tail:
        #     new_path1 = PGMPath.connect(path_1_1, path_2_2, path_l1)
        #     new_path2 = PGMPath.connect(path_2_1, path_1_2, path_l2)
        #
        # if index_from > self.gene_number * 2 >= tail:
        #     new_path1 = PGMPath.connect(path_2_1, path_1_2, path_l1)
        #     new_path2 = PGMPath.connect(path_1_1, path_2_2, path_l2)
        #
        # path_l3: Dict[str, int] = PGMPath.create_pgm_path(f1, t1)
        # new_path3: Dict[str, int] = PGMPath.connect(path_3_1, path_3_2, path_l3)

        temp: List[Optional[Dict[str, Any]]] = [None for _ in range(self.ploidy * 4)]

        for path in new_paths:
            if not self.is_a_cycle(path, index_from, tail):
                temp = self.get_new_choice_structure_based_on_path(path, temp)

        return temp[:len(temp) - temp.count(None)]

        # temp: List[Optional[Dict[str, Any]]] = [None for _ in range(6)]
        #
        # if not self.is_a_cycle(new_path1, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path1, temp, None, 2)
        # if not self.is_a_cycle(new_path2, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path2, temp, None, 2)
        # if not self.is_a_cycle(new_path3, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path3, temp, None, 1)
        #
        # new_choice_structure_number: int = 0
        #
        # for choice_structure in temp:
        #     if choice_structure is not None:
        #         new_choice_structure_number += 1
        #
        # return temp[:new_choice_structure_number]

    def is_a_cycle(self, path1: Dict[str, int], f2: int, t2: int) -> bool:
        """
        Tests if the given path and end nodes form a cycle

        Parameters
        ----------
        path1
            First path to check
        f2
            First end of second path
        t2
            Last end of second path

        Returns
        -------
        bool
            True/False whether it forms a cycle
        """
        f1: int = path1["head"]
        t1: int = path1["tail"]
        f: int = f2
        t: int = t2

        for i in range(self.ploidy - 1, 0, -1):
            if f1 == path1["head"] and f1 > self.gene_number * i * 2:
                f1 -= self.gene_number * i * 2

            if t1 == path1["tail"] and t1 > self.gene_number * i * 2:
                t1 -= self.gene_number * i * 2

            if f == f2 and f > self.gene_number * i * 2:
                f -= self.gene_number * i * 2

            if t == t2 and t > self.gene_number * i * 2:
                t -= self.gene_number * i * 2

        # if f1 > self.gene_number * 2:
        #     f1 -= self.gene_number * 2
        #
        # if t1 > self.gene_number * 2:
        #     t1 -= self.gene_number * 2
        #
        # f: int = f2
        # t: int = t2
        #
        # if f > self.gene_number * 2:
        #     f -= self.gene_number * 2
        #
        # if t > self.gene_number * 2:
        #     t -= self.gene_number * 2

        return (f1 == f and t1 == t) or (f1 == t and t1 == f)

    def get_new_choice_structure_based_on_path(self, new_path1: Dict[str, int],
                                               temp: List[Optional[Dict[str, Any]]],
                                               new_choice_structures: Optional[List[Dict[str, Any]]] = None) \
            -> List[Dict[str, Any]]:
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

        Returns
        -------
        List[Dict[str, Any]]
            temp, modified with the new ChoiceStructure
        """
        from_1: int = new_path1["head"]
        from_small: int = -10000
        tail_1: int = new_path1["tail"]
        tail_small: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if from_small == -10000 and from_1 > self.gene_number * i * 2:
                from_small = from_1 - self.gene_number * i * 2

            if tail_small == -10000 and tail_1 > self.gene_number * i * 2:
                tail_small = tail_1 - self.gene_number * i * 2

        if from_small == -10000:
            from_small = from_1

        if tail_small == -10000:
            tail_small = tail_1

        # if from_1 > self.gene_number * 2:
        #     from_small = from_1 - (self.gene_number * 2)
        # else:
        #     from_small = from_1
        #
        # tail_1: int = new_path1["tail"]
        # tail_small: int
        #
        # if tail_1 > self.gene_number * 2:
        #     tail_small = tail_1 - (self.gene_number * 2)
        # else:
        #     tail_small = tail_1

        temp_index: int = len(temp) - temp.count(None)

        if from_1 > 0:
            small_index: int = check_temp_list(temp, from_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:
                    for choice_structure in new_choice_structures:
                        if choice_structure["index_from"] == from_small:
                            temp[temp_index] = ChoiceStructure.create_cs(choice_structure)
                            ChoiceStructure.set_new_path(temp[temp_index], new_path1, self.ploidy, self.gene_number)
                            find_now = True
                            break

                if not find_now and self.choice_structures[from_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure.create_cs(self.choice_structures[from_small - 1])
                    ChoiceStructure.set_new_path(temp[temp_index], new_path1, self.ploidy, self.gene_number)
                    temp_index += 1
            else:
                ChoiceStructure.set_new_path(temp[small_index], new_path1, self.ploidy, self.gene_number)

        if tail_1 > 0:
            np1: Dict[str, int] = PGMPath.create_pgm_path(tail_1, from_1)
            small_index: int = check_temp_list(temp, tail_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:  # LA2 case
                    for choice_structure in new_choice_structures:
                        if choice_structure["index_from"] == tail_small:
                            temp[temp_index] = ChoiceStructure.create_cs(choice_structure)
                            ChoiceStructure.set_new_path(temp[temp_index], np1, self.ploidy, self.gene_number)
                            find_now = True
                            break

                if not find_now and self.choice_structures[tail_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure.create_cs(self.choice_structures[tail_small - 1])
                    ChoiceStructure.set_new_path(temp[temp_index], np1, self.ploidy, self.gene_number)
                    temp_index += 1
            else:
                ChoiceStructure.set_new_path(temp[small_index], np1, self.ploidy, self.gene_number)

        return temp

    def update_priority(self, cs_index: int):
        """
        Updates priority count

        Parameters
        ----------
        cs_index
            ChoiceStructure index
        """
        if self.choice_structures[cs_index] is not None:
            old_priority: int = self.choice_structures[cs_index]["priority"]
            old_position: int = self.choice_structures[cs_index]["position"]
            priority: int = self.get_priority_count(cs_index)

            if priority != old_priority:
                if old_priority < len(self.priorities):
                    self.priorities[old_priority].remove(old_position)

                new_position: int

                if priority < len(self.priorities):
                    new_position = self.priorities[priority].insert(cs_index)
                else:
                    new_position = -1

                self.choice_structures[cs_index]["priority"] = priority
                self.choice_structures[cs_index]["position"] = new_position

    def get_chromosome_using_start_gene(self, start_index: int, end_index: int) -> Chromosome:
        """
        Gets a chromosome using the first gene

        Parameters
        ----------
        start_index
            Beginning index to search from
        end_index
            Ending index to search to

        Returns
        -------
        Chromosome
            Chromosome as a string
        """
        ancestor_chromosome: str = str()

        # if ploidy == 1:
        #     start: str = self.node_str[start_index - 1]
        #
        #     if start.endswith("h"):
        #         ancestor_chromosome = "-" + start[:-1]
        #     else:
        #         ancestor_chromosome = start[:-1]
        #
        #     start_index = get_gene_next_node(start_index)
        #
        #     while True:
        #         next_gene_index: int = self.find_gray_edge_node(start_index, self.gray_edge, 1)
        #
        #         if next_gene_index == -1000:
        #             break
        #
        #         next_gene: str = self.node_str[next_gene_index - 1]
        #         current_start_index: int
        #
        #         if start_index > self.gene_number * 2:
        #             current_start_index = start_index - (self.gene_number * 2)
        #         else:
        #             current_start_index = start_index
        #
        #         if current_start_index == end_index:
        #             break
        #         else:
        #             if next_gene.endswith("h"):
        #                 ancestor_chromosome += " -" + next_gene[:-1]
        #             else:
        #                 ancestor_chromosome += " " + next_gene[:-1]
        #
        #             start_index = get_gene_next_node(next_gene_index)
        # elif ploidy == 2:
        for i in range(self.ploidy - 1, 0, -1):
            if start_index > self.gene_number * i * 2:
                start: str = self.node_str[start_index - self.gene_number * i * 2 - 1]

                if start.endswith("h"):
                    ancestor_chromosome = "-" + start[:-1] + ploidys[i]
                else:
                    ancestor_chromosome = start[:-1] + ploidys[i]

                break

        if ancestor_chromosome == "":
            start: str = self.node_str[start_index - 1]

            if start.endswith("h"):
                ancestor_chromosome = "-" + start[:-1] + "a"
            else:
                ancestor_chromosome = start[:-1] + "a"

        start_index = get_gene_next_node(start_index)

        while True:
            next_gene_index: int = find_gray_edge_node(start_index, self.gray_edge)

            if next_gene_index == -1000:
                break

            next_gene: str = ""

            for i in range(self.ploidy - 1, 0, -1):
                if next_gene_index > self.gene_number * i * 2:
                    next_gene = self.node_str[next_gene_index - self.gene_number * i * 2 - 1]
                    break

            if next_gene == "":
                next_gene = self.node_str[next_gene_index - 1]

            current_start_index: int = -10000

            for i in range(self.ploidy - 1, 0, -1):
                if start_index > self.gene_number * i * 2:
                    current_start_index = start_index - self.gene_number * i * 2
                    break

            if current_start_index == -10000:
                current_start_index = start_index

            if current_start_index == end_index:
                break
            else:
                if next_gene.endswith("h"):
                    ac_changed: bool = False

                    for i in range(self.ploidy - 1, 0, -1):
                        if next_gene_index > self.gene_number * i * 2:
                            ancestor_chromosome += " -" + next_gene[:-1] + ploidys[i]
                            break

                    if not ac_changed:
                        ancestor_chromosome += " -" + next_gene[:-1] + "a"
                else:
                    ac_changed: bool = False

                    for i in range(self.ploidy - 1, 0, -1):
                        if next_gene_index > self.gene_number * i * 2:
                            ancestor_chromosome += " " + next_gene[:-1] + ploidys[i]
                            break

                    if not ac_changed:
                        ancestor_chromosome += " " + next_gene[:-1] + "a"

                start_index = get_gene_next_node(next_gene_index)

        return Chromosome.from_strings(split_at_whitespace(ancestor_chromosome))
    
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
        int
            The index into self.priorities
        """
        if tail < 0:
            return 200

        for i in range(self.ploidy - 1, 0, -1):
            if index_from == tail + self.gene_number * i * 2 or \
               tail == index_from + self.gene_number * i * 2:
                return 200

        cycle_now: int = self.count_cycle_for_edge(choice_structure_index, index_from, tail)

        if cycle_now == 3:
            return 0

        created_fragment: List[PGMFragment] = self.get_created_fragment(index_from, tail)
        self.get_new_fragment_list(created_fragment)
        created_choice_structure: List[Dict[str, Any]] = self.get_new_choice_structure(choice_structure_index,
                                                                                       index_from, tail)

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

    def count_cycle_for_edge(self, choice_structure_index: int, index_from: int, index_tail: int) -> int:
        """
        Counts the cycles for the current edge

        Parameters
        ----------
        choice_structure_index
            Index of the ChoiceStructure
        index_from
            Head node
        index_tail
            Tail node

        Returns
        -------
        int
            Number of cycles
        """
        head: int = index_from
        tail: int = index_tail

        if self.is_circular_chromosome(head, tail):
            return -1

        for i in range(self.ploidy - 1, 0, -1):
            if head == tail + self.gene_number * i * 2 or \
                    tail == head + self.gene_number * i * 2:
                return -1

        result: int = 0

        if tail > 0:
            for i in range(self.ploidy - 1, 0, -1):
                if head == index_from and head > self.gene_number * i * 2:
                    head -= self.gene_number * i * 2

                if tail == index_tail and tail > self.gene_number * i * 2:
                    tail -= self.gene_number * i * 2

            ancestor_choice_structure: ChoiceStructure = self.choice_structures[
                choice_structure_index]

            froms: List[int] = [path["head"] for path in ancestor_choice_structure["genome_paths"]]
            tails: List[int] = [path["tail"] for path in ancestor_choice_structure["genome_paths"]]

            # f1: int = ancestor_choice_structure["genome_1_path"]["head"]
            # t1: int = ancestor_choice_structure["genome_1_path"]["tail"]
            # f2: int = ancestor_choice_structure["genome_2_path"]["head"]
            # t2: int = ancestor_choice_structure["genome_2_path"]["tail"]
            # f3: int = ancestor_choice_structure["genome_3_path"]["head"]
            # t3: int = ancestor_choice_structure["genome_3_path"]["tail"]

            for i in range(len(froms)):
                old_f: int = froms[i]
                old_t: int = tails[i]

                for j in range(self.ploidy - 1, 0, -1):
                    if froms[i] == old_f and froms[i] > self.gene_number * j * 2:
                        froms[i] -= self.gene_number * j * 2

                    if tails[i] == old_t and tails[i] > self.gene_number * j * 2:
                        tails[i] -= self.gene_number * j * 2

                if froms[i] == head and tails[i] == tail:
                    result += 1

            # if f1 > self.gene_number * 2:
            #     f1 -= self.gene_number * 2
            #
            # if t1 > self.gene_number * 2:
            #     t1 -= self.gene_number * 2
            #
            # if f2 > self.gene_number * 2:
            #     f2 -= self.gene_number * 2
            #
            # if t2 > self.gene_number * 2:
            #     t2 -= self.gene_number * 2
            #
            # if f1 == head and t1 == tail:
            #     result += 1
            #
            # if f2 == head and t2 == tail:
            #     result += 1
            #
            # if f3 == head and t3 == tail:
            #     result += 1

        return result

    def count_all_look_ahead_cycles(self, created_choice_structures: List[Dict[str, Any]]) -> List[int]:
        """
        Counts all look ahead cycles for the set of ChoiceStructures

        Parameters
        ----------
        created_choice_structures
            List of ChoiceStructures to count for

        Returns
        -------
        List[int]
            Counts of look ahead cycles for each ChoiceStructure
        """
        result: List[int] = [int(), int(), int()]
        max_cycle: int = 1
        total: int = 0

        for choice_structure in created_choice_structures:
            if choice_structure is not None:
                index_from: int = choice_structure["index_from"]
                old_cycle: int = self.count_cycle_look_ahead(self.choice_structures[index_from - 1])
                new_cycle: int = self.count_cycle_look_ahead(choice_structure)

                if new_cycle > old_cycle and new_cycle > max_cycle:
                    max_cycle = new_cycle

                total += new_cycle - old_cycle

        result[0] = max_cycle
        result[1] = total
        result[2] = self.count_2_steps_look_ahead_cycle(created_choice_structures, max_cycle)

        return result

    def count_2_steps_look_ahead_cycle(self, created_choice_structures: List[Dict[str, Any]], max_cycle: int) -> int:
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
        int
            Number of look ahead cycles
        """
        current_max: int = 1

        for ancestor_choice_structure in created_choice_structures:
            current_count: int = self.count_cycle_look_ahead(ancestor_choice_structure)
            new_count: int

            if current_count == max_cycle:
                froms: List[int] = [path["head"] for path in ancestor_choice_structure["genome_paths"]]
                tails: List[int] = [path["tail"] for path in ancestor_choice_structure["genome_paths"]]

                # from_1: int = ancestor_choice_structure["genome_1_path"]["head"]
                # from_2: int = ancestor_choice_structure["genome_2_path"]["head"]
                # from_3: int = ancestor_choice_structure["genome_3_path"]["head"]
                #
                # tail_1: int = ancestor_choice_structure["genome_1_path"]["tail"]
                # tail_2: int = ancestor_choice_structure["genome_2_path"]["tail"]
                # tail_3: int = ancestor_choice_structure["genome_3_path"]["tail"]

                if max_cycle > 1:
                    current_from: int = -1000
                    current_tail: int = -1000

                    for i in range(len(tails)):
                        for j in range(len(tails)):
                            if i >= j:
                                continue
                            elif tails[i] == tails[j]:
                                current_from = froms[i]
                                current_tail = tails[i]

                    # if tail_1 == tail_2 or tail_1 == tail_3:
                    #     current_from = from_1
                    #     current_tail = tail_1
                    #
                    # if tail_2 == tail_3:
                    #     current_from = from_2
                    #     current_tail = tail_2

                    current_new_choice_structure: List[Dict[str, Any]] = \
                        self.get_new_choice_structure_2_step(current_from,
                                                             current_tail,
                                                             ancestor_choice_structure,
                                                             created_choice_structures)

                    for choice_structure in current_new_choice_structure:
                        current_count = self.count_cycle_look_ahead(choice_structure)

                        if current_count == 3:
                            return current_count

                        if current_max < current_count:
                            current_max = current_count

                else:
                    current_new_choice_structure: List[Dict[str, Any]]

                    for i in range(len(froms)):
                        current_new_choice_structure = self.get_new_choice_structure_2_step(
                            froms[i],
                            tails[i],
                            ancestor_choice_structure,
                            created_choice_structures)

                        for choice_structure in current_new_choice_structure:
                            current_count = self.count_cycle_look_ahead(choice_structure)

                            if current_count == 3:
                                return current_count

                            if current_max < current_count:
                                current_max = current_count

        return current_max

    def get_new_choice_structure_2_step(self,
                                        index_from: int,
                                        tail: int,
                                        ancestor_choice_structure: ChoiceStructure,
                                        new_choice_structures: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Gets a new set of ChoiceStructures in 2 steps

        Parameters
        ----------
        index_from
            The head node
        tail
            Tail node
        ancestor_choice_structure
            Choice structure for ancestor node
        new_choice_structures
            New ChoiceStructures to use

        Returns
        -------
        List[Dict[str, Any]]
            Set of new ChoiceStructures
        """
        if tail <= 0:
            return list()

        paths_x_1: List[Dict[str, int]] = [path for path in ancestor_choice_structure["genome_paths"]]

        # path_1_1: Dict[str, int] = ancestor_choice_structure["genome_1_path"]
        # path_2_1: Dict[str, int] = ancestor_choice_structure["genome_2_path"]
        # path_3_1: Dict[str, int] = ancestor_choice_structure["genome_3_path"]

        # getTheCSWithGandStart(t, allncs)
        new_choice_structure: Optional[Dict[str, Any]] = None
        cont: bool = False

        start: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if tail > self.gene_number * i * 2:
                start = tail - self.gene_number * i * 2
                break

        if start == -10000:
            start = tail

        for choice_structure in new_choice_structures:
            if choice_structure["index_from"] == start:
                new_choice_structure = ChoiceStructure.create_cs(choice_structure)
                cont = True
                break

        if not cont:
            new_choice_structure = ChoiceStructure.create_cs(self.choice_structures[start - 1])

        paths_x_2: List[Dict[str, int]] = [path for path in new_choice_structure["genome_paths"]]

        # path_1_2: Dict[str, int] = new_choice_structure["genome_1_path"]
        # path_2_2: Dict[str, int] = new_choice_structure["genome_2_path"]
        # path_3_2: Dict[str, int] = new_choice_structure["genome_3_path"]

        froms: List[int] = [index_from]
        tails: List[int] = [tail]

        for i in range(1, self.ploidy):
            if froms[i - 1] > self.gene_number * i * 2:
                froms.append(froms[i - 1])
            else:
                froms.append(index_from + self.gene_number * i * 2)

            if tails[i - 1] > self.gene_number * i * 2:
                tails.append(tails[i - 1])
            else:
                tails.append(tail + self.gene_number * i * 2)

        # from_2: int
        # to_2: int
        # from_3: int
        # to_3: int
        #
        # if index_from > self.gene_number * 2:
        #     from_2 = index_from - self.gene_number * 2
        # else:
        #     from_2 = index_from + self.gene_number * 2
        #
        # if tail > self.gene_number * 2:
        #     to_2 = tail - self.gene_number * 2
        # else:
        #     to_2 = tail + self.gene_number * 2
        #
        # if index_from > self.gene_number * 2:
        #     from_3 = from_2
        # else:
        #     from_3 = index_from
        #
        # if tail > self.gene_number * 2:
        #     to_3 = to_2
        # else:
        #     to_3 = tail

        l_paths: List[Dict[str, int]] = [PGMPath.create_pgm_path(f, t)
                                         for f, t in zip(froms, tails)]

        # path_l1: Dict[str, int] = PGMPath.create_pgm_path(index_from, tail)
        # path_l2: Dict[str, int] = PGMPath.create_pgm_path(from_2, to_2)
        # path_l3: Dict[str, int] = PGMPath.create_pgm_path(from_3, to_3)

        new_paths: List[Optional[Dict[str, int]]] = [None for _ in range(len(paths_x_1))]

        for i in range(self.ploidy - 1, 0, -1):
            for j in range(self.ploidy):
                if index_from <= self.gene_number * i * 2 and tail <= self.gene_number * i * 2:
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[j], l_paths[j])

                if index_from > self.gene_number * i * 2 and tail > self.gene_number * i * 2:
                    circ_j: int = 0 if j == len(l_paths) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[j], l_paths[circ_j])

                if index_from <= self.gene_number * i * 2 < tail:
                    circ_j: int = 0 if j == len(paths_x_2) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[j], paths_x_2[circ_j], l_paths[j])

                if index_from > self.gene_number * i * 2 >= tail:
                    circ_j: int = 0 if j == len(paths_x_1) else j
                    new_paths[j] = PGMPath.connect(paths_x_1[circ_j], paths_x_2[j], l_paths[j])

        # new_path1: Optional[Dict[str, int]] = None
        # new_path2: Optional[Dict[str, int]] = None
        #
        # if index_from <= self.gene_number * 2 and tail <= self.gene_number * 2:
        #     new_path1 = PGMPath.connect(path_1_1, path_1_2, path_l1)
        #     new_path2 = PGMPath.connect(path_2_1, path_2_2, path_l2)
        #
        # if index_from > self.gene_number * 2 and tail > self.gene_number * 2:
        #     new_path1 = PGMPath.connect(path_1_1, path_1_2, path_l2)
        #     new_path2 = PGMPath.connect(path_2_1, path_2_2, path_l1)
        #
        # if index_from <= self.gene_number * 2 < tail:
        #     new_path1 = PGMPath.connect(path_1_1, path_2_2, path_l1)
        #     new_path2 = PGMPath.connect(path_2_1, path_1_2, path_l2)
        #
        # if index_from > self.gene_number * 2 >= tail:
        #     new_path1 = PGMPath.connect(path_2_1, path_1_2, path_l1)
        #     new_path2 = PGMPath.connect(path_1_1, path_2_2, path_l2)
        #
        # new_path3: Optional[Dict[str, int]] = PGMPath.connect(path_3_1, path_3_2, path_l3)

        temp: List[Optional[Dict[str, Any]]] = [None for _ in range(self.ploidy * 4)]

        for path in new_paths:
            if not self.is_a_cycle(path, index_from, tail):
                temp = self.get_new_choice_structure_based_on_path(path, temp, new_choice_structures)

        return temp[:len(temp) - temp.count(None)]

        # temp: List[Optional[Dict[str, Any]]] = [None for _ in range(6)]
        #
        # if not self.is_a_cycle(new_path1, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path1, temp, new_choice_structures)
        # if not self.is_a_cycle(new_path2, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path2, temp, new_choice_structures)
        # if not self.is_a_cycle(new_path3, index_from, tail):
        #     temp = self.get_new_choice_structure_based_on_path(new_path3, temp, new_choice_structures)
        #
        # return temp[:len(temp) - temp.count(None)]

    def count_cycle_look_ahead(self, ancestor_choice_structure: Optional[Dict[str, Any]]) -> int:
        """
        Counts cycle look ahead in not-2-steps

        Parameters
        ----------
        ancestor_choice_structure
            The ancestor choice structure

        Returns
        -------
        int
            Count of cycle look aheads
        """
        count: int = 0
        tails: List[int] = list()

        for path in ancestor_choice_structure["genome_paths"]:
            if path["tail"] < 0:
                count += 1

            tails.append(path["tail"])

        if count > 1:
            return count

        # if ancestor_choice_structure["genome_1_path"]["tail"] < 0 and \
        #         ancestor_choice_structure["genome_2_path"]["tail"] < 0 and \
        #         ancestor_choice_structure["genome_3_path"]["tail"] < 0:
        #     return 3
        #
        # if ancestor_choice_structure["genome_1_path"]["tail"] < 0 and \
        #         ancestor_choice_structure["genome_2_path"]["tail"] < 0:
        #     return 2
        #
        # if ancestor_choice_structure["genome_1_path"]["tail"] < 0 and \
        #         ancestor_choice_structure["genome_3_path"]["tail"] < 0:
        #     return 2
        #
        # if ancestor_choice_structure["genome_2_path"]["tail"] < 0 and \
        #         ancestor_choice_structure["genome_3_path"]["tail"] < 0:
        #     return 2
        #
        # tail_1: int
        # tail_2: int
        # tail_3: int

        for i in range(len(tails)):
            for j in range(self.ploidy - 1, 0, -1):
                if tails[i] == ancestor_choice_structure["genome_paths"][i]["tail"] and \
                        tails[i] > self.gene_number * i * 2:
                    tails[i] -= self.gene_number * i * 2

        # if ancestor_choice_structure["genome_1_path"]["tail"] > self.gene_number * 2:
        #     tail_1 = ancestor_choice_structure["genome_1_path"]["tail"] - self.gene_number * 2
        # else:
        #     tail_1 = ancestor_choice_structure["genome_1_path"]["tail"]
        #
        # if ancestor_choice_structure["genome_2_path"]["tail"] > self.gene_number * 2:
        #     tail_2 = ancestor_choice_structure["genome_2_path"]["tail"] - self.gene_number * 2
        # else:
        #     tail_2 = ancestor_choice_structure["genome_2_path"]["tail"]
        #
        # if ancestor_choice_structure["genome_3_path"]["tail"] > self.gene_number * 2:
        #     tail_3 = ancestor_choice_structure["genome_3_path"]["tail"] - self.gene_number * 2
        # else:
        #     tail_3 = ancestor_choice_structure["genome_3_path"]["tail"]

        count = 0

        for i in range(len(tails)):
            for j in range(len(tails)):
                if i >= j:
                    continue
                elif tails[i] == tails[j]:
                    count += 1

        if count <= 1:
            return 1
        else:
            return count

        # if tail_1 == tail_2 and tail_1 == tail_3:
        #     return 3
        #
        # if tail_1 == tail_2 or tail_1 == tail_3 or tail_2 == tail_3:
        #     return 2
        #
        # return 1

    def is_circular_chromosome(self, node_1: int, node_2: int) -> bool:
        """
        Checks if the given nodes form a circular chromosome

        Parameters
        ----------
        node_1
            First node
        node_2
            Second node

        Returns
        -------
        bool
            True/False if it is a circular chromosome
        """
        n1: int = -10000
        n2: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if n1 == -10000 and node_1 > self.gene_number * i * 2:
                n1 = node_1 - self.gene_number * i * 2

            if n2 == -10000 and node_2 > self.gene_number * i * 2:
                n2 = node_2 - self.gene_number * i * 2

        if n1 == -10000:
            n1 = node_1

        if n2 == -10000:
            n2 = node_2

        # if node_1 > self.gene_number * 2:
        #     n1 = node_1 - self.gene_number * 2
        # else:
        #     n1 = node_1
        #
        # if node_2 > self.gene_number * 2:
        #     n2 = node_2 - self.gene_number * 2
        # else:
        #     n2 = node_2

        if 0 < n1 == self.fragments[n1].end1 and \
           self.fragments[n1] is not None and \
           self.fragments[n1].end2 == n2:
            return True

        if 0 < n2 == self.fragments[n2].end1 and \
           self.fragments[n2] is not None and \
           self.fragments[n2].end2 == n1:
            return True

        return False
