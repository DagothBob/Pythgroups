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
            if temp_list[i]["genome_paths"][0]["head"] == f:
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
    node_str
        Nodes in string form
    ideal_ancestor
        A, where d(A1 (+) ... (+) Am, Polyd) is minimal
    replace
        Which genome to replace
    gene_number
        Number of genes in reference genome
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
        self.node_str: List[str] = list()
        self.ideal_ancestor: Genome = Genome(list())

        self.replace: int = replace
        self.gene_number: int = count_gene_number(reference)
        self.ploidy: int = ploidy

        self.priorities: List[Priority] = self.initialize_priorities()

        self.set_nodes(reference)
        self.polyd: List[Optional[Dict[str, int]]] = self.get_pgm_path(polyd)

        self.fragments: List[Optional[PGMFragment]] = [  # Fragments are in pairs (trios for n=3?)
            None for _ in range(self.gene_number * 2 + 1)]

        for i in range(int(len(self.fragments) / 2)):
            self.fragments[2 * i + 1] = PGMFragment(2 * i + 1, 2 * i + 2)
            self.fragments[2 * i + 2] = PGMFragment(2 * i + 2, 2 * i + 1)

        self.gray_edge: List[Optional[Dict[str, int]]] = [None for _ in range(len(self.polyd))]
        self.gray_edge_index: int = 0

        self.choice_structures: List[Optional[Dict[str, Any]]] = [
            None for _ in range(self.gene_number * 2)]  # Scaling for ploidy is in genome_paths

        cs_index: int = 0

        for i in range(1, self.gene_number * 2 + 1):  # self.polyd is 1-indexed
            self.choice_structures[cs_index] = ChoiceStructure.create_cs(ploidy=self.ploidy)
            self.choice_structures[cs_index]["index_from"] = i

            for j in range(self.ploidy):  # For each genome path (0 + i -> gn*2 + i -> gn*4 + i...)
                self.choice_structures[cs_index]["genome_paths"][j] = self.polyd[(self.gene_number * 2 * j) + i]

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
        priority_size: int = self.gene_number * self.ploidy + 1000
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
        self.node_str = [str() for _ in range(self.gene_number * 2)]  # *2 for both nodes

        index: int = 0

        for chromosome in genome.chromosomes:
            for gene in chromosome.genes:
                node_1: str
                node_2: str

                if gene.name[0] == "-":
                    node_1 = str().join([gene.name[1:], "h"])
                    node_2 = str().join([gene.name[1:], "t"])

                    self.node_str[index] = node_2
                    index += 1
                    self.node_str[index] = node_1
                    index += 1
                else:
                    node_1 = str().join([gene.name, "t"])
                    node_2 = str().join([gene.name, "h"])

                    self.node_str[index] = node_1
                    index += 1
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
        path1: List[Optional[Dict[str, int]]] = [None for _ in range((self.gene_number * 2 * self.ploidy) + 1)]

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
                    return (i + 1) + self.gene_number * 2 * (ord(copy_) - 97)

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
            from_1: int = -10000
            tail_1: int = -10000

            for i in range(self.ploidy - 1, 0, -1):
                if from_1 == -10000 and gray_edge["head"] > self.gene_number * 2 * i:
                    from_1 = gray_edge["head"] - self.gene_number * 2 * i

                if tail_1 == -10000 and gray_edge["tail"] > self.gene_number * 2 * i:
                    tail_1 = gray_edge["tail"] - self.gene_number * 2 * i

            if from_1 == -10000:  # Node is already normalized
                from_1 = gray_edge["head"] + self.gene_number * 2 * (self.ploidy - 1)

            if tail_1 == -10000:  # Node is already normalized
                tail_1 = gray_edge["tail"] + self.gene_number * 2 * (self.ploidy - 1)

            self.gray_edge[self.gray_edge_index] = gray_edge
            self.gray_edge_index += 1
            self.gray_edge[self.gray_edge_index] = PGMPath.create_pgm_path(from_1, tail_1)
            self.gray_edge_index += 1

        cs_index: int = gray_edge["head"]

        while cs_index > self.gene_number * 2:
            cs_index -= self.gene_number * 2

        cs_index -= 1

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
            if start_index == -10000 and path_l["head"] > self.gene_number * 2 * i:
                start_index = path_l["head"] - self.gene_number * 2 * i

            if end_index == - 10000 and path_l["tail"] > self.gene_number * 2 * i:
                end_index = path_l["tail"] - self.gene_number * 2 * i

        if start_index == -10000:
            start_index = path_l["head"]

        if end_index == -10000:
            end_index = path_l["tail"]

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
                    if choice_structure["genome_paths"][i]["tail"] > self.gene_number * 2 * j:
                        tails[i] = choice_structure["genome_paths"][i]["tail"] - self.gene_number * 2 * j

                if tails[i] == -10000:
                    tails[i] = choice_structure["genome_paths"][i]["tail"]

            for i in range(self.ploidy):
                if tails[i] > 0:
                    self.update_priority(tails[i] - 1)

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

        for fragment in self.fragments:
            if fragment is not None:
                start_index: int = fragment.end1
                end_index: int = fragment.end2

                for i in range(self.ploidy):
                    self.ideal_ancestor.add_chromosome(
                        self.get_chromosome_using_start_gene(start_index + self.gene_number * i * 2, end_index))

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

        rng: Random = Random()

        if self.replace == 0:
            rng.seed()
            to_replace = rng.randint(0, 1)
        elif self.replace <= self.ploidy:
            to_replace = self.replace - 1

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
                if node_1 == ancestor_1 and node_1 > self.gene_number * 2 * i:
                    node_1 -= self.gene_number * 2 * i

                if node_2 == ancestor_2 and node_2 > self.gene_number * 2 * i:
                    node_2 -= self.gene_number * 2 * i

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
        froms: List[int] = [-10000 for _ in range(self.ploidy)]
        tails: List[int] = [-10000 for _ in range(self.ploidy)]
        froms[0] = index_from
        tails[0] = tail

        norm_t: int = tail
        t_scalar: int = 0

        for i in range(self.ploidy - 1, 0, -1):
            if norm_t > self.gene_number * 2 * i:
                norm_t -= self.gene_number * 2 * i
                t_scalar = i
                break

        norm_f: int = index_from
        f_scalar: int = 0

        for i in range(self.ploidy - 1, 0, -1):
            if index_from > self.gene_number * 2 * i:
                norm_f -= self.gene_number * 2 * i
                f_scalar = i
                break

        paths_x_1: List[Dict[str, int]] = [p for p in self.choice_structures[choice_structure_index]["genome_paths"]]
        paths_x_2: List[Dict[str, int]] = [p for p in self.choice_structures[norm_t - 1]["genome_paths"]]

        current_f_scalar: int = 0
        current_t_scalar: int = 0

        for i in range(self.ploidy):  # List position
            if current_f_scalar == f_scalar:
                current_f_scalar += 1

            if current_t_scalar == t_scalar:
                current_t_scalar += 1

            if i > 0:
                froms[i] = norm_f + self.gene_number * 2 * current_f_scalar
                tails[i] = norm_t + self.gene_number * 2 * current_t_scalar
                current_f_scalar += 1
                current_t_scalar += 1

        l_paths: List[Dict[str, int]] = [PGMPath.create_pgm_path(froms[i], tails[i]) for i in range(self.ploidy)]

        new_paths: List[Optional[Dict[str, int]]] = [None for _ in range(len(l_paths))]

        # "Greater-than" checks are done in descending order from ploidy and
        # "lesser-than" checks are done in ascending order from 0, for best fit
        #
        # This has to be significantly more complex in aliquoting vs. GGH to solve
        # the problem of choosing a best-fit new path
        for i in range(len(new_paths)):
            if norm_f == norm_t:          # Parallel for all
                new_paths[i] = PGMPath.connect(paths_x_1[i], paths_x_2[i], l_paths[i])
            else:
                path_1: Dict[str, int] = paths_x_1[i]
                l_path: Dict[str, int]

                if index_from == norm_f:  # Parallel for froms
                    l_path = l_paths[i]     # First choice
                else:                       # Next choice
                    l_path = l_paths[(i + 1) % len(l_paths)]

                for j in range(len(paths_x_2)):
                    path_attempt: Dict[str, int] = PGMPath.connect(path_1, paths_x_2[j], l_path)

                    if path_attempt is not None:
                        new_paths[i] = path_attempt
                        break

        if new_paths.count(None) > 0:
            raise Exception("New paths were not created successfully for new ChoiceStructure.")

        temp: List[Optional[Dict[str, Any]]] = [None for _ in range(self.ploidy * 4)]

        for path in new_paths:
            if not self.is_a_cycle(path, index_from, tail):
                temp = self.get_new_choice_structure_based_on_path(path, temp)

        return temp[:len(temp) - temp.count(None)]

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
            if f1 == path1["head"] and f1 > self.gene_number * 2 * i:
                f1 -= self.gene_number * 2 * i

            if t1 == path1["tail"] and t1 > self.gene_number * 2 * i:
                t1 -= self.gene_number * 2 * i

            if f == f2 and f > self.gene_number * 2 * i:
                f -= self.gene_number * 2 * i

            if t == t2 and t > self.gene_number * 2 * i:
                t -= self.gene_number * 2 * i

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
            if from_small == -10000 and from_1 > self.gene_number * 2 * i:
                from_small = from_1 - self.gene_number * 2 * i

            if tail_small == -10000 and tail_1 > self.gene_number * 2 * i:
                tail_small = tail_1 - self.gene_number * 2 * i

        if from_small == -10000:
            from_small = from_1

        if tail_small == -10000:
            tail_small = tail_1

        temp_index: int = len(temp) - temp.count(None)

        if from_1 > 0:
            small_index: int = check_temp_list(temp, from_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:
                    for choice_structure in new_choice_structures:
                        if choice_structure["index_from"] == from_small:
                            temp[temp_index] = ChoiceStructure.create_cs(choice_structure)
                            ChoiceStructure.set_new_path(
                                temp[temp_index], new_path1, self.ploidy, self.gene_number, True)
                            find_now = True
                            break

                if not find_now and self.choice_structures[from_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure.create_cs(self.choice_structures[from_small - 1])
                    ChoiceStructure.set_new_path(
                        temp[temp_index], new_path1, self.ploidy, self.gene_number, True)
                    temp_index += 1
            else:
                ChoiceStructure.set_new_path(temp[small_index], new_path1, self.ploidy, self.gene_number, True)

        if tail_1 > 0:
            np1: Dict[str, int] = PGMPath.create_pgm_path(tail_1, from_1)
            small_index: int = check_temp_list(temp, tail_small)

            if small_index == -1:
                find_now: bool = False

                if new_choice_structures is not None:  # LA2 case
                    for choice_structure in new_choice_structures:
                        if choice_structure["index_from"] == tail_small:
                            temp[temp_index] = ChoiceStructure.create_cs(choice_structure)
                            ChoiceStructure.set_new_path(temp[temp_index], np1, self.ploidy, self.gene_number, True)
                            find_now = True
                            break

                if not find_now and self.choice_structures[tail_small - 1] is not None:
                    temp[temp_index] = ChoiceStructure.create_cs(self.choice_structures[tail_small - 1])
                    ChoiceStructure.set_new_path(temp[temp_index], np1, self.ploidy, self.gene_number, True)
                    temp_index += 1
            else:
                ChoiceStructure.set_new_path(temp[small_index], np1, self.ploidy, self.gene_number, True)

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

        for i in range(self.ploidy - 1, 0, -1):
            if start_index > self.gene_number * 2 * i:
                start: str = self.node_str[start_index - self.gene_number * 2 * i - 1]

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
                if next_gene_index > self.gene_number * 2 * i:
                    next_gene = self.node_str[next_gene_index - self.gene_number * 2 * i - 1]
                    break

            if next_gene == "":
                next_gene = self.node_str[next_gene_index - 1]

            current_start_index: int = -10000

            for i in range(self.ploidy - 1, 0, -1):
                if start_index > self.gene_number * 2 * i:
                    current_start_index = start_index - self.gene_number * 2 * i
                    break

            if current_start_index == -10000:
                current_start_index = start_index

            if current_start_index == end_index:
                break
            else:
                prefix: str

                if next_gene.endswith("h"):
                    prefix = " -"
                else:
                    prefix = " "

                ac_changed: bool = False

                for i in range(self.ploidy - 1, 0, -1):
                    if next_gene_index > self.gene_number * 2 * i:
                        add: str = prefix + next_gene[:-1] + ploidys[i]
                        ancestor_chromosome += add
                        ac_changed = True
                        break

                if not ac_changed:
                    ancestor_chromosome += prefix + next_gene[:-1] + "a"

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
            if index_from == tail + self.gene_number * 2 * i or \
               tail == index_from + self.gene_number * 2 * i:
                return 200

        cycle_now: int = self.count_cycles_for_edge(choice_structure_index, index_from, tail)

        if cycle_now == 3:
            return 0

        created_fragment: List[PGMFragment] = self.get_created_fragment(index_from, tail)
        self.get_new_fragment_list(created_fragment)
        created_choice_structure: List[Dict[str, Any]] = self.get_new_choice_structure(choice_structure_index,
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

    def count_cycles_for_edge(self, choice_structure_index: int, index_from: int, index_tail: int) -> int:
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
            if head == tail + self.gene_number * 2 * i or \
                    tail == head + self.gene_number * 2 * i:
                return -1

        result: int = 0

        if tail > 0:
            for i in range(self.ploidy - 1, 0, -1):
                if head == index_from and head > self.gene_number * 2 * i:
                    head -= self.gene_number * 2 * i

                if tail == index_tail and tail > self.gene_number * 2 * i:
                    tail -= self.gene_number * 2 * i

            ancestor_choice_structure: ChoiceStructure = self.choice_structures[
                choice_structure_index]

            froms: List[int] = [path["head"] for path in ancestor_choice_structure["genome_paths"]]
            tails: List[int] = [path["tail"] for path in ancestor_choice_structure["genome_paths"]]

            for i in range(len(froms)):
                old_f: int = froms[i]
                old_t: int = tails[i]

                for j in range(self.ploidy - 1, 0, -1):
                    if froms[i] == old_f and froms[i] > self.gene_number * 2 * j:
                        froms[i] -= self.gene_number * 2 * j

                    if tails[i] == old_t and tails[i] > self.gene_number * 2 * j:
                        tails[i] -= self.gene_number * 2 * j

                if froms[i] == head and tails[i] == tail:
                    result += 1

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

                    current_new_choice_structure: List[Dict[str, Any]] = \
                        self.get_new_choice_structure_2_step(current_from,
                                                             current_tail,
                                                             ancestor_choice_structure,
                                                             created_choice_structures)

                    for choice_structure in current_new_choice_structure:
                        current_count = self.count_cycle_look_ahead(choice_structure)

                        if current_count == self.ploidy:
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

        # getTheCSWithGandStart(t, allncs)
        new_choice_structure: Optional[Dict[str, Any]] = None
        cont: bool = False

        start: int = -10000

        for i in range(self.ploidy - 1, 0, -1):
            if tail > self.gene_number * 2 * i:
                start = tail - self.gene_number * 2 * i
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

        froms: List[int] = [-10000 for _ in range(self.ploidy)]
        tails: List[int] = [-10000 for _ in range(self.ploidy)]
        froms[0] = index_from
        tails[0] = tail

        l_paths: List[Dict[str, int]] = [PGMPath.create_pgm_path(froms[0], tails[0])]

        for i in range(1, self.ploidy):
            for j in range(self.ploidy - 1, 0, -1):
                if froms[i] == -10000 and index_from > self.gene_number * 2 * j:
                    froms[i] = index_from - self.gene_number * 2 * j

                if tails[i] == -10000 and tail > self.gene_number * 2 * j:
                    tails[i] = tail - self.gene_number * 2 * j

            if froms[i] == -10000:  # Node is already normalized
                froms[i] = index_from + self.gene_number * 2 * (self.ploidy - i)

            if tails[i] == -10000:  # Node is already normalized
                tails[i] = tail + self.gene_number * 2 * (self.ploidy - i)

            l_paths.append(PGMPath.create_pgm_path(froms[i], tails[i]))

        new_paths: List[Optional[Dict[str, int]]] = [None for _ in range(len(paths_x_1))]

        # "Greater-than" checks are done in descending order from ploidy and
        # "lesser-than" checks are done in ascending order from 0, for best fit
        #
        # This has to be significantly more complex in aliquoting vs. GGH to solve
        # the problem of choosing a best-fit new path
        for i in range(self.ploidy - 1, 0, -1):
            for j in range(1, self.ploidy):
                for k in range(len(new_paths)):
                    if new_paths[k] is not None:
                        continue

                    if index_from <= self.gene_number * 2 * j and tail <= self.gene_number * 2 * j:
                        new_path: Dict[str, int] = PGMPath.connect(paths_x_1[k], paths_x_2[k], l_paths[k])

                        if new_path is not None:
                            new_paths[k] = new_path

                    if index_from > self.gene_number * 2 * i and tail > self.gene_number * 2 * i:
                        for p in range(1, len(l_paths)):
                            new_path: Dict[str, int] = PGMPath.connect(
                                paths_x_1[k], paths_x_2[k], l_paths[(k + p) % len(l_paths)])

                            if new_path is not None:
                                new_paths[k] = new_path
                                break

                    if index_from <= self.gene_number * 2 * j and tail > self.gene_number * 2 * i:
                        for p in range(1, len(paths_x_2)):
                            new_path: Dict[str, int] = PGMPath.connect(
                                paths_x_1[k], paths_x_2[(k + p) % len(paths_x_2)], l_paths[k])

                            if new_path is not None:
                                new_paths[k] = new_path
                                break

                    if index_from > self.gene_number * 2 * i and tail <= self.gene_number * 2 * j:
                        for p in range(1, len(paths_x_1)):
                            new_path: Dict[str, int] = PGMPath.connect(
                                paths_x_1[(k + p) % len(paths_x_1)], paths_x_2[k], l_paths[k])

                            if new_path is not None:
                                new_paths[k] = new_path
                                break

        temp: List[Optional[Dict[str, Any]]] = [None for _ in range(self.ploidy * 4)]

        for path in new_paths:
            if not self.is_a_cycle(path, index_from, tail):
                temp = self.get_new_choice_structure_based_on_path(path, temp, new_choice_structures)

        return temp[:len(temp) - temp.count(None)]

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

        for i in range(len(tails)):
            for j in range(self.ploidy - 1, 0, -1):
                if tails[i] == ancestor_choice_structure["genome_paths"][i]["tail"] and \
                        tails[i] > self.gene_number * 2 * j:
                    tails[i] -= self.gene_number * 2 * j
                    break

        max_equal: int = 0

        for i in range(len(tails)):
            equal_vec: List[int] = list()

            for j in range(len(tails)):
                if i == j:
                    continue
                elif tails[i] == tails[j]:
                    if i not in equal_vec:
                        equal_vec.append(i)

                    if j not in equal_vec:
                        equal_vec.append(j)

            count = len(equal_vec)

            if count > max_equal:
                max_equal = count

        if max_equal > 1:
            return max_equal
        else:
            return 1

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
            if n1 == -10000 and node_1 > self.gene_number * 2 * i:
                n1 = node_1 - self.gene_number * 2 * i

            if n2 == -10000 and node_2 > self.gene_number * 2 * i:
                n2 = node_2 - self.gene_number * 2 * i

        if n1 == -10000:
            n1 = node_1

        if n2 == -10000:
            n2 = node_2

        if 0 < n1 == self.fragments[n1].end1 and \
           self.fragments[n1] is not None and \
           self.fragments[n1].end2 == n2:
            return True

        if 0 < n2 == self.fragments[n2].end1 and \
           self.fragments[n2] is not None and \
           self.fragments[n2].end2 == n1:
            return True

        return False
