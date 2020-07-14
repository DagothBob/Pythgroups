from copy import deepcopy
from typing import List, Optional

from DCJOperation import DCJOperation, OperationTypes, FusionSubtypes, TranslocationSubtypes
from GeneNode import GeneNode
from GenomeInString import GenomeInString

"""                               
 Double-Cut-and-Join Rearrangement class                        
                                                             
 DCJ distance is the minimum number of DCJ operations -      
 inversion, translocation, fission, and fusion - required    
 to transform a genome into another. Also can be expressed   
 in an adjacency graph as the number of adjacencies, minus   
 the number of cycles present.                               
                                                             
 Input: Two adjacencies                                      
 Output: Two new adjacencies created by recombining the four 
         involved extremities                                
                                                             
 Based on DCJOperations.java from C.Zheng & D.Sankoff (2011) 
                                                             
 Author: Holger Jensen                                       
"""


def insert_character(s: str, i: int, c: str) -> str:
    """
    Utility function for modifying strings in a similar way to lists.

    Parameters
    ----------
    s
        String to modify
    i
        Index to insert the character at
    c
        Character to insert

    Returns
    -------
    str
        New string with the inserted character
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


class DCJRearrangement:
    """
    Attributes
    ----------
    genome1 : GenomeInString
        Genome 1 as a list of strings
    genome2 : GenomeInString
        Genome 2 as a list of strings
    gene_count : int
        How many genes are in both genomes (must be the same)
    max_chromosome_number : int
        Highest number of chromosomes in either genome
    min_chromosome_number : int
        Lower number of chromosomes in either genome
    gene_nodes_1 : List[List[int]]
        First gene node for performing DCJ operations
    gene_nodes_2 : List[List[int]]
        Second gene node for performing DCJ operations
    chromosomes_for_gene_node_1 : List[List[int]]
        Chromosome for first gene node
    chromosomes_for_gene_node_2 : List[List[int]]
        Chromosome for second gene node
    gene_node_as_string : List[str]
        String representation for both gene nodes
    touched_chromosome : List[int]
        Current chromosome being updated
    """

    def __init__(self, genome1: List[str], genome2: List[str]):
        """
        Constructor

        Parameters
        ----------
        genome1
            First genome for DCJ operations
        genome2
            Second genome for DCJ operations
        """
        self.genome1: GenomeInString = GenomeInString(genome1)
        self.genome2: GenomeInString = GenomeInString(genome2)
        self.gene_count: int = int()

        # Parameters for choosing an appropriate operation
        self.max_chromosome_number: int = int()
        self.min_chromosome_number: int = int()

        # For each row, use GeneNodeAttributes to index
        self.gene_nodes_1: Optional[List[Optional[GeneNode]]] = None
        self.gene_nodes_2: Optional[List[Optional[GeneNode]]] = None

        self.chromosomes_for_gene_node_1: Optional[List[Optional[List[int]]]] = None
        self.chromosomes_for_gene_node_2: Optional[List[Optional[List[int]]]] = None
        self.gene_node_as_string: Optional[List[List[str]]] = None
        self.touched_chromosome: Optional[List[int]] = None

    def initial_value(self):
        """
        Initializer
        """
        gene_count_1: int = 0
        gene_count_2: int = 0

        # Split genome_1 into genes and save the number of them
        for i in range(len(self.genome1.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome1.chromosomes[i])
            gene_count_1 += len(genes)

        # Split genome_2 into genes and save the number of them
        for i in range(len(self.genome2.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome2.chromosomes[i])
            gene_count_2 += len(genes)

        self.gene_count = gene_count_1
        self.gene_node_as_string = [[str(), str()] for _ in range((self.gene_count * 2))]
        self.gene_nodes_1 = [GeneNode() for _ in range((self.gene_count * 2))]
        self.gene_nodes_2 = [GeneNode() for _ in range((self.gene_count * 2))]
        self.chromosomes_for_gene_node_1 = [[int()] for _ in range(len(self.genome1.chromosomes))]
        self.chromosomes_for_gene_node_2 = [[int()] for _ in range(len(self.genome2.chromosomes))]

        index: int = 0

        # Split each chromosome of genome_1 into genes and save genes into gene_nodes_1 and chromosomes_for_gene_node_1
        for i in range(len(self.genome1.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome1.chromosomes[i])
            self.chromosomes_for_gene_node_1[i] = [int() for _ in range((len(genes) * 2))]

            for j in range(len(genes)):
                sign: str = "+"
                genome_name: str = genes[j]

                if genome_name[0] == "-":
                    sign = "-"
                    genome_name = genes[j][1:]

                if sign == "+":
                    self.gene_node_as_string[index][0] = genome_name
                    self.gene_node_as_string[index][1] = "t"
                    index += 1
                    self.gene_node_as_string[index][0] = genome_name
                    self.gene_node_as_string[index][1] = "h"
                else:
                    self.gene_node_as_string[index][0] = genome_name
                    self.gene_node_as_string[index][1] = "h"
                    index += 1
                    self.gene_node_as_string[index][0] = genome_name
                    self.gene_node_as_string[index][1] = "t"

                index += 1

                if j == 0:
                    self.gene_nodes_1[index - 2].adjacency = -1
                else:
                    self.gene_nodes_1[index - 2].adjacency = index - 3

                if j == len(genes) - 1:
                    self.gene_nodes_1[index - 1].adjacency = -1
                else:
                    self.gene_nodes_1[index - 1].adjacency = index

                self.gene_nodes_1[index - 2].chromosome_index = i
                self.gene_nodes_1[index - 1].chromosome_index = i
                self.gene_nodes_1[index - 2].chromosome_position = j * 2
                self.gene_nodes_1[index - 1].chromosome_position = (j * 2) + 1
                self.chromosomes_for_gene_node_1[i][j * 2] = index - 2
                self.chromosomes_for_gene_node_1[i][(j * 2) + 1] = index - 1

        # Split each chromosome of genome_2 into genes and save genes into gene_nodes_2 and chromosomes_for_gene_node_2
        for i in range(len(self.genome2.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome2.chromosomes[i])
            pre_node_index: int = -1
            self.chromosomes_for_gene_node_2[i] = [int() for _ in range((len(genes) * 2))]

            for j in range(len(genes)):
                sign: str = "+"
                genome_name: str = genes[j]

                if genome_name[0] == "-":
                    sign = "-"
                    genome_name = genes[j][1:]

                genome_node_string1: str = insert_character(genome_name, len(genome_name), "t")
                genome_node_string2: str = insert_character(genome_name, len(genome_name), "h")

                if sign == "-":
                    genome_node_string1 = insert_character(genome_name, len(genome_name), "h")
                    genome_node_string2 = insert_character(genome_name, len(genome_name), "t")

                genome_node_index1: int = self.get_node_index(genome_node_string1)
                genome_node_index2: int = self.get_node_index(genome_node_string2)

                self.gene_nodes_2[genome_node_index1].chromosome_index = i
                self.gene_nodes_2[genome_node_index2].chromosome_index = i
                self.gene_nodes_2[genome_node_index1].chromosome_position = j * 2
                self.gene_nodes_2[genome_node_index2].chromosome_position = (j * 2) + 1
                self.chromosomes_for_gene_node_2[i][j * 2] = genome_node_index1
                self.chromosomes_for_gene_node_2[i][(j * 2) + 1] = genome_node_index2

                self.gene_nodes_2[genome_node_index1].adjacency = pre_node_index

                if j != 0:
                    self.gene_nodes_2[pre_node_index].adjacency = genome_node_index1

                pre_node_index = genome_node_index2

            self.gene_nodes_2[pre_node_index].adjacency = -1

    def get_node_index(self, g_string: str) -> int:  # Renamed from findNodeIndex() for clarity
        """
        Get the node index for given gene pair

        Parameters
        ----------
        g_string
            Genome as a string

        Returns
        -------
        int
            Node index
        """
        for i in range(len(self.gene_node_as_string)):
            if g_string == self.gene_node_as_string[i][0] + self.gene_node_as_string[i][1]:
                return i

        return -1

    def get_genome_in_string(self, node_list: List[List[int]]) -> GenomeInString:
        """
        Creates a GenomeInString from the given list of nodes

        Parameters
        ----------
        node_list
            List of nodes

        Returns
        -------
        GenomeInString
            New GenomeInString
        """
        chromosome: List[str] = [str() for _ in range(len(node_list))]

        for i in range(len(node_list)):
            chromosome[i] = ""

            for j in range(len(node_list[i]) // 2):
                genome1: str = self.gene_node_as_string[node_list[i][j * 2]][0]
                tail_or_head1: str = self.gene_node_as_string[node_list[i][j * 2]][1]
                genome2: str = self.gene_node_as_string[node_list[i][(j * 2) + 1]][0]
                tail_or_head2: str = self.gene_node_as_string[node_list[i][(j * 2) + 1]][1]

                if genome1 == genome2 and tail_or_head1 == "t" and tail_or_head2 == "h":
                    chromosome[i] = chromosome[i] + " " + genome1
                else:
                    if genome1 == genome2 and tail_or_head1 == "h" and tail_or_head2 == "t":
                        chromosome[i] = chromosome[i] + " -" + genome1

        return GenomeInString(chromosome)

    def get_result(self, min_chromosome: int, max_chromosome: int, which_chromosome: int, types_of_operation: List[int],
                   number_of_operations: int) -> List[Optional[GenomeInString]]:
        """
        Performs a DCJ operation

        Parameters
        ----------
        min_chromosome
            Minimum chromosome count
        max_chromosome
            Maximum chromosome count
        which_chromosome
            Which chromosome to perform operation on
        types_of_operation
            Types of operation to perform
        number_of_operations
            How many operations to perform

        Returns
        -------
        [GenomeInString]
            New list of GenomeInStrings after operations performed
        """
        all_steps: List[Optional[GenomeInString]] = [None for _ in range(number_of_operations)]
        step_index: int = 0
        current_operation: int = 0

        self.touched_chromosome = [int() for _ in range(len(self.chromosomes_for_gene_node_1))]
        self.min_chromosome_number = min_chromosome
        self.max_chromosome_number = max_chromosome

        more: bool = True

        while current_operation < number_of_operations and more:
            more = False
            operation: DCJOperation = self.expand_all_operations(types_of_operation, which_chromosome)

            if operation.chromosome_1 != -1:
                if operation.operation_type == OperationTypes.INVERSION:
                    print("Inversion: " + str(operation.chromosome_1) + "\n")
                elif operation.operation_type == OperationTypes.TRANSLOCATION:
                    print("Translocation: " + str(operation.chromosome_1) + " " + str(operation.chromosome_2) + "\n")
                elif operation.operation_type == OperationTypes.FISSION:
                    print("Fission: " + str(operation.chromosome_1) + " " + str(operation.chromosome_2) + "\n")
                elif operation.operation_type == OperationTypes.FUSION:
                    print("Fusion: " + str(operation.chromosome_1) + " " + str(operation.chromosome_2) + "\n")

                self.update_all_values(operation, which_chromosome)

                all_steps[step_index] = self.get_genome_in_string(self.chromosomes_for_gene_node_1)
                step_index += 1
                current_operation += 1
                more = True

        return deepcopy(all_steps)[:step_index]

    def update_all_values(self, operation: DCJOperation, which_chromosome: int):
        """
        Updates instance attributes based on the operation

        Parameters
        ----------
        operation
            Which operation is being run
        which_chromosome
            Which chromosome operating on
        """
        if which_chromosome == -1:
            if operation.chromosome_1 != -1:
                self.touched_chromosome[operation.chromosome_1] = 1

            if operation.chromosome_2 != -1:
                self.touched_chromosome[operation.chromosome_2] = 1

        if operation.operation_type == OperationTypes.INVERSION:
            start_node: int = \
                self.gene_nodes_1[operation.node_3].chromosome_position
            end_node: int = self.gene_nodes_1[operation.node_2].chromosome_position
            chromosome_here: int = operation.chromosome_1
            new_chromosome_here: List[int] = \
                [int() for _ in range(len(self.chromosomes_for_gene_node_1[chromosome_here]))]

            for i in range(start_node):
                new_chromosome_here[i] = self.chromosomes_for_gene_node_1[chromosome_here][i]

            for i in range(start_node, end_node + 1):
                new_chromosome_here[i] = self.chromosomes_for_gene_node_1[chromosome_here][end_node - i + start_node]

            for i in range(end_node + 1, len(new_chromosome_here)):
                new_chromosome_here[i] = self.chromosomes_for_gene_node_1[chromosome_here][i]

            for i in range(start_node, end_node + 1):
                self.gene_nodes_1[new_chromosome_here[i]].chromosome_position = i

            self.chromosomes_for_gene_node_1[chromosome_here] = new_chromosome_here
        elif operation.operation_type == OperationTypes.TRANSLOCATION:
            length_1: int = len(self.chromosomes_for_gene_node_1[operation.chromosome_1])
            length_2: int = len(self.chromosomes_for_gene_node_1[operation.chromosome_2])
            length_a: int = 0

            if operation.node_1 != -1:
                length_a = \
                    self.gene_nodes_1[operation.node_1].chromosome_position + 1

            length_b: int = length_1 - length_a
            new_chromosome1: List[int] = list()
            new_chromosome2: List[int] = list()

            if operation.operation_subtype == TranslocationSubtypes.ADCB:
                length_c: int = 0

                if operation.node_4 != -1:
                    length_c = self.gene_nodes_1[operation.node_4].chromosome_position + 1

                length_d: int = length_2 - length_c
                new_chromosome1 = [int() for _ in range(length_a + length_d)]
                new_chromosome2 = [int() for _ in range(length_c + length_b)]

                for i in range(length_a):
                    new_chromosome1[i] = self.chromosomes_for_gene_node_1[operation.chromosome_1][i]

                for i in range(length_d):
                    new_chromosome1[i + length_a] = \
                        self.chromosomes_for_gene_node_1[operation.chromosome_2][i + length_c]

                for i in range(length_c):
                    new_chromosome2[i] = self.chromosomes_for_gene_node_1[operation.chromosome_2][i]

                for i in range(length_b):
                    new_chromosome2[i + length_c] = \
                        self.chromosomes_for_gene_node_1[operation.chromosome_1][i + length_a]
            elif operation.operation_subtype == TranslocationSubtypes.AnCnBD:
                length_c: int = 0

                if operation.node_2 != -1:
                    length_c = self.gene_nodes_1[operation.node_2].chromosome_position + 1

                length_d: int = length_2 - length_c
                new_chromosome1 = [int() for _ in range(length_a + length_c)]
                new_chromosome2 = [int() for _ in range(length_b + length_d)]

                for i in range(length_a):
                    new_chromosome1[i] = deepcopy(
                        self.chromosomes_for_gene_node_1[operation.chromosome_1][i])

                for i in range(length_c):
                    new_chromosome1[i + length_a] = deepcopy(
                        self.chromosomes_for_gene_node_1[operation.chromosome_2][length_c - i - 1])

                for i in range(length_b):
                    new_chromosome2[i] = deepcopy(
                        self.chromosomes_for_gene_node_1[operation.chromosome_1][length_1 - i - 1])

                for i in range(length_d):
                    new_chromosome2[i + length_b] = deepcopy(
                        self.chromosomes_for_gene_node_1[operation.chromosome_2][i + length_c])

            for i in range(len(new_chromosome1)):
                self.gene_nodes_1[new_chromosome1[i]].chromosome_index = operation.chromosome_1
                self.gene_nodes_1[new_chromosome1[i]].chromosome_position = i

            for i in range(len(new_chromosome2)):
                self.gene_nodes_1[new_chromosome2[i]].chromosome_index = operation.chromosome_2
                self.gene_nodes_1[new_chromosome2[i]].chromosome_position = i

            self.chromosomes_for_gene_node_1[operation.chromosome_1] = new_chromosome1
            self.chromosomes_for_gene_node_1[operation.chromosome_2] = new_chromosome2
        elif operation.operation_type == OperationTypes.FISSION:
            chromosome: int = operation.chromosome_1
            length: int = len(self.chromosomes_for_gene_node_1[chromosome])
            length_1: int = self.gene_nodes_1[operation.node_1].chromosome_position + 1
            length_2: int = length - length_1
            new_chromosome1: List[int] = [int() for _ in range(length_1)]
            new_chromosome2: List[int] = [int() for _ in range(length_2)]

            for i in range(length_1):
                new_chromosome1[i] = self.chromosomes_for_gene_node_1[operation.chromosome_1][i]

            for i in range(length_2):
                new_chromosome2[i] = \
                    self.chromosomes_for_gene_node_1[operation.chromosome_1][length_1 + i]

            chromosome_in_gene_node1_temp: List[List[int]] = \
                [list() for _ in range(len(self.chromosomes_for_gene_node_1) + 1)]

            for i in range(len(self.chromosomes_for_gene_node_1)):
                chromosome_in_gene_node1_temp[i] = self.chromosomes_for_gene_node_1[i]

            chromosome_in_gene_node1_temp[operation.chromosome_1] = new_chromosome1
            chromosome_in_gene_node1_temp[len(self.chromosomes_for_gene_node_1)] = new_chromosome2

            self.chromosomes_for_gene_node_1 = [list() for _ in range(len(chromosome_in_gene_node1_temp))]

            for i in range(len(chromosome_in_gene_node1_temp)):
                self.chromosomes_for_gene_node_1[i] = chromosome_in_gene_node1_temp[i]
        elif operation.operation_type == OperationTypes.FUSION:
            chromosome1: int = operation.chromosome_1
            chromosome2: int = operation.chromosome_2
            length_1: int = len(self.chromosomes_for_gene_node_1[chromosome1])
            length_2: int = len(self.chromosomes_for_gene_node_1[chromosome2])
            length: int = length_1 + length_2
            new_chromosome: List[int] = [int() for _ in range(length)]

            if operation.operation_subtype == FusionSubtypes.TWO_REVERSED_PLUS_ONE:
                for i in range(length_2):
                    new_chromosome[i] = self.chromosomes_for_gene_node_1[chromosome2][length_2 - i]

                new_chromosome[length_2:length_1] = deepcopy(self.chromosomes_for_gene_node_1[chromosome1])
            elif operation.operation_subtype == FusionSubtypes.TWO_PLUS_ONE:
                new_chromosome[:length_2] = deepcopy(self.chromosomes_for_gene_node_1[chromosome2])
                new_chromosome[length_2:length_1] = deepcopy(self.chromosomes_for_gene_node_1[chromosome1])
            elif operation.operation_subtype == FusionSubtypes.ONE_PLUS_TWO:
                new_chromosome[:length_1] = deepcopy(self.chromosomes_for_gene_node_1[chromosome1])
                new_chromosome[length_1:length_2] = deepcopy(self.chromosomes_for_gene_node_1[chromosome2])
            elif operation.operation_subtype == FusionSubtypes.ONE_PLUS_TWO_REVERSED:
                new_chromosome[:length_1] = deepcopy(self.chromosomes_for_gene_node_1[chromosome1])

                for i in range(length_2):
                    new_chromosome[i + length_1] = self.chromosomes_for_gene_node_1[chromosome2][length_2 - i]

            chromosome_small: int = chromosome1

            if chromosome2 < chromosome1:
                chromosome_small = chromosome2

            chromosome_big: int = chromosome2

            if chromosome1 > chromosome2:
                chromosome_big = chromosome1

            self.chromosomes_for_gene_node_1[chromosome_small] = new_chromosome
            self.chromosomes_for_gene_node_1[chromosome_big] = [0]

            for i in range(len(new_chromosome)):
                self.gene_nodes_1[new_chromosome[i]].chromosome_position = i

            for i in range(len(self.gene_nodes_1)):
                if self.gene_nodes_1[i].chromosome_index >= chromosome1 and \
                        self.gene_nodes_1[i].chromosome_index >= chromosome2:
                    self.gene_nodes_1[i].chromosome_index = self.gene_nodes_1[i].chromosome_index - 1

            chromosome_in_gene_node1_temp: List[List[int]] = [
                [int()] for _ in range((len(self.chromosomes_for_gene_node_1) - 1))]

            index: int = 0

            for ints in self.chromosomes_for_gene_node_1:
                if ints != [0]:
                    chromosome_in_gene_node1_temp[index] = ints
                    index += 1

            self.chromosomes_for_gene_node_1 = [[int()] for _ in range(len(chromosome_in_gene_node1_temp))]
            self.chromosomes_for_gene_node_1 = chromosome_in_gene_node1_temp

        if operation.node_1 != -1:
            self.gene_nodes_1[operation.node_1].adjacency = operation.node_2

        if operation.node_2 != -1:
            self.gene_nodes_1[operation.node_2].adjacency = operation.node_1

        if operation.node_3 != -1:
            self.gene_nodes_1[operation.node_3].adjacency = operation.node_4

        if operation.node_4 != -1:
            self.gene_nodes_1[operation.node_4].adjacency = operation.node_3

    # Renamed from findAOperation() for accuracy
    def expand_all_operations(self, types_of_operation: List[int], which_chromosome: int) -> DCJOperation:
        """
        Turn the list of given operation types into a list of operation attributes

        Parameters
        ----------
        types_of_operation
            Operation types to expand
        which_chromosome
            Which chromosome is being operated on

        Returns
        -------
        [int]
            List of operation attributes
        """
        result: DCJOperation = DCJOperation()
        result.chromosome_1 = -1

        for operation in types_of_operation:
            result = self.get_operation_attributes(which_chromosome, operation)

            if result.chromosome_1 != -1:
                return result

        return result

    # Renamed from findATypeOperation() for accuracy
    def get_operation_attributes(self, which_chromosome: int, operation_type: int) -> DCJOperation:
        """
        Gets a list of operation attributes for given operation type

        Parameters
        ----------
        which_chromosome
            Which chromosome being operated on
        operation_type
            Operation type to check

        Returns
        -------
        [int]
            List of operation attributes
        """
        result: DCJOperation = DCJOperation()
        result.chromosome_1 = -1

        if operation_type == OperationTypes.FUSION and \
           len(self.chromosomes_for_gene_node_1) <= self.min_chromosome_number or \
           operation_type == OperationTypes.FISSION and \
           len(self.chromosomes_for_gene_node_1) >= self.max_chromosome_number:

            return result

        if operation_type == OperationTypes.FISSION:
            for i in range(len(self.chromosomes_for_gene_node_2) - 1):
                node1: int = self.chromosomes_for_gene_node_2[i][0]
                node2: int = self.chromosomes_for_gene_node_2[i][len(self.chromosomes_for_gene_node_2[i]) - 1]

                for j in range(i + 1, len(self.chromosomes_for_gene_node_2)):
                    node3: int = self.chromosomes_for_gene_node_2[j][0]
                    node4: int = self.chromosomes_for_gene_node_2[j][len(self.chromosomes_for_gene_node_2[j]) - 1]

                    if self.gene_nodes_1[node1].adjacency == node3 and \
                            ((-1 < which_chromosome == self.gene_nodes_1[node1].chromosome_index) or
                             (which_chromosome < 0 and
                              self.touched_chromosome[self.gene_nodes_1[node1].chromosome_index] == 0)):
                        result.chromosome_1 = self.gene_nodes_1[node1].chromosome_index
                        result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                        result.node_1 = node1
                        result.node_2 = -1
                        result.node_3 = node3
                        result.node_4 = -1
                        result.operation_type = OperationTypes.FISSION
                        result.operation_subtype = 1

                        return result
                    elif self.gene_nodes_1[node1].adjacency == node4 and \
                            ((-1 < which_chromosome == self.gene_nodes_1[node1].chromosome_index) or
                             (which_chromosome < 0 and
                              self.touched_chromosome[self.gene_nodes_1[node1].chromosome_index] == 0)):
                        result.chromosome_1 = self.gene_nodes_1[node1].chromosome_index
                        result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                        result.node_1 = node1
                        result.node_2 = -1
                        result.node_3 = node4
                        result.node_4 = -1
                        result.operation_type = OperationTypes.FISSION
                        result.operation_subtype = 1

                        return result
                    elif self.gene_nodes_1[node2].adjacency == node3 and \
                            ((-1 < which_chromosome == self.gene_nodes_1[node2].chromosome_index) or
                             (which_chromosome < 0 and
                              self.touched_chromosome[self.gene_nodes_1[node2].chromosome_index] == 0)):
                        result.chromosome_1 = self.gene_nodes_1[node2].chromosome_index
                        result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                        result.node_1 = node2
                        result.node_2 = -1
                        result.node_3 = node3
                        result.node_4 = -1
                        result.operation_type = OperationTypes.FISSION
                        result.operation_subtype = 1

                        return result
                    elif self.gene_nodes_1[node2].adjacency == node4 and \
                            ((-1 < which_chromosome == self.gene_nodes_1[node2].chromosome_index) or
                             (which_chromosome < 0 and
                              self.touched_chromosome[self.gene_nodes_1[node2].chromosome_index] == 0)):
                        result.chromosome_1 = self.gene_nodes_1[node2].chromosome_index
                        result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                        result.node_1 = node2
                        result.node_2 = -1
                        result.node_3 = node4
                        result.node_4 = -1
                        result.operation_type = OperationTypes.FISSION
                        result.operation_subtype = 1

                        return result

            for i in range(len(self.chromosomes_for_gene_node_2)):
                node1 = self.chromosomes_for_gene_node_2[i][0]
                node2 = self.chromosomes_for_gene_node_2[i][len(self.chromosomes_for_gene_node_2[i]) - 1]

                chromosome_node1: int = self.gene_nodes_1[node1].chromosome_index
                node3: int = self.gene_nodes_1[node1].adjacency

                if node3 != -1 and (chromosome_node1 == which_chromosome or
                                    which_chromosome < 0 and
                                    self.touched_chromosome[chromosome_node1] == 0):
                    result.chromosome_1 = self.gene_nodes_1[node1].chromosome_index
                    result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                    result.node_1 = node1
                    result.node_2 = -1
                    result.node_3 = node3
                    result.node_4 = -1
                    result.operation_type = OperationTypes.FISSION
                    result.operation_subtype = 1

                    return result

                chromosome_node2: int = self.gene_nodes_1[node2].chromosome_index
                node4: int = self.gene_nodes_1[node2].adjacency

                if node4 != -1 and (chromosome_node2 == which_chromosome or
                                    which_chromosome < 0 and
                                    self.touched_chromosome[chromosome_node2] == 0):
                    result.chromosome_1 = self.gene_nodes_1[node2].chromosome_index
                    result.chromosome_2 = len(self.chromosomes_for_gene_node_1)
                    result.node_1 = node2
                    result.node_2 = -1
                    result.node_3 = node4
                    result.node_4 = -1
                    result.operation_type = OperationTypes.FISSION
                    result.operation_subtype = 1

                    return result
        elif operation_type == OperationTypes.FUSION:
            for i in range(len(self.chromosomes_for_gene_node_1) - 1):
                node1: int = self.chromosomes_for_gene_node_1[i][0]
                node2: int = self.chromosomes_for_gene_node_1[i][len(self.chromosomes_for_gene_node_1[i]) - 1]

                for j in range(i + 1, len(self.chromosomes_for_gene_node_1)):
                    node3: int = self.chromosomes_for_gene_node_1[j][0]
                    node4: int = self.chromosomes_for_gene_node_1[j][len(self.chromosomes_for_gene_node_1[j]) - 1]

                    if (which_chromosome < 0 and self.touched_chromosome[i] == 0 and
                        self.touched_chromosome[j] == 0) or (which_chromosome > -1 and (i == which_chromosome or
                                                                                        j == which_chromosome)):
                        if self.gene_nodes_2[node1].adjacency == node3:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node1
                            result.node_2 = node3
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.TWO_REVERSED_PLUS_ONE

                            return result
                        elif self.gene_nodes_2[node1].adjacency == node4:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node1
                            result.node_2 = node4
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.TWO_PLUS_ONE

                            return result
                        elif self.gene_nodes_2[node2].adjacency == node3:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node2
                            result.node_2 = node3
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.ONE_PLUS_TWO

                            return result
                        elif self.gene_nodes_2[node2].adjacency == node4:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node2
                            result.node_2 = node4
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.ONE_PLUS_TWO_REVERSED

                            return result

            for i in range(len(self.chromosomes_for_gene_node_1) - 1):
                node1: int = self.chromosomes_for_gene_node_1[i][0]
                node2: int = self.chromosomes_for_gene_node_1[i][len(self.chromosomes_for_gene_node_1[i]) - 1]

                for j in range(i + 1, len(self.chromosomes_for_gene_node_1)):
                    node3: int = self.chromosomes_for_gene_node_1[j][0]
                    node4: int = self.chromosomes_for_gene_node_1[j][len(self.chromosomes_for_gene_node_1[j]) - 1]

                    if (which_chromosome < 0 and self.touched_chromosome[i] == 0 and
                        self.touched_chromosome[j] == 0) or (which_chromosome > -1 and (i == which_chromosome or
                                                                                        j == which_chromosome)):
                        if self.gene_nodes_2[node1].adjacency != -1 and self.gene_nodes_2[node3].adjacency != -1:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node1
                            result.node_2 = node3
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.TWO_REVERSED_PLUS_ONE

                            return result
                        elif self.gene_nodes_2[node1].adjacency != -1 and self.gene_nodes_2[node4].adjacency != -1:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node1
                            result.node_2 = node4
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.TWO_PLUS_ONE

                            return result
                        elif self.gene_nodes_2[node2].adjacency != -1 and self.gene_nodes_2[node3].adjacency != -1:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node2
                            result.node_2 = node3
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.ONE_PLUS_TWO

                            return result
                        elif self.gene_nodes_2[node2].adjacency != -1 and self.gene_nodes_2[node4].adjacency != -1:
                            result.chromosome_1 = i
                            result.chromosome_2 = j
                            result.node_1 = node2
                            result.node_2 = node4
                            result.node_3 = -1
                            result.node_4 = -1
                            result.operation_type = OperationTypes.FUSION
                            result.operation_subtype = FusionSubtypes.ONE_PLUS_TWO_REVERSED

                            return result

        elif operation_type == OperationTypes.INVERSION:
            if which_chromosome > -1:
                operation: DCJOperation = self.get_reversal(which_chromosome)

                if operation.chromosome_1 != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome)):
                    if self.touched_chromosome[c] == 0:
                        operation: DCJOperation = self.get_reversal(c)

                        if operation.chromosome_1 != -1:
                            return operation
        elif operation_type == OperationTypes.TRANSLOCATION:
            if which_chromosome > -1:
                operation: DCJOperation = self.get_translocation(which_chromosome, which_chromosome)

                if operation.chromosome_1 != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome) - 1):
                    if self.touched_chromosome[c] == 0:
                        operation: DCJOperation = self.get_translocation(c, which_chromosome)

                        if operation.chromosome_1 != -1:
                            return operation

        return result

    def get_reversal(self, chromosome: int) -> DCJOperation:  # Renamed from findReversal() for clarity
        """
        Gets the reversal of a given chromosome

        Parameters
        ----------
        chromosome
            Chromosome to perform reversal on

        Returns
        -------
        [int]
            List of result attributes
        """
        result: DCJOperation = DCJOperation()
        result.chromosome_1 = -1

        for i in range(int(len(self.chromosomes_for_gene_node_1[chromosome]) / 2)):
            node1: int

            if i == 0:
                node1 = -1
            else:
                node1 = self.chromosomes_for_gene_node_1[chromosome][(i * 2) - 1]

            node2: int = self.chromosomes_for_gene_node_1[chromosome][i * 2]

            if node1 != 1 and self.gene_nodes_2[node1].adjacency != node2 or \
               node2 != -1 and self.gene_nodes_2[node2].adjacency != node1:
                node3: int

                if node1 == -1:
                    node3 = -1
                else:
                    node3 = self.gene_nodes_2[node1].adjacency

                node4: int = self.gene_nodes_2[node2].adjacency

                if node3 != -1 and \
                   self.gene_nodes_1[node3].chromosome_index == chromosome and \
                   self.gene_nodes_1[node3].chromosome_position > int(i * 2):
                    node5: int = self.gene_nodes_1[node3].adjacency

                    if (node5 == -1) or (self.gene_nodes_1[node3].chromosome_position <
                                         self.gene_nodes_1[node5].chromosome_position):
                        result.chromosome_1 = chromosome
                        result.chromosome_2 = chromosome
                        result.node_1 = node1
                        result.node_2 = node3
                        result.node_3 = node2
                        result.node_4 = node5
                        result.operation_type = OperationTypes.INVERSION
                        result.operation_subtype = 1

                        return result
                elif node4 != -1 and \
                        self.gene_nodes_1[node4].chromosome_index == chromosome and \
                        self.gene_nodes_1[node4].chromosome_position > int(i * 2):
                    node6: int = self.gene_nodes_1[node4].adjacency

                    if node6 != -1 and \
                       self.gene_nodes_1[node6].chromosome_position < self.gene_nodes_1[node4].chromosome_position:
                        result.chromosome_1 = chromosome
                        result.chromosome_2 = chromosome
                        result.node_1 = node1
                        result.node_2 = node6
                        result.node_3 = node2
                        result.node_4 = node4
                        result.operation_type = OperationTypes.INVERSION
                        result.operation_subtype = 1

                        return result

        return result

    # Renamed from findTranslocation() for clarity
    def get_translocation(self, chromosome: int, which_chromosome: int) -> DCJOperation:
        """
        Gets the translocation of a given chromosome

        Parameters
        ----------
        chromosome
            Chromosome to perform translocation on
        which_chromosome
            Which chromosome is being operated on

        Returns
        -------
        [int]
            List of result attributes
        """
        result: DCJOperation = DCJOperation()
        result.chromosome_1 = -1

        for i in range((len(self.chromosomes_for_gene_node_1[chromosome]) // 2) + 1):
            node1: int = -1

            if i != 0:
                node1 = self.chromosomes_for_gene_node_1[chromosome][(2 * i) - 1]

            node2: int = -1

            if i != len(self.chromosomes_for_gene_node_1[chromosome]) // 2:
                node2 = self.chromosomes_for_gene_node_1[chromosome][2 * i]

            if node1 != -1 and self.gene_nodes_2[node1].adjacency != node2 or \
               node2 != -1 and self.gene_nodes_2[node2].adjacency != node1:
                node3: int = -1

                if node1 != -1:
                    node3 = self.gene_nodes_2[node1].adjacency

                node4: int = -1

                if node2 != -1:
                    node4 = self.gene_nodes_2[node2].adjacency

                if (node3 != -1 and self.gene_nodes_1[node3].chromosome_index != chromosome) and \
                   (which_chromosome > -1 or self.touched_chromosome[self.gene_nodes_1[node3].chromosome_index] == 0):
                    node5: int = self.gene_nodes_1[node3].adjacency

                    if node2 + node5 != -2:
                        sub_type: int = TranslocationSubtypes.AnCnBD

                        if (self.chromosomes_for_gene_node_1[self.gene_nodes_1[node3].chromosome_index][0] == node3) or\
                           (node5 != -1 and self.gene_nodes_1[node5].chromosome_position <
                                self.gene_nodes_1[node3].chromosome_position):
                            sub_type = TranslocationSubtypes.ADCB

                        result.chromosome_1 = chromosome
                        result.chromosome_2 = self.gene_nodes_1[node3].chromosome_index
                        result.node_1 = node1
                        result.node_2 = node3
                        result.node_3 = node2
                        result.node_4 = node5
                        result.operation_type = OperationTypes.TRANSLOCATION
                        result.operation_subtype = sub_type

                        return result

                if (node4 != -1 and self.gene_nodes_1[node4].chromosome_index != chromosome) and \
                   (which_chromosome > -1 or self.touched_chromosome[self.gene_nodes_1[node4].chromosome_index] == 0):
                    node6: int = self.gene_nodes_1[node4].adjacency

                    if node1 + node6 != -2:
                        sub_type: int = TranslocationSubtypes.ADCB

                        if (self.chromosomes_for_gene_node_1[self.gene_nodes_1[node4].chromosome_index][0] == node4) or\
                           (node6 != -1 and self.gene_nodes_1[node6].chromosome_position <
                                self.gene_nodes_1[node4].chromosome_position):
                            sub_type = TranslocationSubtypes.AnCnBD

                        result.chromosome_1 = chromosome
                        result.chromosome_2 = self.gene_nodes_1[node4].chromosome_index
                        result.node_1 = node1
                        result.node_2 = node6
                        result.node_3 = node2
                        result.node_4 = node4
                        result.operation_type = OperationTypes.TRANSLOCATION
                        result.operation_subtype = sub_type

                        return result

        return result