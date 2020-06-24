from copy import deepcopy
from enum import IntEnum
from typing import List, Optional

from GenomeInString import GenomeInString

"""                               
 Double-Cut-and-Join Operations class                        
                                                             
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


class GeneNodeAttributes(IntEnum):
    ADJACENCY: int = 0
    CHROMOSOMES_INDEX: int = 1
    CHROMOSOME_POSITION: int = 2
    LINKED_ADJACENCY: int = 3


class OperationTypes(IntEnum):
    INVERSION: int = 1
    TRANSLOCATION: int = 2
    FISSION: int = 3
    FUSION: int = 4


class TranslocationSubtypes(IntEnum):
    ADCB: int = 1
    AnCnBD: int = 2


class OperationOptions(IntEnum):
    CHROMOSOME1: int = 0
    CHROMOSOME2: int = 1
    NODE1: int = 2
    NODE2: int = 3
    NODE3: int = 4
    NODE4: int = 5
    TYPE_OF_OPERATION: int = 6
    OPERATION_SUBTYPE: int = 7


class DCJOperations:
    """
    Attributes
    ----------
    genome1 : GenomeInString
        Genome 1 as a list of strings
    genome2 : GenomeInString
        Genome 2 as a list of strings
    gene_number : int
        How many genes are in both genomes (must be the same)
    max_chromosome_number : int
        Highest number of chromosomes in either genome
    min_chromosome_number : int
        Lower number of chromosomes in either genome
    gene_node1 : List[List[int]]
        First gene node for performing DCJ operations
    gene_node2 : List[List[int]]
        Second gene node for performing DCJ operations
    chromosome_in_gene_node1 : List[List[int]]
        Chromosome for first gene node
    chromosome_in_gene_node2 : List[List[int]]
        Chromosome for second gene node
    gene_node_in_string : List[str]
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
        self.gene_number: int = int()

        # Parameters for choosing an appropriate operation
        self.max_chromosome_number: int = int()
        self.min_chromosome_number: int = int()

        # For each row, use GeneNodeAttributes to index
        self.gene_node1: List[Optional[List[int]]] = list(list())
        self.gene_node2: List[Optional[List[int]]] = list(list())

        self.chromosome_in_gene_node1: List[Optional[List[int]]] = list(list())
        self.chromosome_in_gene_node2: List[Optional[List[int]]] = list(list())
        self.gene_node_in_string: List[List[str]] = list(list(str()))
        self.touched_chromosome: List[int] = list()

    def initial_value(self):
        """
        Initializer
        """
        gene_number1: int = 0
        gene_number2: int = 0

        # Split genome1 into genes and save the number of them
        for i in range(len(self.genome1.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome1.chromosomes[i])
            gene_number1 += len(genes)

        # Split genome2 into genes and save the number of them
        for i in range(len(self.genome2.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome2.chromosomes[i])
            gene_number2 += len(genes)

        self.gene_number = gene_number1
        self.gene_node_in_string = [[str(), str()] for _ in range((self.gene_number * 2))]
        self.gene_node1 = [[int(), int(), int()] for _ in range((self.gene_number * 2))]
        self.gene_node2 = [[int(), int(), int()] for _ in range((self.gene_number * 2))]
        self.chromosome_in_gene_node1 = [[int()] for _ in range(len(self.genome1.chromosomes))]
        self.chromosome_in_gene_node2 = [[int()] for _ in range(len(self.genome2.chromosomes))]

        index: int = 0

        # Split each chromosome of genome1 into genes and save genes into gene_node1 and chromosome_in_gene_node1
        for i in range(len(self.genome1.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome1.chromosomes[i])
            self.chromosome_in_gene_node1[i] = [int() for _ in range((len(genes) * 2))]

            for j in range(len(genes)):
                sign: str = "+"
                genome_name: str = genes[j]

                if genome_name[0] == "-":
                    sign = "-"
                    genome_name = genes[j][1:]

                if sign == "+":
                    self.gene_node_in_string[index][0] = genome_name
                    self.gene_node_in_string[index][1] = "t"
                    index += 1
                    self.gene_node_in_string[index][0] = genome_name
                    self.gene_node_in_string[index][1] = "h"
                else:
                    self.gene_node_in_string[index][0] = genome_name
                    self.gene_node_in_string[index][1] = "h"
                    index += 1
                    self.gene_node_in_string[index][0] = genome_name
                    self.gene_node_in_string[index][1] = "t"

                index += 1

                if j == 0:
                    self.gene_node1[index - 2][GeneNodeAttributes.ADJACENCY] = -1
                else:
                    self.gene_node1[index - 2][GeneNodeAttributes.ADJACENCY] = index - 3

                if j == len(genes) - 1:
                    self.gene_node1[index - 1][GeneNodeAttributes.ADJACENCY] = -1
                else:
                    self.gene_node1[index - 1][GeneNodeAttributes.ADJACENCY] = index

                self.gene_node1[index - 2][GeneNodeAttributes.CHROMOSOMES_INDEX] = i
                self.gene_node1[index - 1][GeneNodeAttributes.CHROMOSOMES_INDEX] = i
                self.gene_node1[index - 2][GeneNodeAttributes.CHROMOSOME_POSITION] = j * 2
                self.gene_node1[index - 1][GeneNodeAttributes.CHROMOSOME_POSITION] = (j * 2) + 1
                self.chromosome_in_gene_node1[i][j * 2] = index - 2
                self.chromosome_in_gene_node1[i][(j * 2) + 1] = index - 1

        # Split each chromosome of genome2 into genes and save genes into gene_node2 and chromosome_in_gene_node2
        for i in range(len(self.genome2.chromosomes)):
            genes: List[str] = split_at_whitespace(self.genome2.chromosomes[i])
            pre_node_index: int = -1
            self.chromosome_in_gene_node2[i] = [int() for _ in range((len(genes) * 2))]

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

                self.gene_node2[genome_node_index1][GeneNodeAttributes.CHROMOSOMES_INDEX] = i
                self.gene_node2[genome_node_index2][GeneNodeAttributes.CHROMOSOMES_INDEX] = i
                self.gene_node2[genome_node_index1][GeneNodeAttributes.CHROMOSOME_POSITION] = j * 2
                self.gene_node2[genome_node_index2][GeneNodeAttributes.CHROMOSOME_POSITION] = (j * 2) + 1
                self.chromosome_in_gene_node2[i][j * 2] = genome_node_index1
                self.chromosome_in_gene_node2[i][(j * 2) + 1] = genome_node_index2

                self.gene_node2[genome_node_index1][0] = pre_node_index

                if j != 0:
                    self.gene_node2[pre_node_index][0] = genome_node_index1

                pre_node_index = genome_node_index2

            self.gene_node2[pre_node_index][0] = -1

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
        for i in range(len(self.gene_node_in_string)):
            if g_string == self.gene_node_in_string[i][0] + self.gene_node_in_string[i][1]:
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
                genome1: str = self.gene_node_in_string[node_list[i][j * 2]][0]
                tail_or_head1: str = self.gene_node_in_string[node_list[i][j * 2]][1]
                genome2: str = self.gene_node_in_string[node_list[i][(j * 2) + 1]][0]
                tail_or_head2: str = self.gene_node_in_string[node_list[i][(j * 2) + 1]][1]

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

        self.touched_chromosome = [int() for _ in range(len(self.chromosome_in_gene_node1))]
        self.min_chromosome_number = min_chromosome
        self.max_chromosome_number = max_chromosome

        more: bool = True

        while current_operation < number_of_operations and more:
            more = False
            operation: List[int] = self.expand_all_operations(types_of_operation, which_chromosome)

            if operation[OperationOptions.CHROMOSOME1] != -1:
                if operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.INVERSION:
                    print("Inversion: " + str(operation[OperationOptions.CHROMOSOME1]) + "\n")
                elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.TRANSLOCATION:
                    print("Translocation: " + str(operation[OperationOptions.CHROMOSOME1]) + " " + str(
                        operation[OperationOptions.CHROMOSOME2]) + "\n")
                elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.FISSION:
                    print("Fission: " + str(operation[OperationOptions.CHROMOSOME1]) + " " + str(
                        operation[OperationOptions.CHROMOSOME2]) + "\n")
                elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.FUSION:
                    print("Fusion: " + str(operation[OperationOptions.CHROMOSOME1]) + " " + str(
                        operation[OperationOptions.CHROMOSOME2]) + "\n")

                self.update_all_values(operation, which_chromosome)

                all_steps[step_index] = self.get_genome_in_string(self.chromosome_in_gene_node1)
                step_index += 1
                current_operation += 1
                more = True

        return deepcopy(all_steps)[:step_index]

    def update_all_values(self, operation: List[int], which_chromosome: int):
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
            if operation[OperationOptions.CHROMOSOME1] != -1:
                self.touched_chromosome[operation[OperationOptions.CHROMOSOME1]] = 1

            if operation[OperationOptions.CHROMOSOME2] != -1:
                self.touched_chromosome[operation[OperationOptions.CHROMOSOME2]] = 1

        if operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.INVERSION:
            start_node: int = self.gene_node1[operation[OperationOptions.NODE3]][GeneNodeAttributes.CHROMOSOME_POSITION]
            end_node: int = self.gene_node1[operation[OperationOptions.NODE2]][GeneNodeAttributes.CHROMOSOME_POSITION]
            chromosome_here: int = operation[OperationOptions.CHROMOSOME1]
            new_chromosome_here: List[int] = [int() for _ in range(len(self.chromosome_in_gene_node1[chromosome_here]))]

            for i in range(start_node):
                new_chromosome_here[i] = self.chromosome_in_gene_node1[chromosome_here][i]

            for i in range(start_node, end_node + 1):
                new_chromosome_here[i] = self.chromosome_in_gene_node1[chromosome_here][end_node - i + start_node]

            for i in range(end_node + 1, len(new_chromosome_here)):
                new_chromosome_here[i] = self.chromosome_in_gene_node1[chromosome_here][i]

            for i in range(start_node, end_node + 1):
                self.gene_node1[new_chromosome_here[i]][GeneNodeAttributes.CHROMOSOME_POSITION] = i

            self.chromosome_in_gene_node1[chromosome_here] = new_chromosome_here
        elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.TRANSLOCATION:
            length1: int = len(self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]])
            length2: int = len(self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]])
            lengtha: int = 0

            if operation[OperationOptions.NODE1] != -1:
                lengtha = self.gene_node1[operation[OperationOptions.NODE1]][GeneNodeAttributes.CHROMOSOME_POSITION] + 1

            lengthb: int = length1 - lengtha
            new_chromosome1: List[int] = list()
            new_chromosome2: List[int] = list()

            if operation[OperationOptions.OPERATION_SUBTYPE] == TranslocationSubtypes.ADCB:
                lengthc: int = 0

                if operation[OperationOptions.NODE4] != -1:
                    lengthc = self.gene_node1[operation[OperationOptions.NODE4]][
                                  GeneNodeAttributes.CHROMOSOME_POSITION] + 1

                lengthd: int = length2 - lengthc
                new_chromosome1 = [int() for _ in range(lengtha + lengthd)]
                new_chromosome2 = [int() for _ in range(lengthc + lengthb)]

                for i in range(lengtha):
                    new_chromosome1[i] = self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]][i]

                for i in range(lengthd):
                    new_chromosome1[i + lengtha] = \
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]][i + lengthc]

                for i in range(lengthc):
                    new_chromosome2[i] = self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]][i]

                for i in range(lengthb):
                    new_chromosome2[i + lengthc] = \
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]][i + lengtha]
            elif operation[OperationOptions.OPERATION_SUBTYPE] == TranslocationSubtypes.AnCnBD:
                lengthc: int = 0

                if operation[OperationOptions.NODE2] != -1:
                    lengthc = self.gene_node1[operation[OperationOptions.NODE2]][
                                  GeneNodeAttributes.CHROMOSOME_POSITION] + 1

                lengthd: int = length2 - lengthc
                new_chromosome1 = [int() for _ in range(lengtha + lengthc)]
                new_chromosome2 = [int() for _ in range(lengthb + lengthd)]

                for i in range(lengtha):
                    new_chromosome1[i] = deepcopy(
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]][i])

                for i in range(lengthc):
                    new_chromosome1[i + lengtha] = deepcopy(
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]][lengthc - i - 1])

                for i in range(lengthb):
                    new_chromosome2[i] = deepcopy(
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]][length1 - i - 1])

                for i in range(lengthd):
                    new_chromosome2[i + lengthb] = deepcopy(
                        self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]][i + lengthc])

            for i in range(len(new_chromosome1)):
                self.gene_node1[new_chromosome1[i]][GeneNodeAttributes.CHROMOSOMES_INDEX] = operation[
                    OperationOptions.CHROMOSOME1]
                self.gene_node1[new_chromosome1[i]][GeneNodeAttributes.CHROMOSOME_POSITION] = i

            for i in range(len(new_chromosome2)):
                self.gene_node1[new_chromosome2[i]][GeneNodeAttributes.CHROMOSOMES_INDEX] = operation[
                    OperationOptions.CHROMOSOME2]
                self.gene_node1[new_chromosome2[i]][GeneNodeAttributes.CHROMOSOME_POSITION] = i

            self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME1]] = new_chromosome1
            self.chromosome_in_gene_node1[operation[OperationOptions.CHROMOSOME2]] = new_chromosome2
        elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.FISSION:
            chromosome: int = operation[OperationOptions.CHROMOSOME1]
            length: int = len(self.chromosome_in_gene_node1[chromosome])
            length1: int = self.gene_node1[operation[OperationOptions.NODE1]][
                               GeneNodeAttributes.CHROMOSOME_POSITION] + 1
            length2: int = length - length1
            new_chromosome1: List[int] = [int() for _ in range(length1)]
            new_chromosome2: List[int] = [int() for _ in range(length2)]

            if length1 >= 0:
                new_chromosome1[:length1] = deepcopy(self.chromosome_in_gene_node1[
                                                         operation[OperationOptions.CHROMOSOME1]])

            if length2 >= 0:
                new_chromosome2[:length2] = deepcopy(self.chromosome_in_gene_node1[
                                                         operation[OperationOptions.CHROMOSOME1]])[length1:]

            chromosome_in_gene_node1_temp: List[List[int]] = [
                [int()] for _ in range((len(self.chromosome_in_gene_node1) + 1))]

            chromosome_in_gene_node1_temp[:len(self.chromosome_in_gene_node1)] = deepcopy(self.chromosome_in_gene_node1)

            chromosome_in_gene_node1_temp[operation[OperationOptions.CHROMOSOME1]] = new_chromosome1
            chromosome_in_gene_node1_temp[len(self.chromosome_in_gene_node1)] = new_chromosome2

            self.chromosome_in_gene_node1 = [[int()] for _ in range(len(chromosome_in_gene_node1_temp))]

            self.chromosome_in_gene_node1[:len(chromosome_in_gene_node1_temp)] = deepcopy(chromosome_in_gene_node1_temp)
        elif operation[OperationOptions.TYPE_OF_OPERATION] == OperationTypes.FUSION:
            chromosome1: int = operation[OperationOptions.CHROMOSOME1]
            chromosome2: int = operation[OperationOptions.CHROMOSOME2]
            length1: int = len(self.chromosome_in_gene_node1[chromosome1])
            length2: int = len(self.chromosome_in_gene_node1[chromosome2])
            length: int = length1 + length2
            new_chromosome: List[int] = [int() for _ in range(length)]

            if operation[OperationOptions.OPERATION_SUBTYPE] == 1:
                for i in range(length2):
                    new_chromosome[i] = self.chromosome_in_gene_node1[chromosome2][length2 - i]

                new_chromosome[length2:length1] = deepcopy(self.chromosome_in_gene_node1[chromosome1])
            elif operation[OperationOptions.OPERATION_SUBTYPE] == 2:
                new_chromosome[:length2] = deepcopy(self.chromosome_in_gene_node1[chromosome2])
                new_chromosome[length2:length1] = deepcopy(self.chromosome_in_gene_node1[chromosome1])
            elif operation[OperationOptions.OPERATION_SUBTYPE] == 3:
                new_chromosome[:length1] = deepcopy(self.chromosome_in_gene_node1[chromosome1])
                new_chromosome[length1:length2] = deepcopy(self.chromosome_in_gene_node1[chromosome2])
            elif operation[OperationOptions.OPERATION_SUBTYPE] == 4:
                new_chromosome[:length1] = deepcopy(self.chromosome_in_gene_node1[chromosome1])

                for i in range(length2):
                    new_chromosome[i + length1] = self.chromosome_in_gene_node1[chromosome2][length2 - i]

            chromosome_small: int = chromosome1

            if chromosome2 < chromosome1:
                chromosome_small = chromosome2

            chromosome_big: int = chromosome2

            if chromosome1 > chromosome2:
                chromosome_big = chromosome1

            self.chromosome_in_gene_node1[chromosome_small] = new_chromosome
            self.chromosome_in_gene_node1[chromosome_big] = [0]

            for i in range(len(new_chromosome)):
                self.gene_node1[new_chromosome[i]][GeneNodeAttributes.CHROMOSOME_POSITION] = i

            for i in range(len(self.gene_node1)):
                if self.gene_node1[i][1] >= chromosome1 and \
                        self.gene_node1[i][GeneNodeAttributes.CHROMOSOMES_INDEX] >= chromosome2:
                    self.gene_node1[i][1] = self.gene_node1[i][GeneNodeAttributes.CHROMOSOMES_INDEX] - 1

            chromosome_in_gene_node1_temp: List[List[int]] = [
                [int()] for _ in range((len(self.chromosome_in_gene_node1) - 1))]
            index: int = 0

            for ints in self.chromosome_in_gene_node1:
                if ints != [0]:
                    chromosome_in_gene_node1_temp[index] = ints
                    index += 1

            self.chromosome_in_gene_node1 = [[int()] for _ in range(len(chromosome_in_gene_node1_temp))]
            self.chromosome_in_gene_node1 = chromosome_in_gene_node1_temp

        if operation[OperationOptions.NODE1] != -1:
            self.gene_node1[operation[OperationOptions.NODE1]][GeneNodeAttributes.ADJACENCY] = operation[
                OperationOptions.NODE2]

        if operation[OperationOptions.NODE2] != -1:
            self.gene_node1[operation[OperationOptions.NODE2]][GeneNodeAttributes.ADJACENCY] = operation[
                OperationOptions.NODE1]

        if operation[OperationOptions.NODE3] != -1:
            self.gene_node1[operation[OperationOptions.NODE3]][GeneNodeAttributes.ADJACENCY] = operation[
                OperationOptions.NODE4]

        if operation[OperationOptions.NODE4] != -1:
            self.gene_node1[operation[OperationOptions.NODE4]][GeneNodeAttributes.ADJACENCY] = operation[
                OperationOptions.NODE3]

    # Renamed from findAOperation() for accuracy
    def expand_all_operations(self, types_of_operation: List[int], which_chromosome: int) -> List[int]:
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
        result: List[int] = [int() for _ in range(8)]
        result[OperationOptions.CHROMOSOME1] = -1

        for operation in types_of_operation:
            result = self.get_operation_attributes(which_chromosome, operation)

            if result[OperationOptions.CHROMOSOME1] != -1:
                return result

        return result

    # Renamed from findATypeOperation() for accuracy
    def get_operation_attributes(self, which_chromosome: int, operation_type: int) -> List[int]:
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
        result: List[int] = [int() for _ in range(8)]
        result[OperationOptions.CHROMOSOME1] = -1

        if (operation_type == OperationTypes.FUSION and
            len(self.chromosome_in_gene_node1) <= self.min_chromosome_number) or \
                (operation_type == OperationTypes.FISSION and
                 len(self.chromosome_in_gene_node1) >= self.max_chromosome_number):
            return result

        if operation_type == OperationTypes.FISSION:
            for i in range(len(self.chromosome_in_gene_node2) - 1):
                node1: int = self.chromosome_in_gene_node2[i][0]
                node2: int = self.chromosome_in_gene_node2[i][len(self.chromosome_in_gene_node2[i]) - 1]

                for j in range(i + 1, len(self.chromosome_in_gene_node2)):
                    node3: int = self.chromosome_in_gene_node2[j][0]
                    node4: int = self.chromosome_in_gene_node2[j][len(self.chromosome_in_gene_node2[j]) - 1]

                    if self.gene_node1[node1][GeneNodeAttributes.ADJACENCY] == node3 or \
                            self.gene_node1[node1][GeneNodeAttributes.ADJACENCY] == node4 or \
                            self.gene_node1[node2][GeneNodeAttributes.ADJACENCY] == node3 or \
                            self.gene_node1[node2][GeneNodeAttributes.ADJACENCY] == node4:
                        if (-1 < which_chromosome == self.gene_node1[node1][GeneNodeAttributes.CHROMOSOMES_INDEX]) \
                                or (which_chromosome < 0 and
                                    self.touched_chromosome[
                                        self.gene_node1[node1][GeneNodeAttributes.CHROMOSOMES_INDEX]] == 0):
                            result[OperationOptions.CHROMOSOME1] = self.gene_node1[node1][
                                GeneNodeAttributes.CHROMOSOMES_INDEX]
                            result[OperationOptions.NODE1] = node1
                        elif (-1 < which_chromosome == self.gene_node1[node2][GeneNodeAttributes.CHROMOSOMES_INDEX]) \
                                or (which_chromosome < 0 and self.touched_chromosome[self.gene_node1[node2]
                                    [GeneNodeAttributes.CHROMOSOMES_INDEX]] == 0):
                            result[OperationOptions.CHROMOSOME1] = self.gene_node1[node2][
                                GeneNodeAttributes.CHROMOSOMES_INDEX]
                            result[OperationOptions.NODE1] = node2

                        if self.gene_node1[node1][GeneNodeAttributes.ADJACENCY] == node3 or \
                                self.gene_node1[node2][GeneNodeAttributes.ADJACENCY] == node3:
                            result[OperationOptions.NODE3] = node3
                        elif self.gene_node1[node1][GeneNodeAttributes.ADJACENCY] == node4 or \
                                self.gene_node1[node2][GeneNodeAttributes.ADJACENCY] == node4:
                            result[OperationOptions.NODE3] = node4

                        result[OperationOptions.CHROMOSOME2] = len(self.chromosome_in_gene_node1)
                        result[OperationOptions.NODE2] = -1
                        result[OperationOptions.NODE4] = -1
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.FISSION
                        result[OperationOptions.OPERATION_SUBTYPE] = 1

                        return result

            for i in range(len(self.chromosome_in_gene_node2)):
                node1 = self.chromosome_in_gene_node2[i][0]
                node2 = self.chromosome_in_gene_node2[i][len(self.chromosome_in_gene_node2[i]) - 1]

                chromosome_node1: int = self.gene_node1[node1][GeneNodeAttributes.CHROMOSOMES_INDEX]
                node3: int = self.gene_node1[node1][GeneNodeAttributes.ADJACENCY]

                if node3 != -1 and (chromosome_node1 == which_chromosome or
                                    which_chromosome < 0 and self.touched_chromosome[chromosome_node1] == 0):
                    result[OperationOptions.CHROMOSOME1] = self.gene_node1[node1][GeneNodeAttributes.CHROMOSOMES_INDEX]
                    result[OperationOptions.CHROMOSOME2] = len(self.chromosome_in_gene_node1)
                    result[OperationOptions.NODE1] = node1
                    result[OperationOptions.NODE2] = -1
                    result[OperationOptions.NODE3] = node3
                    result[OperationOptions.NODE4] = -1
                    result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.FISSION
                    result[OperationOptions.OPERATION_SUBTYPE] = 1

                    return result

                chromosome_node2: int = self.gene_node1[node2][GeneNodeAttributes.CHROMOSOMES_INDEX]
                node4: int = self.gene_node1[node2][GeneNodeAttributes.ADJACENCY]

                if node4 != -1 and (chromosome_node2 == which_chromosome or
                                    which_chromosome < 0 and self.touched_chromosome[chromosome_node2] == 0):
                    result[OperationOptions.CHROMOSOME1] = self.gene_node1[node2][GeneNodeAttributes.CHROMOSOMES_INDEX]
                    result[OperationOptions.CHROMOSOME2] = len(self.chromosome_in_gene_node1)
                    result[OperationOptions.NODE1] = node2
                    result[OperationOptions.NODE2] = -1
                    result[OperationOptions.NODE3] = node4
                    result[OperationOptions.NODE4] = -1
                    result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.FISSION
                    result[OperationOptions.OPERATION_SUBTYPE] = 1

                    return result
        elif operation_type == OperationTypes.FUSION:
            for i in range(len(self.chromosome_in_gene_node1) - 1):
                node1: int = self.chromosome_in_gene_node1[i][0]
                node2: int = self.chromosome_in_gene_node1[i][len(self.chromosome_in_gene_node1[i]) - 1]

                for j in range(i + 1, len(self.chromosome_in_gene_node1)):
                    node3: int = self.chromosome_in_gene_node1[j][0]
                    node4: int = self.chromosome_in_gene_node1[j][len(self.chromosome_in_gene_node1[j]) - 1]

                    if (which_chromosome < 0 and self.touched_chromosome[i] == 0 and
                        self.touched_chromosome[j] == 0) or \
                            (which_chromosome > -1 and (i == which_chromosome or j == which_chromosome)):
                        if self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] == node3:
                            result[OperationOptions.NODE1] = node1
                            result[OperationOptions.NODE2] = node3
                            result[OperationOptions.OPERATION_SUBTYPE] = 1
                        elif self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] == node4:
                            result[OperationOptions.NODE1] = node1
                            result[OperationOptions.NODE2] = node4
                            result[OperationOptions.OPERATION_SUBTYPE] = 2
                        elif self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] == node3:
                            result[OperationOptions.NODE1] = node2
                            result[OperationOptions.NODE2] = node3
                            result[OperationOptions.OPERATION_SUBTYPE] = 3
                        elif self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] == node4:
                            result[OperationOptions.NODE1] = node2
                            result[OperationOptions.NODE2] = node4
                            result[OperationOptions.OPERATION_SUBTYPE] = 4

                        result[OperationOptions.CHROMOSOME1] = i
                        result[OperationOptions.CHROMOSOME2] = j
                        result[OperationOptions.NODE3] = -1
                        result[OperationOptions.NODE4] = -1
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.FUSION

                        return result

            for i in range(len(self.chromosome_in_gene_node1) - 1):
                node1: int = self.chromosome_in_gene_node1[i][0]
                node2: int = self.chromosome_in_gene_node1[i][len(self.chromosome_in_gene_node1[i]) - 1]

                for j in range(i + 1, len(self.chromosome_in_gene_node1)):
                    node3: int = self.chromosome_in_gene_node1[j][0]
                    node4: int = self.chromosome_in_gene_node1[j][len(self.chromosome_in_gene_node1[j]) - 1]

                    if (which_chromosome < 0 and self.touched_chromosome[i] == 0 and
                        self.touched_chromosome[j] == 0) or \
                            (which_chromosome > -1 and (i == which_chromosome or j == which_chromosome)):
                        if self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] != -1 and \
                                self.gene_node2[node3][GeneNodeAttributes.ADJACENCY] != -1:
                            result[OperationOptions.NODE1] = node1
                            result[OperationOptions.NODE2] = node3
                            result[OperationOptions.OPERATION_SUBTYPE] = 1
                        elif self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] != -1 and \
                                self.gene_node2[node4][GeneNodeAttributes.ADJACENCY] != -1:
                            result[OperationOptions.NODE1] = node1
                            result[OperationOptions.NODE2] = node4
                            result[OperationOptions.OPERATION_SUBTYPE] = 2
                        elif self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] != -1 and \
                                self.gene_node2[node3][GeneNodeAttributes.ADJACENCY] != -1:
                            result[OperationOptions.NODE1] = node2
                            result[OperationOptions.NODE2] = node3
                            result[OperationOptions.OPERATION_SUBTYPE] = 3
                        elif self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] != -1 and \
                                self.gene_node2[node4][GeneNodeAttributes.ADJACENCY] != -1:
                            result[OperationOptions.NODE1] = node2
                            result[OperationOptions.NODE2] = node4
                            result[OperationOptions.OPERATION_SUBTYPE] = 4

                        result[OperationOptions.CHROMOSOME1] = i
                        result[OperationOptions.CHROMOSOME2] = j
                        result[OperationOptions.NODE3] = -1
                        result[OperationOptions.NODE4] = -1
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.FUSION

                        return result

        elif operation_type == OperationTypes.INVERSION:
            if which_chromosome > -1:
                operation: List[int] = self.get_reversal(which_chromosome)

                if operation[OperationOptions.CHROMOSOME1] != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome)):
                    if self.touched_chromosome[c] == 0:
                        operation: List[int] = self.get_reversal(c)

                        if operation[OperationOptions.CHROMOSOME1] != -1:
                            return operation
        elif operation_type == OperationTypes.TRANSLOCATION:
            if which_chromosome > -1:
                operation: List[int] = self.get_translocation(which_chromosome, which_chromosome)

                if operation[OperationOptions.CHROMOSOME1] != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome) - 1):
                    if self.touched_chromosome[c] == 0:
                        operation: List[int] = self.get_translocation(c, which_chromosome)

                        if operation[OperationOptions.CHROMOSOME1] != -1:
                            return operation

        return result

    def get_reversal(self, chromosome: int) -> List[int]:  # Renamed from findReversal() for clarity
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
        result: List[int] = [int() for _ in range(8)]
        result[OperationOptions.CHROMOSOME1] = -1

        for i in range(len(self.chromosome_in_gene_node1[chromosome]) // 2):
            node1: int = -1

            if i != 0:
                node1 = self.chromosome_in_gene_node1[chromosome][(2 * i) - 1]

            node2: int = self.chromosome_in_gene_node1[chromosome][2 * i]

            if node1 != -1 and self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] != node2 or \
                    node2 != -1 and self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] != node1:
                node3: int = -1

                if node1 != -1:
                    node3 = self.gene_node2[node1][GeneNodeAttributes.ADJACENCY]

                node4: int = self.gene_node2[node2][GeneNodeAttributes.ADJACENCY]

                if (node3 != -1 and self.gene_node1[node3][GeneNodeAttributes.CHROMOSOMES_INDEX] == chromosome) and \
                        (self.gene_node1[node3][GeneNodeAttributes.CHROMOSOME_POSITION] > (2 * i)):
                    node5: int = self.gene_node1[node3][GeneNodeAttributes.ADJACENCY]

                    if node5 == -1 or self.gene_node1[node3][GeneNodeAttributes.CHROMOSOME_POSITION] < \
                            self.gene_node1[node5][GeneNodeAttributes.CHROMOSOME_POSITION]:
                        result[OperationOptions.CHROMOSOME1] = chromosome
                        result[OperationOptions.CHROMOSOME2] = chromosome
                        result[OperationOptions.NODE1] = node1
                        result[OperationOptions.NODE2] = node3
                        result[OperationOptions.NODE3] = node2
                        result[OperationOptions.NODE4] = node5
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.INVERSION
                        result[OperationOptions.OPERATION_SUBTYPE] = 1

                        return result

                if (node4 != -1 and self.gene_node1[node4][GeneNodeAttributes.CHROMOSOMES_INDEX] == chromosome) and \
                        (self.gene_node1[node4][GeneNodeAttributes.CHROMOSOME_POSITION] > (2 * i)):
                    node6: int = self.gene_node1[node4][GeneNodeAttributes.ADJACENCY]

                    if node6 != -1 or self.gene_node1[node6][GeneNodeAttributes.CHROMOSOME_POSITION] < \
                            self.gene_node1[node4][GeneNodeAttributes.CHROMOSOME_POSITION]:
                        result[OperationOptions.CHROMOSOME1] = chromosome
                        result[OperationOptions.CHROMOSOME2] = chromosome
                        result[OperationOptions.NODE1] = node1
                        result[OperationOptions.NODE2] = node6
                        result[OperationOptions.NODE3] = node2
                        result[OperationOptions.NODE4] = node4
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.INVERSION
                        result[OperationOptions.OPERATION_SUBTYPE] = 1

                        return result

        return result

    # Renamed from findTranslocation() for clarity
    def get_translocation(self, chromosome: int, which_chromosome: int) -> List[int]:
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
        result: List[int] = [int() for _ in range(8)]
        result[OperationOptions.CHROMOSOME1] = -1

        for i in range((len(self.chromosome_in_gene_node1[chromosome]) // 2) + 1):
            node1: int = -1

            if i != 0:
                node1 = self.chromosome_in_gene_node1[chromosome][(2 * i) - 1]

            node2: int = -1

            if i != len(self.chromosome_in_gene_node1[chromosome]) // 2:
                node2 = self.chromosome_in_gene_node1[chromosome][2 * i]

            if node1 != -1 and self.gene_node2[node1][GeneNodeAttributes.ADJACENCY] != node2 or \
                    node2 != -1 and self.gene_node2[node2][GeneNodeAttributes.ADJACENCY] != node1:
                node3: int = -1

                if node1 != -1:
                    node3 = self.gene_node2[node1][GeneNodeAttributes.ADJACENCY]

                node4: int = -1

                if node2 != -1:
                    node4 = self.gene_node2[node2][GeneNodeAttributes.ADJACENCY]

                if (node3 != -1 and self.gene_node1[node3][GeneNodeAttributes.CHROMOSOMES_INDEX] != chromosome) and \
                        (which_chromosome > -1 or self.touched_chromosome[
                            self.gene_node1[node3][GeneNodeAttributes.CHROMOSOMES_INDEX]] == 0):
                    node5: int = self.gene_node1[node3][GeneNodeAttributes.ADJACENCY]

                    if node2 + node5 != -2:
                        sub_type: int = 2

                        if (self.chromosome_in_gene_node1[self.gene_node1[node3][1]][
                                GeneNodeAttributes.ADJACENCY] == node3) or \
                                (node5 != -1 and self.gene_node1[node5][GeneNodeAttributes.CHROMOSOME_POSITION] <
                                 self.gene_node1[node3][GeneNodeAttributes.CHROMOSOME_POSITION]):
                            sub_type = 1

                        result[OperationOptions.CHROMOSOME1] = chromosome
                        result[OperationOptions.CHROMOSOME2] = self.gene_node1[node3][
                            GeneNodeAttributes.CHROMOSOMES_INDEX]
                        result[OperationOptions.NODE1] = node1
                        result[OperationOptions.NODE2] = node3
                        result[OperationOptions.NODE3] = node2
                        result[OperationOptions.NODE4] = node5
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.TRANSLOCATION
                        result[OperationOptions.OPERATION_SUBTYPE] = sub_type

                        return result

                if (node4 != -1 and self.gene_node1[node4][GeneNodeAttributes.CHROMOSOMES_INDEX] != chromosome) and \
                        (which_chromosome > -1 or self.touched_chromosome[
                            self.gene_node1[node4][GeneNodeAttributes.CHROMOSOMES_INDEX]] == 0):
                    node6: int = self.gene_node1[node4][GeneNodeAttributes.ADJACENCY]

                    if node1 + node6 != -2:
                        sub_type: int = 1

                        if (self.chromosome_in_gene_node1[self.gene_node1[node4][1]][
                                GeneNodeAttributes.ADJACENCY] == node4) or \
                                (node6 != -1 and self.gene_node1[node6][GeneNodeAttributes.CHROMOSOME_POSITION] <
                                 self.gene_node1[node4][GeneNodeAttributes.CHROMOSOME_POSITION]):
                            sub_type = 2

                        result[OperationOptions.CHROMOSOME1] = chromosome
                        result[OperationOptions.CHROMOSOME2] = self.gene_node1[node4][
                            GeneNodeAttributes.CHROMOSOMES_INDEX]
                        result[OperationOptions.NODE1] = node1
                        result[OperationOptions.NODE2] = node6
                        result[OperationOptions.NODE3] = node2
                        result[OperationOptions.NODE4] = node4
                        result[OperationOptions.TYPE_OF_OPERATION] = OperationTypes.TRANSLOCATION
                        result[OperationOptions.OPERATION_SUBTYPE] = sub_type

                        return result

        return result
