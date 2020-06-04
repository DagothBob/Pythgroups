from enum import Enum

from GenomeInString import GenomeInString


def insert_character(s: str, i: int, c: str) -> str:
    return s[:i] + c + s[(i + 1):]


class OperationsEnum(Enum):
    INVERSION = 1
    TRANSLOCATION = 2
    FISSION = 3
    FUSION = 4


class DCJOperations:
    def __init__(self, genome1: [str], genome2: [str]):
        self.genome1: GenomeInString = GenomeInString(genome1)
        self.genome2: GenomeInString = GenomeInString(genome2)
        self.gene_number: int = 0
        self.max_chromosome_number: int = 0
        self.min_chromosome_number: int = 0
        self.gene_node1: [[int]] = [[0]]
        self.gene_node2: [[int]] = [[0]]
        self.chromosome_in_gene_node1: [[int]] = [[0]]
        self.chromosome_in_gene_node2: [[int]] = [[0]]
        self.gene_node_in_string: [str] = [" "]
        self.touched_chromosome: [int] = [0]

    def initial_value(self):
        gene_number1: int = 0
        gene_number2: int = 0

        for i in range(len(self.genome1.chromosomes)):
            genes: [str] = self.genome1.chromosomes.split(' ')
            gene_number1 += len(genes)

        for i in range(len(self.genome2.chromosomes)):
            genes: [str] = self.genome2.chromosomes.split(' ')
            gene_number2 += len(genes)

        self.gene_number = gene_number1
        self.gene_node_in_string = [" " * (self.gene_number * 2)] * 2
        self.gene_node1 = [[0] * (self.gene_number * 2)] * 3
        self.gene_node2 = [[0] * (self.gene_number * 2)] * 3
        self.chromosome_in_gene_node1 = [[0] * len(self.genome1.chromosomes)]
        self.chromosome_in_gene_node2 = [[0] * len(self.genome2.chromosomes)]

        index: int = 0

        for i in range(len(self.genome1.chromosomes)):
            genes: [str] = self.genome1.chromosomes[i].split(' ')
            self.chromosome_in_gene_node1[i] = [0 * (len(genes) * 2)]

            for j in range(len(genes)):
                sign: str = "+"
                genome_name: str = genes[j]

                if genome_name[0] == "-":
                    sign = "-"
                    genome_name = genes[j][1:]

                if sign == "+":
                    insert_character(self.gene_node_in_string[index], 0, genome_name)
                    insert_character(self.gene_node_in_string[index], 1, "t")
                    index += 1
                    insert_character(self.gene_node_in_string[index], 0, genome_name)
                    insert_character(self.gene_node_in_string[index], 1, "h")
                else:
                    insert_character(self.gene_node_in_string[index], 0, genome_name)
                    insert_character(self.gene_node_in_string[index], 1, "h")
                    index += 1
                    insert_character(self.gene_node_in_string[index], 0, genome_name)
                    insert_character(self.gene_node_in_string[index], 1, "t")

                index += 1

                if j == 0:
                    self.gene_node1[index - 2][0] = -1
                else:
                    self.gene_node1[index - 2][0] = index - 3

                if j == len(genes) - 1:
                    self.gene_node1[index - 1][0] = -1
                else:
                    self.gene_node1[index - 1][0] = index

                self.gene_node1[index - 2][1] = i
                self.gene_node1[index - 1][1] = i
                self.gene_node1[index - 2][2] = j * 2
                self.gene_node1[index - 1][2] = (j * 2) + 1
                self.chromosome_in_gene_node1[i][j * 2] = index - 2
                self.chromosome_in_gene_node1[i][(j * 2) + 1] = index - 1

        for i in range(len(self.genome2.chromosomes)):
            genes: [str] = self.genome2.chromosomes[i].split(' ')
            pre_node_index: int = -1
            self.chromosome_in_gene_node2[i] = [0] * (len(genes) * 2)

            for j in range(len(genes)):
                sign: str = "+"
                genome_name: str = genes[j]

                if genome_name[0] == "-":
                    sign = "-"
                    genome_name = genes[j][1:]

                genome_node_string1: str = insert_character(genome_name, len(genome_name) - 1, "t")
                genome_node_string2: str = insert_character(genome_name, len(genome_name) - 1, "h")

                if sign == "-":
                    genome_node_string1 = insert_character(genome_name, len(genome_name) - 1, "h")
                    genome_node_string2 = insert_character(genome_name, len(genome_name) - 1, "t")

                genome_node_index1: int = self.find_node_index(genome_node_string1)
                genome_node_index2: int = self.find_node_index(genome_node_string2)

                self.gene_node2[genome_node_index1 - 2][1] = i
                self.gene_node2[genome_node_index2 - 1][1] = i
                self.gene_node2[genome_node_index1 - 2][2] = j * 2
                self.gene_node2[genome_node_index2 - 1][2] = (j * 2) + 1
                self.chromosome_in_gene_node2[i][j * 2] = genome_node_index1 - 2
                self.chromosome_in_gene_node2[i][(j * 2) + 1] = genome_node_index2 - 1

                self.gene_node2[genome_node_index1][0] = pre_node_index

                if j != 0:
                    self.gene_node2[pre_node_index][0] = genome_node_index1

                pre_node_index = genome_node_index2

            self.gene_node2[pre_node_index][0] = -1

    def find_node_index(self, g_string: str) -> int:
        for i in range(len(self.gene_node_in_string)):
            if g_string == self.gene_node_in_string[i][0] + self.gene_node_in_string[i][1]:
                return i

        return -1

    def get_genome_in_string(self, node_list: [[int]]) -> GenomeInString:
        chromosome: [str] = [" " * len(node_list)]

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

    def get_result(self, min_chromosome: int, max_chromosome: int, which_chromosome: int, types_of_operation: [int],
                   number_of_operations: int):
        all_steps: [GenomeInString] = [GenomeInString(None)] * number_of_operations
        step_index: int = 0
        current_operation: int = 0

        self.touched_chromosome = [0] * len(self.chromosome_in_gene_node1)
        self.min_chromosome_number = min_chromosome
        self.max_chromosome_number = max_chromosome

        more: bool = True

        while current_operation < number_of_operations and more:
            more = False
            operation: [int] = self.find_an_operation(types_of_operation, which_chromosome)

            if operation[0] != -1:
                if operation[6] == 1:
                    print("Inversion: " + str(operation[0]) + "\n")
                if operation[6] == 2:
                    print("Translocation: " + str(operation[0]) + " " + str(operation[1]) + "\n")
                if operation[6] == 3:
                    print("Fission: " + str(operation[0]) + " " + str(operation[1]) + "\n")
                if operation[6] == 4:
                    print("Fusion: " + str(operation[0]) + " " + str(operation[1]) + "\n")

                self.update_all_values(operation, which_chromosome)

                genome_here: GenomeInString = self.get_genome_in_string(self.chromosome_in_gene_node1)
                all_steps[step_index] = genome_here
                step_index += 1
                current_operation += 1
                more = True

        result: [GenomeInString] = [GenomeInString(None)] * step_index

        if len(result) >= 0:
            result[:len(result)] = all_steps.copy()

        return result

    def update_all_values(self, operation: [int], which_chromosome: int):
        if which_chromosome == -1:
            if operation[0] != -1:
                self.touched_chromosome[operation[0]] = 1

            if operation[1] != -1:
                self.touched_chromosome[operation[1]] = 1

        if operation[6] == 1:
            start_node: int = self.gene_node1[operation[4]][2]
            end_node: int = self.gene_node1[operation[3]][2]
            chromosome_here: int = operation[0]
            new_chromosome_here: [int] = [0] * len(self.chromosome_in_gene_node1[chromosome_here])

            if start_node >= 0:
                self.chromosome_in_gene_node1[chromosome_here] = new_chromosome_here.copy()[:start_node]

            for i in range(start_node, end_node + 1):
                new_chromosome_here[i] = self.chromosome_in_gene_node1[chromosome_here][end_node - i + start_node]

            if len(new_chromosome_here) - end_node + 1 >= 0:
                new_chromosome_here[(end_node + 1):(len(new_chromosome_here) - end_node + 1)] = \
                    self.chromosome_in_gene_node1[chromosome_here].copy()[(end_node + 1):]

            for i in range(start_node, end_node + 1):
                self.gene_node1[new_chromosome_here[i]][2] = i

            self.chromosome_in_gene_node1[chromosome_here] = new_chromosome_here
        elif operation[6] == 2:
            length1: int = len(self.chromosome_in_gene_node1[operation[0]])
            length2: int = len(self.chromosome_in_gene_node1[operation[1]])
            lengtha: int = 0

            if operation[2] != -1:
                lengtha = self.gene_node1[operation[2]][2] + 1

            lengthb: int = length1 - lengtha
            new_chromosome1: [int] = [0]
            new_chromosome2: [int] = [0]

            if operation[7] == 1:
                lengthc: int = 0

                if operation[5] != -1:
                    lengthc = self.gene_node1[operation[5]][2] + 1

                lengthd: int = length2 - lengthc
                new_chromosome1 = [0] * (lengtha + lengthd)
                new_chromosome2 = [0] * (lengthc + lengthb)

                if lengtha >= 0:
                    new_chromosome1[:lengtha] = self.chromosome_in_gene_node1[operation[0]].copy()

                if lengthd >= 0:
                    new_chromosome1[lengtha:lengthd] = self.chromosome_in_gene_node1[operation[1]].copy()[lengthc:]

                if lengthc >= 0:
                    new_chromosome2[:lengthc] = self.chromosome_in_gene_node1[operation[1]].copy()

                if lengthb >= 0:
                    new_chromosome2[lengthc:lengthb] = self.chromosome_in_gene_node1[operation[0]].copy()[lengtha:]
            elif operation[7] == 2:
                lengthc: int = 0

                if operation[3] != -1:
                    lengthc = self.gene_node1[operation[3]][2] + 1

                lengthd: int = length2 - lengthc
                new_chromosome1 = [0] * (lengtha + lengthc)
                new_chromosome2 = [0] * (lengthb + lengthd)

                if lengtha >= 0:
                    new_chromosome1[:lengtha] = self.chromosome_in_gene_node1[operation[0]].copy()

                for i in range(lengthc):
                    new_chromosome1[i + lengtha] = self.chromosome_in_gene_node1[operation[1]][lengthc - i - 1]

                for i in range(lengthb):
                    new_chromosome2[i] = self.chromosome_in_gene_node1[operation[0]][length1 - i - 1]

                if lengthd >= 0:
                    new_chromosome2[lengthb:lengthd] = self.chromosome_in_gene_node1[operation[1]].copy()[lengthc:]

            for i in range(len(new_chromosome1)):
                self.gene_node1[new_chromosome1[i]][1] = operation[0]
                self.gene_node1[new_chromosome1[i]][2] = i

            for i in range(len(new_chromosome2)):
                self.gene_node1[new_chromosome2[i]][1] = operation[1]
                self.gene_node1[new_chromosome2[i]][2] = i

            self.chromosome_in_gene_node1[operation[0]] = new_chromosome1
            self.chromosome_in_gene_node1[operation[1]] = new_chromosome2
        elif operation[6] == 3:
            chromosome: int = operation[0]
            length: int = len(self.chromosome_in_gene_node1[chromosome])
            length1: int = self.gene_node1[operation[2]][2] + 1
            length2: int = length - length1
            new_chromosome1: [int] = [0] * length1
            new_chromosome2: [int] = [0] * length2

            if length1 >= 0:
                new_chromosome1[:length1] = self.chromosome_in_gene_node1[operation[0]].copy()

            if length2 >= 0:
                new_chromosome2[:length2] = self.chromosome_in_gene_node1[operation[0]].copy()[length1:]

            chromosome_in_gene_node1_temp: [[int]] = [[0] * (len(self.chromosome_in_gene_node1) + 1)]

            chromosome_in_gene_node1_temp[:len(self.chromosome_in_gene_node1)] = self.chromosome_in_gene_node1.copy()

            chromosome_in_gene_node1_temp[operation[0]] = new_chromosome1
            chromosome_in_gene_node1_temp[len(self.chromosome_in_gene_node1)] = new_chromosome2

            self.chromosome_in_gene_node1 = [[0] * len(chromosome_in_gene_node1_temp)]

            self.chromosome_in_gene_node1[:len(chromosome_in_gene_node1_temp)] = chromosome_in_gene_node1_temp.copy()
        elif operation[6] == 4:
            chromosome1: int = operation[0]
            chromosome2: int = operation[1]
            length1: int = len(self.chromosome_in_gene_node1[chromosome1])
            length2: int = len(self.chromosome_in_gene_node1[chromosome2])
            length: int = length1 + length2
            new_chromosome: [int] = [0] * length

            if operation[7] == 1:
                for i in range(length2):
                    new_chromosome[i] = self.chromosome_in_gene_node1[chromosome2][length2 - i]

                new_chromosome[length2:length1] = self.chromosome_in_gene_node1[chromosome1].copy()
            elif operation[7] == 2:
                new_chromosome[:length2] = self.chromosome_in_gene_node1[chromosome2].copy()
                new_chromosome[length2:length1] = self.chromosome_in_gene_node1[chromosome1].copy()
            elif operation[7] == 3:
                new_chromosome[:length1] = self.chromosome_in_gene_node1[chromosome1].copy()
                new_chromosome[length1:length2] = self.chromosome_in_gene_node1[chromosome2].copy()
            elif operation[7] == 4:
                new_chromosome[:length1] = self.chromosome_in_gene_node1[chromosome1].copy()

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
                self.gene_node1[new_chromosome[i]][2] = i

            for i in range(len(self.gene_node1)):
                if self.gene_node1[i][1] >= chromosome1 and self.gene_node1[i][1] >= chromosome2:
                    self.gene_node1[i][1] = self.gene_node1[i][1] - 1

            chromosome_in_gene_node1_temp: [[int]] = [[0] * (len(self.chromosome_in_gene_node1) - 1)]
            index: int = 0

            for ints in self.chromosome_in_gene_node1:
                if ints != [0]:
                    chromosome_in_gene_node1_temp[index] = ints
                    index += 1

            self.chromosome_in_gene_node1 = [[0] * len(chromosome_in_gene_node1_temp)]
            self.chromosome_in_gene_node1 = chromosome_in_gene_node1_temp

        if operation[2] != -1:
            self.gene_node1[operation[2]][0] = operation[3]

        if operation[3] != -1:
            self.gene_node1[operation[3]][0] = operation[2]

        if operation[4] != -1:
            self.gene_node1[operation[4]][0] = operation[5]

        if operation[5] != -1:
            self.gene_node1[operation[5]][0] = operation[4]

    def find_an_operation(self, types_of_operation: [int], which_chromosome: int) -> [int]:
        result: [int] = [0] * 8
        result[0] = -1

        for operation in types_of_operation:
            result = self.find_a_type_operation(which_chromosome, operation)

            if result[0] != -1:
                return result

        return result

    def find_a_type_operation(self, which_chromosome: int, operation_type: int) -> [int]:
        result: [int] = [0] * 8
        result[0] = -1

        if (operation_type == OperationsEnum.FUSION and
                len(self.chromosome_in_gene_node1) <= self.min_chromosome_number) or \
           (operation_type == OperationsEnum.FISSION and
                len(self.chromosome_in_gene_node1) >= self.max_chromosome_number):
            return result

        if operation_type == OperationsEnum.FISSION:
            for i in range(len(self.chromosome_in_gene_node2) - 1):
                node1: int = self.chromosome_in_gene_node2[i][0]
                node2: int = self.chromosome_in_gene_node2[i][len(self.chromosome_in_gene_node2[i]) - 1]

                for j in range(i + 1, len(self.chromosome_in_gene_node2)):
                    node3: int = self.chromosome_in_gene_node2[j][0]
                    node4: int = self.chromosome_in_gene_node2[j][len(self.chromosome_in_gene_node2[j]) - 1]

                    if self.gene_node1[node1][0] == node3 or \
                       self.gene_node1[node1][0] == node4 or \
                       self.gene_node1[node2][0] == node3 or \
                       self.gene_node1[node2][0] == node4:
                        if (-1 < which_chromosome == self.gene_node1[node1][1]) or \
                           (which_chromosome < 0 and self.touched_chromosome[self.gene_node1[node1][1]] == 0):
                            result[0] = self.gene_node1[node1][1]
                            result[2] = node1
                        elif (-1 < which_chromosome == self.gene_node1[node2][1]) or \
                             (which_chromosome < 0 and self.touched_chromosome[self.gene_node1[node2][1]] == 0):
                            result[0] = self.gene_node1[node2][1]
                            result[2] = node2

                        if self.gene_node1[node1][0] == node3 or self.gene_node1[node2][0] == node3:
                            result[4] = node3
                        elif self.gene_node1[node1][0] == node4 or self.gene_node1[node2][0] == node4:
                            result[4] = node4

                        result[1] = len(self.chromosome_in_gene_node1)
                        result[3] = -1
                        result[5] = -1
                        result[6] = 3
                        result[7] = 1

                        return result

            for i in range(len(self.chromosome_in_gene_node2)):
                node1 = self.chromosome_in_gene_node2[i][0]
                node2 = self.chromosome_in_gene_node2[i][len(self.chromosome_in_gene_node2[i]) - 1]

                chromosome_node1: int = self.gene_node1[node1][1]
                node3: int = self.gene_node1[node1][0]

                if node3 != -1 and (chromosome_node1 == which_chromosome or
                                    which_chromosome < 0 and self.touched_chromosome[chromosome_node1] == 0):
                    result[0] = self.gene_node1[node1][1]
                    result[1] = len(self.chromosome_in_gene_node1)
                    result[2] = node1
                    result[3] = -1
                    result[4] = node3
                    result[5] = -1
                    result[6] = 3
                    result[7] = 1

                    return result

                chromosome_node2: int = self.gene_node1[node2][1]
                node4: int = self.gene_node1[node2][0]

                if node4 != -1 and (chromosome_node2 == which_chromosome or
                                    which_chromosome < 0 and self.touched_chromosome[chromosome_node2] == 0):
                    result[0] = self.gene_node1[node2][1]
                    result[1] = len(self.chromosome_in_gene_node1)
                    result[2] = node2
                    result[3] = -1
                    result[4] = node4
                    result[5] = -1
                    result[6] = 3
                    result[7] = 1

                    return result
        elif operation_type == OperationsEnum.FUSION:
            for i in range(len(self.chromosome_in_gene_node1) - 1):
                node1: int = self.chromosome_in_gene_node1[i][0]
                node2: int = self.chromosome_in_gene_node1[i][len(self.chromosome_in_gene_node1[i]) - 1]

                for j in range(i + 1, len(self.chromosome_in_gene_node1)):
                    node3: int = self.chromosome_in_gene_node1[j][0]
                    node4: int = self.chromosome_in_gene_node1[j][len(self.chromosome_in_gene_node1[j]) - 1]

                    if (which_chromosome < 0 and self.touched_chromosome[i] == 0 and
                            self.touched_chromosome[j] == 0) or \
                       (which_chromosome > -1 and (i == which_chromosome or j == which_chromosome)):
                        if self.gene_node2[node1][0] == node3:
                            result[2] = node1
                            result[3] = node3
                            result[7] = 1
                        elif self.gene_node2[node1][0] == node4:
                            result[2] = node1
                            result[3] = node4
                            result[7] = 2
                        elif self.gene_node2[node2][0] == node3:
                            result[2] = node2
                            result[3] = node3
                            result[7] = 3
                        elif self.gene_node2[node2][0] == node4:
                            result[2] = node2
                            result[3] = node4
                            result[7] = 4

                        result[0] = i
                        result[1] = j
                        result[4] = -1
                        result[5] = -1
                        result[6] = 4

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
                        if self.gene_node2[node1][0] != -1 and self.gene_node2[node3][0] != -1:
                            result[2] = node1
                            result[3] = node3
                            result[7] = 1
                        elif self.gene_node2[node1][0] != -1 and self.gene_node2[node4][0] != -1:
                            result[2] = node1
                            result[3] = node4
                            result[7] = 2
                        elif self.gene_node2[node2][0] != -1 and self.gene_node2[node3][0] != -1:
                            result[2] = node2
                            result[3] = node3
                            result[7] = 3
                        elif self.gene_node2[node2][0] != -1 and self.gene_node2[node4][0] != -1:
                            result[2] = node2
                            result[3] = node4
                            result[7] = 4

                        result[0] = i
                        result[1] = j
                        result[4] = -1
                        result[5] = -1
                        result[6] = 4

                        return result

        elif operation_type == OperationsEnum.INVERSION:
            if which_chromosome > -1:
                operation: [int] = self.find_reversal(which_chromosome)

                if operation[0] != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome)):
                    if self.touched_chromosome[c] == 0:
                        operation: [int] = self.find_reversal(c)

                        if operation[0] != -1:
                            return operation
        elif operation_type == OperationsEnum.TRANSLOCATION:
            if which_chromosome > -1:
                operation: [int] = self.find_translocation(which_chromosome, which_chromosome)

                if operation[0] != -1:
                    return operation

            if which_chromosome < 0:
                for c in range(len(self.touched_chromosome) - 1):
                    if self.touched_chromosome[c] == 0:
                        operation: [int] = self.find_translocation(c, which_chromosome)

                        if operation[0] != -1:
                            return operation

        return result

    def find_reversal(self, chromosome: int) -> [int]:
        result: [int] = [0] * 8
        result[0] = -1

        for i in range(len(self.chromosome_in_gene_node1[chromosome]) // 2):
            node1: int = -1

            if i != 0:
                node1 = self.chromosome_in_gene_node1[chromosome][(2 * i) - 1]

            node2: int = self.chromosome_in_gene_node1[chromosome][2 * i]

            if node1 != -1 and self.gene_node2[node1][0] != node2 or \
               node2 != -1 and self.gene_node2[node2][0] != node1:
                node3: int = -1

                if node1 != -1:
                    node3 = self.gene_node2[node1][0]

                node4: int = self.gene_node2[node2][0]

                if (node3 != -1 and self.gene_node1[node3][1] == chromosome) and \
                   (self.gene_node1[node3][2] > (2 * i)):
                    node5: int = self.gene_node1[node3][0]

                    if node5 == -1 or self.gene_node1[node3][2] < self.gene_node1[node5][2]:
                        result[0] = chromosome
                        result[1] = chromosome
                        result[2] = node1
                        result[3] = node3
                        result[4] = node2
                        result[5] = node5
                        result[6] = 1
                        result[7] = 1

                        return result

                if (node4 != -1 and self.gene_node1[node4][1] == chromosome) and \
                   (self.gene_node1[node4][2] > (2 * i)):
                    node6: int = self.gene_node1[node4][0]

                    if node6 != -1 or self.gene_node1[node6][2] < self.gene_node1[node4][2]:
                        result[0] = chromosome
                        result[1] = chromosome
                        result[2] = node1
                        result[3] = node6
                        result[4] = node2
                        result[5] = node4
                        result[6] = 1
                        result[7] = 1

                        return result

        return result

    def find_translocation(self, chromosome: int, which_chromosome: int) -> [int]:
        result: [int] = [0] * 8
        result[0] = -1

        for i in range((len(self.chromosome_in_gene_node1[chromosome]) // 2) + 1):
            node1: int = -1

            if i != 0:
                node1 = self.chromosome_in_gene_node1[chromosome][(2 * i) - 1]

            node2: int = -1

            if i != len(self.chromosome_in_gene_node1[chromosome]) // 2:
                node2 = self.chromosome_in_gene_node1[chromosome][2 * i]

            if node1 != -1 and self.gene_node2[node1][0] != node2 or \
               node2 != -1 and self.gene_node2[node2][0] != node1:
                node3: int = -1

                if node1 != -1:
                    node3 = self.gene_node2[node1][0]

                node4: int = -1

                if node2 != -1:
                    node4 = self.gene_node2[node2][0]

                if (node3 != -1 and self.gene_node1[node3][1] != chromosome) and \
                   (which_chromosome > -1 or self.touched_chromosome[self.gene_node1[node3][1]] == 0):
                    node5: int = self.gene_node1[node3][0]

                    if node2 + node5 != -2:
                        sub_type: int = 2

                        if (self.chromosome_in_gene_node1[self.gene_node1[node3][1]][0] == node3) or \
                           (node5 != -1 and self.gene_node1[node5][2] < self.gene_node1[node3][2]):
                            sub_type = 1

                        result[0] = chromosome
                        result[1] = self.gene_node1[node3][1]
                        result[2] = node1
                        result[3] = node3
                        result[4] = node2
                        result[5] = node5
                        result[6] = 2
                        result[7] = sub_type

                        return result

                if (node4 != -1 and self.gene_node1[node4][1] != chromosome) and \
                   (which_chromosome > -1 or self.touched_chromosome[self.gene_node1[node4][1]] == 0):
                    node6: int = self.gene_node1[node4][0]

                    if node1 + node6 != -2:
                        sub_type: int = 1

                        if (self.chromosome_in_gene_node1[self.gene_node1[node4][1]][0] == node4) or \
                           (node6 != -1 and self.gene_node1[node6][2] < self.gene_node1[node4][2]):
                            sub_type = 2

                        result[0] = chromosome
                        result[1] = self.gene_node1[node4][1]
                        result[2] = node1
                        result[3] = node6
                        result[4] = node2
                        result[5] = node4
                        result[6] = 2
                        result[7] = sub_type

                        return result

        return result
