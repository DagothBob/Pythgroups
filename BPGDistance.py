import BPGPath
from BPGPath import BPGPath


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Calculate the distance between two genomes using          #
# a breakpoint graph and counting the number of cycles.     #
#                                                           #
# Based on BPGDistance.java from C.Zheng & D.Sankoff (2011) #
#                                                           #
# Author: Holger Jensen                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def insert_character(s: str, i: int, c: str) -> str:
    return s[:i] + c + s[(i + 1):]


def is_valid_cycle(edge1: int, edge2: int):
    if edge1 >= 0 or edge2 >= 0:
        return False

    if (edge1 / 2) * 2 != edge1:
        return True

    return (edge2 / 2) * 2 != edge2


class BPGDistance:
    def __init__(self, genome1: [str], genome2: [str]):
        index1: int = 0
        index2: int = 0
        num_genes1: int = 0
        num_genes2: int = 0
        self.genome1: str = " " * len(genome1)                # Genome 1 as String
        self.genome2: str = " " * len(genome2)                # Genome 2 as String

        for gene in genome1:
            genes: [str] = gene.split(' ')

            if len(genes) != 0:
                num_genes1 += len(genes)
                self.genome1 = insert_character(self.genome1, index1, gene)
                index1 += 1

        for gene in genome2:
            genes: [str] = gene.split(' ')

            if len(genes) != 0:
                num_genes2 += len(genes)
                self.genome2 = insert_character(self.genome2, index2, gene)
                index2 += 1

        if num_genes1 == num_genes2:                          # Algorithm requires genomes are equal length
            self.gene_number: int = num_genes1                # Number of genes

        self.gene_ints: [int] = [0] * (self.gene_number * 2)  # Integer values for genes
        self.gene_strs1: str = " " * (self.gene_number * 2)   # String representation of genes (genome 1)
        self.gene_strs2: str = " " * (self.gene_number * 2)   # String representation of genes (genome 2)
        self.genome_paths1: [BPGPath] = []                    # Paths for combining into cycles (genome 1)
        self.genome_paths2: [BPGPath] = []                    # Paths for combining into cycles (genome 2)
        self.distance: int = 0                                # Distance

    # Initialize BPG graph to the state which contains the set of gene/telomere paths
    def graph_init(self):
        index1: int = 0

        for g_string in self.genome1:
            genes: [str] = g_string.split(' ')

            for gene in genes:
                first_character: str = gene[0]  # Sign indicating gene is head-tail or tail-head

                if first_character == "-":
                    node1: str = gene[1:]
                    node1.join("h")
                    node2: str = gene[1:]
                    node2.join("t")
                    self.gene_ints[index1] = index1 + 1
                    self.gene_strs1 = insert_character(self.gene_strs1, index1, node2)

                    index1 += 1
                    self.gene_ints[index1] = index1 + 1
                    self.gene_strs1 = insert_character(self.gene_strs1, index1, node1)
                else:
                    node1: str = gene
                    node1.join("t")
                    node2: str = gene
                    node2.join("h")
                    self.gene_ints[index1] = index1 + 1
                    self.gene_strs1 = insert_character(self.gene_strs1, index1, node1)

                    index1 += 1
                    self.gene_ints[index1] = index1 + 1
                    self.gene_strs1 = insert_character(self.gene_strs1, index1, node2)

                index1 += 1

        self.gene_strs2 = self.gene_strs1

        self.genome_paths1 = self.get_path(self.genome1, 1)
        self.genome_paths2 = self.get_path(self.genome2, 2)

    # Gets a gene/telomere path
    def get_path(self, genome: [str], node: int) -> [BPGPath]:
        null_node: int = -node
        path1: [BPGPath] = [BPGPath(None, None, None, None)] * ((self.gene_number * 2) + 1)

        for gene in genome:
            genes: [str] = gene.split(' ')
            pre_node: int = 0

            for i in range(len(genes)):
                first_character: str = genes[i][0, 1]  # Sign indicating gene is head-tail or tail-head

                if first_character == "-":
                    node1: str = genes[i][1:]
                    node1.join('h')
                    node2: str = genes[i][1:]
                    node2.join('t')
                else:
                    node1: str = genes[i]
                    node1.join('t')
                    node2: str = genes[i]
                    node2.join('h')

                node1_int: int = self.find_node_int(node1, node)
                node2_int: int = self.find_node_int(node2, node)

                if i == 0:
                    path1[node1_int] = BPGPath(node1_int, null_node, node, node)
                    pre_node: int = node2_int
                    null_node -= 2

                if i != 0 and i != len(genes) - 1:
                    path1[node1_int] = BPGPath(node1_int, pre_node, node, node)
                    path1[pre_node] = BPGPath(pre_node, node1_int, node, node)
                    pre_node: int = node2_int

                if i == len(genes) - 1:
                    if len(genes) != 1:
                        path1[pre_node] = BPGPath(pre_node, node1_int, node, node)
                        path1[node1_int] = BPGPath(node1_int, pre_node, node, node)

                    path1[node2_int] = BPGPath(node2_int, null_node, node, node)
                    null_node -= 2

        return path1

    # Gets the node's int value from the full string representation
    def find_node_int(self, gene_str: str, gene_int: int) -> int:
        if gene_int == 1:
            for i in range(len(self.gene_strs1)):
                if self.gene_strs1[i] != " " and self.gene_strs1 == gene_str:
                    self.gene_strs1 = insert_character(self.gene_strs1, i, " ")
                    return i + 1

        if gene_int == 2:
            for i in range(len(self.gene_strs2)):
                if self.gene_strs2[i] != " " and self.gene_strs2 == gene_str:
                    self.gene_strs2 = insert_character(self.gene_strs2, i, " ")
                    return i + 1

        return 0

    # Produce cycles from the gene/telomere paths and count them to calculate genome distance
    def calculate_distance(self):
        self.graph_init()

        cycle_number: int = 0
        good_path_number: int = 0

        # Get the first path out of a list of paths with None values filtered out
        a_path1: BPGPath = list(filter(lambda ps: True if ps is not None else False, self.genome_paths1))[0]

        while a_path1 is not None:
            node1: int = a_path1.head
            node2: int = a_path1.tail
            start: int = node1

            if node1 > 0:
                self.genome_paths1[node1] = None

            if node2 > 0:
                self.genome_paths1[node2] = None

            if node1 < node2:
                node_big: int = node2
                node_small: int = node1
            else:
                node_big: int = node1
                node_small: int = node2

            more: bool = True

            while more:
                if a_path1.head < 0 and a_path1.tail < 0 and is_valid_cycle(a_path1.head, a_path1.tail):
                    good_path_number += 1
                    break

                l_path: BPGPath = self.genome_paths2[node_big]
                l_node1: int = l_path.head
                l_node2: int = l_path.tail

                if l_node2 > 0 and l_node2 == node_small:
                    cycle_number += 1
                    more = False
                    self.genome_paths2[l_node1] = None
                    self.genome_paths2[l_node2] = None
                if l_node2 > 0 and l_node2 != node_small:
                    a_path2: BPGPath = self.genome_paths1[l_node2]
                    ap2_tail: int = a_path2.tail

                    a_path1 = BPGPath(ap2_tail, node_small, 1, 1)
                    node1 = a_path1.head
                    node2 = a_path1.tail

                    if node1 < node2:
                        node_big = node2
                        node_small = node1
                    else:
                        node_big = node1
                        node_small = node2

                    paths1_tail: int = self.genome_paths1[l_node2].tail

                    if paths1_tail > 0:
                        self.genome_paths1[paths1_tail] = None

                    self.genome_paths2[l_node1] = None
                    self.genome_paths2[l_node2] = None
                    self.genome_paths1[l_node2] = None
                if l_node2 < 0 and node_small < 0:
                    if is_valid_cycle(l_node2, node_small):
                        good_path_number += 1

                    more = False
                    self.genome_paths2[l_node1] = None
                if l_node2 < 0 < node_small:
                    a_path1 = BPGPath(node_small, l_node2, 1, 2)
                    node1 = a_path1.head
                    node2 = a_path1.tail

                    if node1 < node2:
                        node_big = node2
                        node_small = node1
                    else:
                        node_big = node1
                        node_small = node2

                    self.genome_paths2[l_node1] = None

            a_path1 = list(filter(lambda ps: True if ps is not None else False, self.genome_paths1))[start]

        self.distance = self.gene_number + len(self.genome1) - cycle_number - good_path_number
