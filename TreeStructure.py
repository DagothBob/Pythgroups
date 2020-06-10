from typing import Optional

from GenomeInString import GenomeInString
from MedianData import MedianData
from PGMPath import PGMPath
from PGMPathForAGenome import PGMPathForAGenome


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


class TreeStructure:
    """
    Attributes
    ----------
    number_of_ancestors
        Number of ancestor nodes in the tree structure
    number_of_leaves
        Number of leaf nodes in the tree structure
    gene_number
        Number of genes in the structure
    leaves
        List of leaf nodes
    medians
        List of MedianData information
    node_int
        Current node as a list of integers
    node_string
        Current node as a list of strings
    all_paths
        All paths as PGMForAGenomes
    """
    def __init__(self,
                 number_of_ancestors: int,
                 number_of_leaves: int,
                 gene_number: int,
                 paths: [PGMPathForAGenome],
                 node_strings: [str],
                 node_ints: [int],
                 ancestor_genome_string: [Optional[GenomeInString]]):
        """
        Constructor

        Parameters
        ----------
        number_of_ancestors
            Number of ancestor nodes in the tree structure
        number_of_leaves
            Number of leaf nodes in the tree structure
        gene_number
            Number of genes in the structure
        paths
            All paths as PGMForAGenomes
        node_strings
            Current node as a list of strings
        node_ints
            Current node as a list of integers
        ancestor_genome_string
            Ancestor genomes as a list of GenomeInStrings
        """
        self.number_of_ancestors: int = number_of_ancestors
        self.number_of_leaves: int = number_of_leaves
        self.gene_number: int = gene_number
        self.leaves: [[int]] = [[-1] * (self.number_of_ancestors + self.number_of_leaves)] * 3
        self.medians: [MedianData] = [MedianData(None, None, None, 0, 0, 0, 0, 0, None)] * self.number_of_ancestors
        self.node_int: [int] = [0] * (self.gene_number * 2)
        self.node_string: [str] = [" "] * (self.gene_number * 2)

        if ancestor_genome_string is None:
            self.node_int[:] = node_ints.copy()[:]
            self.node_string[:] = node_strings.copy()[:]

            self.all_paths: [PGMPathForAGenome] = [PGMPathForAGenome(None)] * len(paths)

            for i in range(len(paths)):
                self.all_paths[i] = PGMPathForAGenome(paths[i].paths)

            self.all_paths[3] = PGMPathForAGenome(self.get_pgm_path(None, 3))
            self.set_tree_structure(3, 0, 1, 2)
        else:
            genome1: [str] = ancestor_genome_string[0].chromosomes
            index1: int = 0

            for chromosome in genome1:
                genes: [str] = chromosome.split(' ')

                for gene in genes:
                    first_character: str = gene[0:1]
                    node1: str
                    node2: str

                    if first_character == '-':
                        node1 = gene[1:]
                        node1.join('h')
                        node2 = gene[1:]
                        node2.join('t')

                        self.node_int[index1] = index1 + 1
                        self.node_string[index1] = node2
                        index1 += 1
                        self.node_int[index1] = index1 + 1
                        self.node_string[index1] = node1
                    else:
                        node1 = gene
                        node1.join('t')
                        node2 = gene
                        node2.join('h')

                        self.node_int[index1] = index1 + 1
                        self.node_string[index1] = node1
                        index1 += 1
                        self.node_int[index1] = index1 + 1
                        self.node_string[index1] = node2

                    index1 += 1

            self.all_genomes: [GenomeInString] = [GenomeInString(None)] * (
                        self.number_of_leaves + self.number_of_ancestors)

            for i in range(len(ancestor_genome_string)):
                self.all_genomes[i] = GenomeInString(ancestor_genome_string[i].chromosomes)

            self.all_paths: [PGMPathForAGenome] = [PGMPathForAGenome(None)] * (
                        self.number_of_leaves + self.number_of_ancestors)

            for i in range(len(ancestor_genome_string)):
                self.all_paths[i] = PGMPathForAGenome(self.get_pgm_path(self.all_genomes[i].chromosomes, i))

            for i in range(len(ancestor_genome_string), len(self.all_genomes)):
                self.all_paths[i] = PGMPathForAGenome(self.get_pgm_path(None, i))

    def set_tree_structure(self, which_genome: int, genome1: int, genome2: int, genome3: int):
        """
        Sets the leaves and medians of the tree structure

        Parameters
        ----------
        which_genome
            Which genome this belongs to
        genome1
            First genome to fill tree leaves with
        genome2
            Second genome to fill tree leaves with
        genome3
            Third genome to fill tree leaves with
        """
        self.leaves[which_genome][0] = genome1
        self.leaves[which_genome][1] = genome2
        self.leaves[which_genome][2] = genome3
        self.medians[which_genome - self.number_of_leaves] = MedianData(self.all_paths[genome1],
                                                                        self.all_paths[genome2],
                                                                        self.all_paths[genome3],
                                                                        self.gene_number,
                                                                        which_genome,
                                                                        genome1,
                                                                        genome2,
                                                                        genome3,
                                                                        self.node_string)

    def get_relation(self) -> [[int]]:
        """
        Returns
        -------
        Relation between leaves
        """
        relation: [[int]] = [[0] * (self.number_of_leaves + self.number_of_ancestors)] * (
                    self.number_of_leaves + self.number_of_ancestors)

        for i in range(len(self.leaves)):
            for j in range(len(self.leaves[i])):
                if self.leaves[i][j] != -1:
                    relation[i][self.leaves[i][j]] = 2

        return relation

    def get_pgm_path(self, genome: [Optional[str]], which_genome: int) -> [PGMPath]:
        """
        Gets list of PGMPaths for a genome

        Parameters
        ----------
        genome
            Genome to get the PGMPath for
        which_genome
            Identifier for the genome

        Returns
        -------
        [PGMPath]
            List of PGMPaths for the genome
        """
        if genome is None:
            path2: [PGMPath] = [PGMPath(0, 0, 0, 0)] * ((2 * self.gene_number) + 1)

            for i in range(1, len(path2)):
                path2[i] = PGMPath(i, 0, which_genome, -1)

            return path2

        path1: [PGMPath] = [PGMPath(0, 0, 0, 0)] * ((2 * self.gene_number) + 1)
        null_node: int = -1

        for chromosome in genome:
            genes: [str] = chromosome.split(' ')
            pre_node: int = 0

            for j in range(len(genes)):
                first_character: str = genes[j][0:1]
                node1: str
                node2: str

                if first_character == '-':
                    node1 = genes[j][1:]
                    node1.join('h')
                    node2 = genes[j][1:]
                    node2.join('t')
                else:
                    node1 = genes[j]
                    node1.join('t')
                    node2 = genes[j]
                    node2.join('h')

                node1_int: int = self.find_node_int(node1)
                node2_int: int = self.find_node_int(node2)

                if node1_int == 0 or node2_int == 0:
                    print("Gene ", str(genes[j], " does not exist in the other genome.\n"))

                if j == 0:
                    path1[node1_int] = PGMPath(node1_int, null_node, which_genome, which_genome)
                    pre_node = node2_int
                    null_node -= 1
                elif j != 0 and j != len(genes) - 1:
                    path1[node1_int] = PGMPath(node1_int, pre_node, which_genome, which_genome)
                    path1[pre_node] = PGMPath(pre_node, node1_int, which_genome, which_genome)
                    pre_node = node2_int
                elif j == len(genes) - 1:
                    if len(genes) != 1:
                        path1[pre_node] = PGMPath(pre_node, node1_int, which_genome, which_genome)
                        path1[node1_int] = PGMPath(node1_int, pre_node, which_genome, which_genome)

                    path1[node2_int] = PGMPath(node2_int, null_node, which_genome, which_genome)
                    null_node -= 1

        return path1

    def find_node_int(self, ancestor_string: str) -> int:
        """
        Gets the integer node identifier for an ancestor node

        Parameters
        ----------
        ancestor_string
            Ancestor node represented as a string

        Returns
        -------
        int
            Identifier for the node
        """
        for i in range(len(self.node_string)):
            if self.node_string[i] == ancestor_string:
                return i + 1

        return 0
