from typing import List

from Gene import Gene

"""
Chromosome object as part of a Genome

Based on Chromosome.java from C.Zheng & D.Sankoff (2011)

Author: Holger Jensen, Oskar Jensen
"""


class Chromosome:
    """
    Attributes
    ----------
    genes: List[Gene]
        Genes which make up this chromosome
    """
    def __init__(self, genes: List[Gene]):
        """
        Constructor

        Parameters
        ----------
        genes
            Genes which make up this chromosome
        """
        self.genes: List[Gene] = genes

    @classmethod
    def from_chromosome(cls, chromosome: "Chromosome") -> "Chromosome":
        """
        Copy constructor

        Parameters
        ----------
        chromosome
            Chromosome object to copy attribute values from

        Returns
        -------
        Chromosome
            New Chromosome with data copied from supplied object instance
        """
        return cls(chromosome.genes)

    @classmethod
    def from_strings(cls, gene_strings: List[str]) -> "Chromosome":
        """
        Constructor from list of gene names

        Parameters
        ----------
        gene_strings
            List of gene names to construct Gene objects with

        Returns
        -------
        Chromosome
            New Chromosome with constructed genes as object data
        """
        return cls([Gene.with_name(s) for s in gene_strings])

    def to_strings(self) -> List[str]:
        """
        Convert Chromosome object into string representation for backwards compatibility

        Returns
        -------
        List[str]
            String list representation of the chromosome
        """
        return [gene.name for gene in self.genes]

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        String representation of the Chromosome object
        """
        ret: str = "chr: \n"

        for gene in self.genes:
            ret += str(gene) + " "

        return ret
