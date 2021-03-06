from typing import List

from Chromosome import Chromosome

"""
Object for representing a genome as a list of Chromosomes

Based on Genome.java from C.Zheng & D.Sankoff (2011)

Author: Holger Jensen, Oskar Jensen
"""


def split_at_whitespace(strings: str) -> List[str]:
    """
    Strips, then splits, a string, then strips the substrings again

    Parameters
    ----------
    strings
        List of strings to operate on

    Returns
    -------
    List[str]
        Set of cleaned-up strings
    """
    return [s.strip() for s in strings.strip().split(" ") if s.strip() != ""]


class Genome:
    """
    Attributes
    ----------
    chromosomes: List[Chromosome]
        List of chromosomes belonging to this genome
    """
    def __init__(self, chromosomes: List[Chromosome]):
        """
        Constructor

        Parameters
        ----------
        chromosomes
            List of chromosomes making up the genome
        """
        self.chromosomes: List[Chromosome] = chromosomes

    @classmethod
    def from_strings(cls, g_string: List[str]) -> "Genome":
        """
        Construct a Genome from a list of strings representing chromosomes

        Parameters
        ----------
        g_string
            String representing the genome

        Returns
        -------
        Genome
            Genome instance with the specified chromosomes
        """
        return cls([Chromosome.from_strings(split_at_whitespace(s)) for s in g_string])

    def add_chromosome(self, chromosome: Chromosome):
        """
        Adds a new chromosome to the genome object

        Parameters
        ----------
        chromosome
            Chromosome to be appended to the genome
        """
        temp: List[Chromosome] = list(self.chromosomes)
        temp.append(chromosome)

        self.chromosomes = temp

    def remove_chromosome_at_index(self, index: int):
        """
        Remove chromosome at the given index

        Parameters
        ----------
        index
            Index of the chromosome to be removed
        """
        if index < len(self.chromosomes):
            del(self.chromosomes[index])

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        str
            String representation of the Genome object
        """
        ret: str = str()

        for chromosome in self.chromosomes:
            ret += str(chromosome) + "\n"

        return ret
