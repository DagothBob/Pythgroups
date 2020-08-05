from __future__ import annotations

from enum import IntEnum

"""
Gene object as part of a Chromosome

Based on Gene.java from C.Zheng & D.Sankoff (2011)

Author: Holger Jensen, Oskar Jensen
"""


class TailOrHead(IntEnum):
    NULL: int = -1
    TAIL: int = 1
    HEAD: int = 2


class Gene:
    """
    Attributes
    ----------
    name: str
        Name of the gene
    node_1: int
        First node (Tail or head)
    node_2: int
        Second node (Tail or head)
    """
    def __init__(self, name: str, node1: int, node2: int):
        """
        Constructor

        Parameters
        ----------
        name
            Name of the gene
        node1
            First node (Tail or head)
        node2
            Second node (Tail or head)
        """
        self.name: str = name
        self.node_1: int = node1
        self.node_2: int = node2

    @classmethod
    def default(cls) -> Gene:
        """
        Construct with default values

        Returns
        -------
        Gene
            Gene object with default values for name, node_1, node_2
        """
        return cls("", TailOrHead.NULL, TailOrHead.NULL)

    @classmethod
    def with_name(cls, name: str) -> Gene:
        """
        Construct with name

        Parameters
        ----------
        name
            Name of the gene

        Returns
        -------
        Gene
            Gene object with name and node values based on parity of the gene
        """
        if name[0] == "-":
            return cls(name, TailOrHead.HEAD, TailOrHead.TAIL)
        else:
            return cls(name, TailOrHead.HEAD, TailOrHead.TAIL)

    @classmethod
    def from_gene(cls, gene: Gene) -> Gene:
        """
        Copy constructor from existing Gene object

        Parameters
        ----------
        gene
            Gene to copy attributes from

        Returns
        -------
        Gene
            New Gene instance copy
        """
        return cls(gene.name, gene.node_1, gene.node_2)

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        str
            String representation of the object's attributes
        """
        return self.name  # + "\t" + str(self.node_1) + "\t" + str(self.node_2)
