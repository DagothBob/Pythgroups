from __future__ import annotations

from typing import Optional, Dict, Any

"""                                 
 Used in MedianData for solving the rearrangement median problem.
                                                                 
 Based on ChoiceStructure.java from C.Zheng & D.Sankoff (2011)   
                                                                 
 Author: Holger Jensen, Oskar Jensen                             
"""


class ChoiceStructure:
    """
    Attributes
    ----------
    index_from : int
        Object instance identifier
    for_which_genome : int
        For which genome this structure refers to
    priority : int
        Priority level for makeCycles algorithm
    position : int
        Position of the choice structure in Priority's cs_index
    genome_1_path : Optional[Dict[str, int]]
        Path for genome 1
    genome_2_path : Optional[Dict[str, int]]
        Path for genome 2
    genome_3_path : Optional[Dict[str, int]]
        Path for genome 3
    gray_edge : Optional[Dict[str, int]]
        Used when updating paths and priorities in the makeCycles algorithm
    """

    def __init__(self):
        """
        Constructor
        """
        self.index_from: int = int()
        self.for_which_genome: int = int()
        self.priority: int = int()
        self.position: int = int()
        self.genome_1_path: Optional[Dict[str, int]] = None
        self.genome_2_path: Optional[Dict[str, int]] = None
        self.genome_3_path: Optional[Dict[str, int]] = None
        self.gray_edge: Optional[Dict[str, int]] = None

    def from_cs(self, cs: ChoiceStructure):
        """
        Constructs a new ChoiceStructure from an existing given one

        Parameters
        ----------
        cs
            ChoiceStructure to copy from
        """
        self.index_from = cs.index_from
        self.for_which_genome = cs.for_which_genome
        self.priority = cs.priority
        self.position = cs.position
        self.genome_1_path = cs.genome_1_path
        self.genome_2_path = cs.genome_2_path
        self.genome_3_path = cs.genome_3_path
        self.gray_edge = cs.gray_edge

    def set_new_path(self, path: Dict[str, int], ploidy: Optional[int] = None, gene_number: Optional[int] = None):
        """
        Sets instance paths to the given PGMPath if its genome matches

        Parameters
        ----------
        path
            PGMPath to copy from
        ploidy
            Monoploid or diploid
        gene_number
            Number of genes
        """
        genome_here: int

        if ploidy is None:
            genome_here = path["genome_head"]

            if genome_here == self.genome_1_path["genome_head"]:
                self.genome_1_path = path

            if genome_here == self.genome_2_path["genome_head"]:
                self.genome_2_path = path

            if genome_here == self.genome_3_path["genome_head"]:
                self.genome_3_path = path
        else:
            genome_here = path["head"]

            if genome_here > gene_number * 2:
                genome_here -= gene_number * 2

            if self.index_from != genome_here:
                raise Exception(
                    "Object instance attribute index_from is not equal to from in ChoiceStructure.set_new_path().\n")

            if ploidy == 1:
                self.genome_3_path = path
            elif ploidy == 2:
                if self.index_from == path["head"]:
                    self.genome_1_path = path

                if self.index_from + gene_number * 2 == path["head"]:
                    self.genome_2_path = path


def create_cs(source: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """ Creates a dictionary representation of a Choice Structure, either with default values or based on an existing
    Choice Structure. Implemented in dictionary form to improve performance.

    Parameters
    ----------
    source : Optional[Dict[str, Any]]
        Optional choice structure to base this copy off of

    Returns
    -------
    Dict[str, Any]
        Newly created Choice Structure
    """
    cs: Dict[str, Any] = dict()
    if source is None:
        cs["index_from"]: int = int()
        cs["for_which_genome"]: int = int()
        cs["priority"]: int = int()
        cs["position"]: int = int()
        cs["genome_1_path"]: Optional[Dict[str, int]] = None
        cs["genome_2_path"]: Optional[Dict[str, int]] = None
        cs["genome_3_path"]: Optional[Dict[str, int]] = None
        cs["gray_edge"]: Optional[Dict[str, int]] = None
    else:
        cs["index_from"]: int = source["index_from"]
        cs["for_which_genome"]: int = source["for_which_genome"]
        cs["priority"]: int = source["priority"]
        cs["position"]: int = source["position"]
        cs["genome_1_path"]: Optional[Dict[str, int]] = source["genome_1_path"]
        cs["genome_2_path"]: Optional[Dict[str, int]] = source["genome_2_path"]
        cs["genome_3_path"]: Optional[Dict[str, int]] = source["genome_3_path"]
        cs["gray_edge"]: Optional[Dict[str, int]] = source["gray_edge"]

    return cs


def set_new_path(source: Dict[str, Any], path: Dict[str, int],
                 ploidy: Optional[int] = None, gene_number: Optional[int] = None) -> Optional[Dict[str, Any]]:
    """
    Sets instance paths to the given PGMPath if its genome matches

    Parameters
    ----------
    source
        The source choice structure to modify
    path
        PGMPath to copy from
    ploidy
        Monoploid or diploid
    gene_number
        Number of genes
    """
    genome_here: int

    if ploidy is None:
        genome_here = path["genome_head"]

        if genome_here == source["genome_1_path"]["genome_head"]:
            source["genome_1_path"] = path

        if genome_here == source["genome_2_path"]["genome_head"]:
            source["genome_2_path"] = path

        if genome_here == source["genome_3_path"]["genome_head"]:
            source["genome_3_path"] = path
    else:
        genome_here = path["head"]

        if genome_here > gene_number * 2:
            genome_here -= gene_number * 2

        if source["index_from"] != genome_here:
            raise Exception(
                "Object instance attribute index_from is not equal to from in ChoiceStructure.set_new_path().\n")

        if ploidy == 1:
            source["genome_3_path"] = path
        elif ploidy == 2:
            if source["index_from"] == path["head"]:
                source["genome_1_path"] = path

            if source["index_from"] + gene_number * 2 == path["head"]:
                source["genome_2_path"] = path

    return source
