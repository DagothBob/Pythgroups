from __future__ import annotations

from typing import Optional, Dict, Any

"""                                 
 Used in MedianData for solving the rearrangement median problem.
                                                                 
 Based on ChoiceStructure.java from C.Zheng & D.Sankoff (2011)   
                                                                 
 Author: Holger Jensen, Oskar Jensen                             
"""


def create_cs(source: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """ Creates a dictionary representation of a Choice Structure, originally a class in the Java implementation.
    Implemented in dictionary form to improve performance.
    Can be created either with default values or based on an existing Choice Structure.

    Parameters
    ----------
    source : Optional[Dict[str, Any]]
        Optional choice structure to base this copy off of. If left blank, will instantiate with default values.

    Returns
    -------
    Dict[str, Any]
        Newly created Choice Structure
    """
    cs: Dict[str, Any] = dict()
    cs["index_from"]: int
    cs["for_which_genome"]: int
    cs["priority"]: int
    cs["position"]: int
    cs["genome_1_path"]: Optional[Dict[str, int]]
    cs["genome_2_path"]: Optional[Dict[str, int]]
    cs["genome_3_path"]: Optional[Dict[str, int]]
    cs["gray_edge"]: Optional[Dict[str, int]]

    if source is None:
        cs["index_from"], \
            cs["for_which_genome"], \
            cs["priority"], \
            cs["position"], \
            cs["genome_1_path"], \
            cs["genome_2_path"], \
            cs["genome_3_path"], \
            cs["gray_edge"] = int(), int(), int(), int(), None, None, None, None
    else:
        cs["index_from"], \
            cs["for_which_genome"], \
            cs["priority"], \
            cs["position"], \
            cs["genome_1_path"], \
            cs["genome_2_path"], \
            cs["genome_3_path"], \
            cs["gray_edge"] = source["index_from"], source["for_which_genome"], source["priority"], source["position"],\
            source["genome_1_path"], source["genome_2_path"], source["genome_3_path"], source["gray_edge"]

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
