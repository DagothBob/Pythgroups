from __future__ import annotations

from typing import Optional, Dict

"""                              
 Path for use in the PathGroups algorithm (Median problem) 
                                                           
 Based on PGMPath.java from C.Zheng & D.Sankoff (2011)     
                                                           
 Author: Holger Jensen, Oskar Jensen                                   
"""


def create_pgm_path(head: int, tail: int, genome_head: Optional[int] = None, genome_tail: Optional[int] = None) -> \
                                                                                                        Dict[str, int]:
    """ Creates a dictionary representation of a PGMPath, originally a class in the Java implementation.
    Implemented in dictionary form to improve performance.

    Parameters
    ----------
    head : int
        Head gene
    tail : int
        Tail gene
    genome_head : Optional[int]
        Genome for the head gene
    genome_tail : Optional[int]
        Genome for the tail gene

    Returns
    -------
    Dict[str, int]
        A dictionary representation of a PGMPath
    """
    return {"head": head, "tail": tail, "genome_head": genome_head, "genome_tail": genome_tail}


def connect(path1: Dict[str, int],
            path2: Dict[str, int],
            path_l: Dict[str, int],
            which_genome: Optional[int] = None) -> Optional[Dict[str, int]]:
    """
    Connect two PGMPaths

    Parameters
    ----------
    path1
        First path to connect
    path2
        Second path to connect
    path_l
        Path used for connecting two paths
    which_genome
        Which genome this path belongs to (Omit for GGH)

    Returns
    -------
    Optional[Dict[str, int]]
        New PGMPath from the given paths or None if they are unable to be connected
    """
    head: int
    tail: int
    genome_head: int
    genome_tail: int

    if which_genome is not None:
        if path1["head"] == path_l["head"] and path1["genome_head"] == path_l["genome_head"] and \
                path2["head"] == path_l["tail"] and path2["genome_head"] == path_l["genome_tail"]:
            head = path1["tail"]
            genome_head = path1["genome_tail"]
            tail = path2["tail"]
            genome_tail = path2["genome_tail"]

            if head == 0 and genome_head == -1:
                head = path1["head"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["head"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["head"] == path_l["tail"] and path1["genome_head"] == path_l["genome_tail"] and \
                path2["head"] == path_l["head"] and path2["genome_head"] == path_l["genome_head"]:
            head = path1["tail"]
            genome_head = path1["genome_tail"]
            tail = path2["tail"]
            genome_tail = path2["genome_tail"]

            if head == 0 and genome_head == -1:
                head = path1["head"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["head"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["head"] == path_l["tail"] and path1["genome_head"] == path_l["genome_tail"] and \
                path2["tail"] == path_l["head"] and path2["genome_tail"] == path_l["genome_head"]:
            head = path1["tail"]
            genome_head = path1["genome_tail"]
            tail = path2["head"]
            genome_tail = path2["genome_head"]

            if head == 0 and genome_head == -1:
                head = path1["head"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["tail"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["head"] == path_l["head"] and path1["genome_head"] == path_l["genome_head"] and \
                path2["tail"] == path_l["tail"] and path2["genome_tail"] == path_l["genome_tail"]:
            head = path1["tail"]
            genome_head = path1["genome_tail"]
            tail = path2["head"]
            genome_tail = path2["genome_head"]

            if head == 0 and genome_head == -1:
                head = path1["head"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["tail"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["tail"] == path_l["head"] and path1["genome_tail"] == path_l["genome_head"] and \
                path2["head"] == path_l["tail"] and path2["genome_head"] == path_l["genome_tail"]:
            head = path1["head"]
            genome_head = path1["genome_head"]
            tail = path2["tail"]
            genome_tail = path2["genome_tail"]

            if head == 0 and genome_head == -1:
                head = path1["tail"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["head"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["tail"] == path_l["tail"] and path1["genome_tail"] == path_l["genome_tail"] and \
                path2["head"] == path_l["head"] and path2["genome_head"] == path_l["genome_head"]:
            head = path1["head"]
            genome_head = path1["genome_head"]
            tail = path2["tail"]
            genome_tail = path2["genome_tail"]

            if head == 0 and genome_head == -1:
                head = path1["tail"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["head"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["tail"] == path_l["head"] and path1["genome_tail"] == path_l["genome_head"] and \
                path2["tail"] == path_l["tail"] and path2["genome_tail"] == path_l["genome_tail"]:
            head = path1["head"]
            genome_head = path1["genome_head"]
            tail = path2["head"]
            genome_tail = path2["genome_head"]

            if head == 0 and genome_head == -1:
                head = path1["tail"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["tail"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif path1["tail"] == path_l["tail"] and path1["genome_tail"] == path_l["genome_tail"] and \
                path2["tail"] == path_l["head"] and path2["genome_tail"] == path_l["genome_head"]:
            head = path1["head"]
            genome_head = path1["genome_head"]
            tail = path2["head"]
            genome_tail = path2["genome_head"]

            if head == 0 and genome_head == -1:
                head = path1["tail"]
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = path2["tail"]
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)
    else:  # GGH version
        if path1["head"] == path_l["head"] and path2["head"] == path_l["tail"]:
            return create_pgm_path(path1["tail"], path2["tail"])
        elif path1["head"] == path_l["tail"] and path2["head"] == path_l["head"]:
            return create_pgm_path(path1["tail"], path2["tail"])
        elif path1["head"] == path_l["tail"] and path2["tail"] == path_l["head"]:
            return create_pgm_path(path1["tail"], path2["head"])
        elif path1["head"] == path_l["head"] and path2["tail"] == path_l["tail"]:
            return create_pgm_path(path1["tail"], path2["head"])
        elif path1["tail"] == path_l["head"] and path2["head"] == path_l["tail"]:
            return create_pgm_path(path1["head"], path2["tail"])
        elif path1["tail"] == path_l["tail"] and path2["head"] == path_l["head"]:
            return create_pgm_path(path1["head"], path2["tail"])
        elif path1["tail"] == path_l["head"] and path2["tail"] == path_l["tail"]:
            return create_pgm_path(path1["head"], path2["head"])
        elif path1["tail"] == path_l["tail"] and path2["tail"] == path_l["head"]:
            return create_pgm_path(path1["head"], path2["head"])

    return None
