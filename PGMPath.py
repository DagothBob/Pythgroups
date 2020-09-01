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
    Connect two PGMPaths by placing path1 at the head of the new path
    and path2 at the tail, while maintaining path1's relative facing
    direction to path2

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
    
    p1_head: int = path1["head"]
    p1_tail: int = path1["tail"]
    p2_head: int = path2["head"]
    p2_tail: int = path2["tail"]
    p_l_head: int = path_l["head"]
    p_l_tail: int = path_l["tail"]

    if which_genome is not None:
        p1_ghead: int = path1["genome_head"]
        p1_gtail: int = path1["genome_tail"]
        p2_ghead: int = path2["genome_head"]
        p2_gtail: int = path2["genome_tail"]
        pl_ghead: int = path_l["genome_head"]
        pl_gtail: int = path_l["genome_tail"]

        if p1_head == p_l_head and p1_ghead == pl_ghead and \
                p2_head == p_l_tail and p2_ghead == pl_gtail:
            head = p1_tail
            genome_head = p1_gtail
            tail = p2_tail
            genome_tail = p2_gtail

            if head == 0 and genome_head == -1:
                head = p1_head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_head
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_head == p_l_tail and p1_ghead == pl_gtail and \
                p2_head == p_l_head and p2_ghead == pl_ghead:
            head = p1_tail
            genome_head = p1_gtail
            tail = p2_tail
            genome_tail = p2_gtail

            if head == 0 and genome_head == -1:
                head = p1_head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_head
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_head == p_l_tail and p1_ghead == pl_gtail and \
                p2_tail == p_l_head and p2_gtail == pl_ghead:
            head = p1_tail
            genome_head = p1_gtail
            tail = p2_head
            genome_tail = p2_ghead

            if head == 0 and genome_head == -1:
                head = p1_head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_tail
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_head == p_l_head and p1_ghead == pl_ghead and \
                p2_tail == p_l_tail and p2_gtail == pl_gtail:
            head = p1_tail
            genome_head = p1_gtail
            tail = p2_head
            genome_tail = p2_ghead

            if head == 0 and genome_head == -1:
                head = p1_head
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_tail
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_tail == p_l_head and p1_gtail == pl_ghead and \
                p2_head == p_l_tail and p2_ghead == pl_gtail:
            head = p1_head
            genome_head = p1_ghead
            tail = p2_tail
            genome_tail = p2_gtail

            if head == 0 and genome_head == -1:
                head = p1_tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_head
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_tail == p_l_tail and p1_gtail == pl_gtail and \
                p2_head == p_l_head and p2_ghead == pl_ghead:
            head = p1_head
            genome_head = p1_ghead
            tail = p2_tail
            genome_tail = p2_gtail

            if head == 0 and genome_head == -1:
                head = p1_tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_head
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_tail == p_l_head and p1_gtail == pl_ghead and \
                p2_tail == p_l_tail and p2_gtail == pl_gtail:
            head = p1_head
            genome_head = p1_ghead
            tail = p2_head
            genome_tail = p2_ghead

            if head == 0 and genome_head == -1:
                head = p1_tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_tail
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)

        elif p1_tail == p_l_tail and p1_gtail == pl_gtail and \
                p2_tail == p_l_head and p2_gtail == pl_ghead:
            head = p1_head
            genome_head = p1_ghead
            tail = p2_head
            genome_tail = p2_ghead

            if head == 0 and genome_head == -1:
                head = p1_tail
                genome_head = which_genome

            if tail == 0 and genome_tail == -1:
                tail = p2_tail
                genome_tail = which_genome

            return create_pgm_path(head, tail, genome_head, genome_tail)
    else:  # GGH/ALQ version
        if p1_head == p_l_head and p2_head == p_l_tail:
            return create_pgm_path(p1_tail, p2_tail)
        elif p1_head == p_l_tail and p2_head == p_l_head:
            return create_pgm_path(p1_tail, p2_tail)
        elif p1_head == p_l_tail and p2_tail == p_l_head:
            return create_pgm_path(p1_tail, p2_head)
        elif p1_head == p_l_head and p2_tail == p_l_tail:
            return create_pgm_path(p1_tail, p2_head)
        elif p1_tail == p_l_head and p2_head == p_l_tail:
            return create_pgm_path(p1_head, p2_tail)
        elif p1_tail == p_l_tail and p2_head == p_l_head:
            return create_pgm_path(p1_head, p2_tail)
        elif p1_tail == p_l_head and p2_tail == p_l_tail:
            return create_pgm_path(p1_head, p2_head)
        elif p1_tail == p_l_tail and p2_tail == p_l_head:
            return create_pgm_path(p1_head, p2_head)

    return None
