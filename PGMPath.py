from __future__ import annotations

from typing import Optional, Dict

"""                              
 Path for use in the PathGroups algorithm (Median problem) 
                                                           
 Based on PGMPath.java from C.Zheng & D.Sankoff (2011)     
                                                           
 Author: Holger Jensen, Oskar Jensen                                   
"""


# class PGMPath:
#     """
#     Attributes
#     ----------
#     head : int
#         Head gene
#     tail : int
#         Tail gene
#     genome_head : int
#         Genome for the head gene
#     genome_tail : int
#         Genome for the tail gene
#     """
#
#     def __init__(self, head: int, tail: int, ghead: Optional[int] = None, gtail: Optional[int] = None):
#         """
#         Constructor
#
#         Parameters
#         ----------
#         head
#             Head gene
#         tail
#             Tail gene
#         ghead
#             Genome for the head gene
#         gtail
#             Genome for the tail gene
#         """
#         self.head: int = head
#         self.tail: int = tail
#         self.genome_head: Optional[int] = ghead
#         self.genome_tail: Optional[int] = gtail
#
#     def __str__(self):
#         return "h: {}, t: {}, gh: {}, gt: {}".format(self.head, self.tail, self.genome_head, self.genome_tail)
#
#     @staticmethod
#     def connect(path1: PGMPath,
#                 path2: PGMPath,
#                 path_l: PGMPath,
#                 which_genome: Optional[int] = None) -> Optional[PGMPath]:
#         """
#         Connect two PGMPaths
#
#         Parameters
#         ----------
#         path1
#             First path to connect
#         path2
#             Second path to connect
#         path_l
#             Path used for connecting two paths
#         which_genome
#             Which genome this path belongs to (Omit for GGH)
#
#         Returns
#         -------
#         Optional[PGMPath]
#             New PGMPath from the given paths or None if they are unable to be connected
#         """
#         head: int
#         tail: int
#         genome_head: int
#         genome_tail: int
#
#         if which_genome is not None:
#             if path1.head == path_l.head and path1.genome_head == path_l.genome_head and \
#                     path2.head == path_l.tail and path2.genome_head == path_l.genome_tail:
#                 head = path1.tail
#                 genome_head = path1.genome_tail
#                 tail = path2.tail
#                 genome_tail = path2.genome_tail
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.head
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.head
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.head == path_l.tail and path1.genome_head == path_l.genome_tail and \
#                     path2.head == path_l.head and path2.genome_head == path_l.genome_head:
#                 head = path1.tail
#                 genome_head = path1.genome_tail
#                 tail = path2.tail
#                 genome_tail = path2.genome_tail
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.head
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.head
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.head == path_l.tail and path1.genome_head == path_l.genome_tail and \
#                     path2.tail == path_l.head and path2.genome_tail == path_l.genome_head:
#                 head = path1.tail
#                 genome_head = path1.genome_tail
#                 tail = path2.head
#                 genome_tail = path2.genome_head
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.head
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.tail
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.head == path_l.head and path1.genome_head == path_l.genome_head and \
#                     path2.tail == path_l.tail and path2.genome_tail == path_l.genome_tail:
#                 head = path1.tail
#                 genome_head = path1.genome_tail
#                 tail = path2.head
#                 genome_tail = path2.genome_head
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.head
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.tail
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.tail == path_l.head and path1.genome_tail == path_l.genome_head and \
#                     path2.head == path_l.tail and path2.genome_head == path_l.genome_tail:
#                 head = path1.head
#                 genome_head = path1.genome_head
#                 tail = path2.tail
#                 genome_tail = path2.genome_tail
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.tail
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.head
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.tail == path_l.tail and path1.genome_tail == path_l.genome_tail and \
#                     path2.head == path_l.head and path2.genome_head == path_l.genome_head:
#                 head = path1.head
#                 genome_head = path1.genome_head
#                 tail = path2.tail
#                 genome_tail = path2.genome_tail
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.tail
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.head
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.tail == path_l.head and path1.genome_tail == path_l.genome_head and \
#                     path2.tail == path_l.tail and path2.genome_tail == path_l.genome_tail:
#                 head = path1.head
#                 genome_head = path1.genome_head
#                 tail = path2.head
#                 genome_tail = path2.genome_head
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.tail
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.tail
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#
#             elif path1.tail == path_l.tail and path1.genome_tail == path_l.genome_tail and \
#                     path2.tail == path_l.head and path2.genome_tail == path_l.genome_head:
#                 head = path1.head
#                 genome_head = path1.genome_head
#                 tail = path2.head
#                 genome_tail = path2.genome_head
#
#                 if head == 0 and genome_head == -1:
#                     head = path1.tail
#                     genome_head = which_genome
#
#                 if tail == 0 and genome_tail == -1:
#                     tail = path2.tail
#                     genome_tail = which_genome
#
#                 return PGMPath(head, tail, genome_head, genome_tail)
#         else:  # GGH version
#             if path1.head == path_l.head and path2.head == path_l.tail:
#                 return PGMPath(path1.tail, path2.tail)
#             elif path1.head == path_l.tail and path2.head == path_l.head:
#                 return PGMPath(path1.tail, path2.tail)
#             elif path1.head == path_l.tail and path2.tail == path_l.head:
#                 return PGMPath(path1.tail, path2.head)
#             elif path1.head == path_l.head and path2.tail == path_l.tail:
#                 return PGMPath(path1.tail, path2.head)
#             elif path1.tail == path_l.head and path2.head == path_l.tail:
#                 return PGMPath(path1.head, path2.tail)
#             elif path1.tail == path_l.tail and path2.head == path_l.head:
#                 return PGMPath(path1.head, path2.tail)
#             elif path1.tail == path_l.head and path2.tail == path_l.tail:
#                 return PGMPath(path1.head, path2.head)
#             elif path1.tail == path_l.tail and path2.tail == path_l.head:
#                 return PGMPath(path1.head, path2.head)
#
#         return None


def create_pgm_path(head: int, tail: int, genome_head: Optional[int] = None, genome_tail: Optional[int] = None) -> \
                                                                                                        Dict[str, int]:
    """ Dictionary representation of a PGMPath, originally represented as a class in Java. This change was made for
    optimization purposes

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
