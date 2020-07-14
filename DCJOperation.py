from enum import IntEnum

"""
Operation object for use in DCJRearrangement algorithm

Author: Holger Jensen
"""


class OperationTypes(IntEnum):
    NULL: int = 0
    INVERSION: int = 1
    TRANSLOCATION: int = 2
    FISSION: int = 3
    FUSION: int = 4


class TranslocationSubtypes(IntEnum):
    NULL: int = 0
    ADCB: int = 1  # (A,D,C,B)
    AnCnBD: int = 2  # (-A,-C,B,D)


class FusionSubtypes(IntEnum):
    NULL: int = 0
    TWO_REVERSED_PLUS_ONE: int = 1
    TWO_PLUS_ONE: int = 2
    ONE_PLUS_TWO: int = 3
    ONE_PLUS_TWO_REVERSED: int = 4


class DCJOperation:
    """
    Attributes
    ----------
    chromosome_1: int
        First chromosome to operate on
    chromosome_2: int
        Second chromosome to operate on
    node_1: int
        First node for operations (a)
    node_2: int
        Second node for operations (b)
    node_3: int
        Third node for operations (c)
    node_4: int
        Fourth node for operations (d)
    operation_type: int
        Inversion, Translocation, Fission, or Fusion
    operation_subtype: int
        Subtype for Translocation or Fusion
    """
    def __init__(self):
        """
        Constructor
        """
        self.chromosome_1: int = 0
        self.chromosome_2: int = 0
        self.node_1: int = 0
        self.node_2: int = 0
        self.node_3: int = 0
        self.node_4: int = 0
        self.operation_type: int = OperationTypes.NULL
        self.operation_subtype: int = TranslocationSubtypes.NULL

    def __str__(self) -> str:
        """
        String override

        Returns
        -------
        String representation of the object
        """
        return "Chromosome 1: " + str(self.chromosome_1) + "\n" + \
               "Chromosome 2: " + str(self.chromosome_2) + "\n" + \
               "Node 1: " + str(self.node_1) + "\n" + \
               "Node 2: " + str(self.node_2) + "\n" + \
               "Node 3: " + str(self.node_3) + "\n" + \
               "Node 4: " + str(self.node_4) + "\n" + \
               "Operation type: " + str(self.operation_type) + "\n" + \
               "Operation subtype: " + str(self.operation_subtype) + "\n"
