import sys
import yaml

from BPGDistance import BPGDistance
from DCJOperations import DCJOperations, OperationTypes
from GenomeInString import GenomeInString

from Bio import Phylo
from io import StringIO
from typing import List, Dict, TextIO
from enum import Enum


"""
 Driver program for Pythgroups
 
 Command line usage: 
 $ python [directory of GenomeReconstruction.py] [algorithm]
 where [algorithm] is the desired genome reconstruction algorithm:
    - SmallPhylogeny
    - GenomeAliquoting
    - DCJRearrangements
 
 Based on the small phylogeny program developed by C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""


def small_phylogeny() -> str:
    """
    Reconstructs the ancestor(s) of known modern genomes given an unrooted binary phylogenetic tree
    using the Pathgroups algorithm.

    Returns
    -------

    """
    config_file = open("config.yaml", "r")
    config_data = yaml.safe_load(config_file)
    tree = Phylo.read(StringIO(config_data.get("tree_structure")), "newick")
    genome_file = open(config_data.get("genome_file"), "r")

    return "--Tree structure--\n" + str(tree) + "\n\n--Genome file--\n" + genome_file.read()


def genome_aliquoting() -> str:
    """
    See the 2010 paper, section 2.5

    Returns
    -------

    """

    return "Hello yes this is genome aliquoting"


def parse_genomes(genome_file: TextIO) -> Dict[str, List[str]]:
    """
    Parses the genome data from a file into a dictionary containing each genome and their chromosomes

    Parameters
    ----------
    genome_file : TextIO
        text file containing all the genomes and their chromosomes

    Returns
    -------
    Dict[str, List[str]]
        A dictionary containing each genome's chromosomes,
        where keys are genome headers and values are lists of chromosomes
    """
    genomes: Dict[str, List[str]] = {}
    line: str = genome_file.readline()
    header: str = ""
    while len(line) != 0:

        # ">" indicates a header, sets that as the current genome to add chromosomes to
        if line.startswith(">"):
            header = line.replace(">", "").replace("\n", "")

        # Add chromosomes to the current genome until an empty line, another header, or the end of file is found
        else:
            chromosomes: List[str] = []
            while len(line) != 0 and not line.startswith(">") and line != "\n":
                chromosome: str = line.replace("$", "").replace("\n", "")
                chromosomes.append(chromosome)
                line = genome_file.readline()
            genomes[header] = chromosomes
        line = genome_file.readline()

    return genomes


def dcj_rearrangements():
    """
    Performs double cut-and-join (DCJ) operations until the source genome matches the target genome,
    records the state of the intermediate genomes and their distances to the target genome at each step

    Performs DCJ operations on the first 2 genomes found in the genome file
    """
    # Create the dictionary of genomes from the input file
    config_file = open("config.yaml", "r")
    config_data = yaml.safe_load(config_file)
    genome_file = open(config_data.get("genome_file"), "r")
    genomes: Dict[str, List[str]] = parse_genomes(genome_file)

    # Get the first 2 genomes from the input file
    values_view = genomes.values()
    value_iterator = iter(values_view)
    genome1: List[str] = next(value_iterator)
    genome2: List[str] = next(value_iterator)

    # Calculate the initial BPG distance
    bpg_dist: BPGDistance = BPGDistance(genome1, genome2)
    bpg_dist.calculate_distance()
    cur_dist: int = bpg_dist.distance

    # Perform DCJ operations until distance == 0 or there are no more DCJ operations to perform (?)
    operation_types: List[int] = [OperationTypes.INVERSION, OperationTypes.TRANSLOCATION]
    dcj: DCJOperations = DCJOperations(genome1, genome2)
    dcj.initial_value()
    more: bool = True
    while cur_dist > 0 and more:
        more = False
        rearrange_state: List[GenomeInString] = dcj.get_result(1, 30, -2, operation_types, 10)
        if len(rearrange_state) > 0:
            more = True
            for genome in rearrange_state:
                print("*******")
                for i in range(0, len(genome.chromosomes)):
                    print("Chromosome " + str(i) + "\n" + genome.chromosomes[i])

            new_genome: GenomeInString = rearrange_state[len(rearrange_state) - 1]
            bpg_dist = BPGDistance(new_genome.chromosomes, genome2)
            bpg_dist.calculate_distance()
            cur_dist = bpg_dist.distance
            print("a run, steps: " + str(len(rearrange_state)) + ", cur_dist: " + str(cur_dist))
            dcj = DCJOperations(new_genome.chromosomes, genome2)
            dcj.initial_value()


def get_algorithm(alg: str) -> str:
    """
    Calls the appropriate function based on the user input.

    Parameters
    ----------
    alg : str
        User input, the name of the desired algorithm

    Returns
    -------
    function
        The appropriate function based on the user input, or an error if the algorithm is invalid
    """
    # Essentially a switch statement
    return {
        "SmallPhylogeny": small_phylogeny(),
        "GenomeAliquoting": genome_aliquoting(),
        "DCJRearrangements": dcj_rearrangements()
    }.get(alg, "Algorithm doesn't exist")


def main():
    # algorithm = sys.argv[1]  # The first argument when calling the program
    output = get_algorithm("DCJRearrangements")


if __name__ == '__main__':
    main()
