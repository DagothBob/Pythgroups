import sys
import yaml
from Bio import Phylo
from io import StringIO


"""
 Driver program for Pythgroups
 
 Command line usage: 
 $ python [directory of GenomeReconstruction.py] [algorithm]
 where [algorithm] is the desired genome reconstruction algorithm:
    - SmallPhylogeny
    - GenomeAliquoting
 
 Based on the small phylogeny program developed by C.Zheng & D.Sankoff (2011)

 Author: Oskar Jensen
"""

# TODO: setup.py?


def small_phylogeny() -> str:
    config_file = open("config.yaml", "r")
    config_data = yaml.safe_load(config_file)
    tree = Phylo.read(StringIO(config_data.get("tree_structure")), "newick")
    genome_file = open(config_data.get("genome_file"), "r")

    return "--Tree structure--\n" + str(tree) + "\n\n--Genome file--\n" + genome_file.read()


def genome_aliquoting() -> str:
    return "Hello yes this is genome aliquoting"


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
        "GenomeAliquoting": genome_aliquoting()
    }.get(alg, "Algorithm doesn't exist")


def main():
    algorithm = sys.argv[1]  # The first argument when calling the program
    print(get_algorithm(algorithm))


if __name__ == '__main__':
    main()
