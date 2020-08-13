from typing import List, Dict

import yaml

CONFIG_DIR = "config.yaml"
CONFIG_GENOME_FILE = "genome_file"
CONFIG_TREE_STRUCTURE = "tree_structure"


def parse_gene_file(config_dir: str) -> List[List[str]]:
    """
    Parse the raw gene data, extracting the individual data points for each gene

    Parameters
    ----------
    config_dir
        Directory of the config file

    Returns
    -------
    List[List[str]]
        2D list of genes each represented as lists of data.
        Rows: Genes
        Columns: familyID, chr, start, end, strand, genome
    """
    parsed_data: List[List[str]] = []

    config_file = open(config_dir, "r")
    config_data = yaml.safe_load(config_file)

    # Parse the genes, where genes are delimited by new lines and data points by tab characters
    with open(config_data.get(CONFIG_GENOME_FILE)) as genome_file:
        for line in genome_file:
            data: List[str] = line.strip().split("\t")
            if data[0] != "geneName":  # Ignore the header text
                parsed_data.append([data[i] for i in [1, 2, 3, 4, 5, 8]])  # familyID, chr, start, end, strand, genome

    config_file.close()

    return parsed_data


def group_genomes(gene_data: List[List[str]]) -> Dict[str, List[List[str]]]:
    """
    Group all the genes from the parsed data set into their respective genome IDs.
    Formats the data in the format needed for pathgroups

    Parameters
    ----------
    gene_data : List[List[str]]
        A list of all parsed genes in the data set

    Returns
    -------
    Dict[str, List[List[str]]]
        Key: genome ID, value: 2D list of each gene in the genome.
        Rows: Genes
        Columns: strand + family ID, chromosome, start
    """

    # Get each unique genome from the data set
    all_genome_names: List[str] = [row[5] for row in gene_data]
    genome_names: List[str] = list(set(all_genome_names))

    # Group each gene into their respective genomes, formatting the genes which are represented as their family IDs
    genomes: Dict[str, List[List[str]]] = {}
    for gene in gene_data:
        genome_name: str = next(x for x in genome_names if x == gene[5])
        if genome_name not in genomes:
            genomes[genome_name] = []
        if gene[4] == "+":
            genomes[genome_name].append(["{}".format(gene[0]), gene[1], gene[2]])
        else:
            genomes[genome_name].append(["-{}".format(gene[0]), gene[1], gene[2]])

    return genomes


def to_pathgroups_format(genome_data: Dict[str, List[List[str]]]) -> str:
    """
    Format the genome data needed for the pathgroups algorithm.
    Sorts all the genes in each genome based on their start positions and groups them into their respective chromosomes.

    Parameters
    ----------
    genome_data : Dict[str, List[List[str]]]
        Genes grouped into their respective genomes

    Returns
    -------
    str
        Genome data in the format required for the pathgroups algorithm
    """
    output: str = ""
    for genome, genes in genome_data.items():
        output += ">{}\n".format(genome)
        genes.sort(key=lambda x: int(x[2]))

        all_chromosomes: List[int] = [int(gene[1]) for gene in genes]
        chromosomes: List[int] = list(set(all_chromosomes))
        chromosomes.sort()
        for chromosome in chromosomes:
            chromosome_str = ["{} ".format(gene[0]) for gene in genes if int(gene[1]) == chromosome]
            for gene in chromosome_str:
                output += gene
            output += "$\n"
        output += "\n"
    return output
