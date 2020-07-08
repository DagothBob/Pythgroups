import yaml
from typing import List, Dict, Set


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
        A list of genes (rows) where each gene is represented as a list data (columns)
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

    return parsed_data


def parse_pgm_genome_data(gene_data: List[List[str]]) -> Dict[int, List[List[str]]]:
    """
    Group all the genes from the data set into their respective gene family IDs

    Parameters
    ----------
    gene_data : List[List[str]]
        A list of genes (rows) where each gene is represented as a list data (columns)

    Returns
    -------
    Dict[int, List[List[str]]]
        Key: gene family ID, value: list of each gene in the family
    """
    # Get each unique gene family ID from the data set
    all_family_ids: List[int] = [int(row[0]) for row in gene_data]
    min_id: int = min(all_family_ids)
    max_id: int = max(all_family_ids)
    family_ids: List[int] = [i for i in range(min_id, max_id + 1)]

    # Group all genes into their families with their IDs
    gene_families: Dict[int, List[List[str]]] = {}
    for i in range(min_id, max_id + 1):
        gene_families[i] = []
    for gene in gene_data:
        family_id = family_ids[int(gene[0]) - min_id]
        gene_families[family_id].append(gene)

    return gene_families
