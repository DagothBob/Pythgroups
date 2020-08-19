import networkx as nx

from typing import List, Optional, Dict, Set
from Bio.Phylo.Newick import Tree
import matplotlib.pyplot as plt


class NetworkxNode:
    """
    Represents a node in the networkx graph, containing additional information like a unique integer identifier

    Attributes
    ----------
    name : str
        Name of the genome
    neighbors : List[int]
        The genome_ids of this genome's neighbors
    genome_id : int
        Integer identifier for this genome
    """
    def __init__(self, name: str, neighbors: List[int], genome_id: int):
        """
        Constructor

        Attributes
        ----------
        name : str
            Name of the genome
        neighbors : List[int]
            The genome_ids of this genome's neighbors
        genome_id : int
            Integer identifier for this genome
        """
        self.name: str = name
        self.neighbors: List[int] = neighbors
        self.genome_id: int = genome_id

    def __str__(self):
        return "name: {}, neighbors: {}, genome_id: {}".format(self.name, self.neighbors, self.genome_id)


def parse_medians(graph_nodes: List[NetworkxNode]) -> List[NetworkxNode]:
    """
    Get the median nodes from a list of GenomeNodes

    Parameters
    ----------
    graph_nodes : List[NetworkxNode]
        List of GraphNodes to parse the medians from

    Returns
    -------
    Dict[Clade, List[Clade]]
        Key: median to reconstruct, value: list of 3 neighboring clades
    """
    medians: List[NetworkxNode] = list()
    for node in graph_nodes:
        if (len(node.neighbors)) == 3:
            medians.append(NetworkxNode(node.name, node.neighbors, node.genome_id))
        # TODO: Exception handling for illegal tree structures (eg. node with > 3 neighbors)

    return medians


def genome_nodes_from_tree(tree: Tree) -> List[NetworkxNode]:
    """
    Converts a tree into a networkx graph consisting of GenomeNode representations of each node

    Parameters
    ----------
    tree : Tree
        Tree object to convert into a graph

    Returns
    -------
    List[NetworkxNode]
        A list of GenomeNodes representing each node, their neighbors, and their unique IDs
    """
    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()
    graph_indexes = dict()
    node_index = 0
    for t in terminals:
        graph_indexes[t.name] = node_index
        node_index += 1
    for nt in nonterminals:
        if nt.name is not None and len(nt.clades) > 0:
            graph_indexes[nt.name] = node_index
            node_index += 1

    graph: nx.Graph = nx.Graph()
    for nt in nonterminals:
        if nt.name is not None and len(nt.clades) > 0:
            for clade in nt.clades:
                graph.add_edge(nt.name, clade.name)

    graph_nodes = list()
    for node in graph:
        genome_id = graph_indexes[node]
        neighbor_ids = list()
        for neighbor in graph.neighbors(node):
            neighbor_ids.append(graph_indexes[neighbor])
        graph_nodes.append(NetworkxNode(node, neighbor_ids, genome_id))

    return graph_nodes


def get_node(nodes: List[NetworkxNode], target: int) -> Optional[NetworkxNode]:
    """ Get the node with the specified genome_id in a list of NetworkxNodes

    Parameters
    ----------
    nodes : List[NetworkxNode]
        List of nodes to search in
    target : int
        The genome_id of the target node

    Returns
    -------
    NetworkxNode
        Target node
    """
    for node in nodes:
        if node.genome_id == target:
            return node

    return None


def print_tree_structure(nodes: List[NetworkxNode], distances: Dict):
    graph: nx.Graph = nx.Graph()

    edges = set()
    color_map: List[str] = []
    for node in nodes:
        if len(node.neighbors) == 3:
            color_map.append("#ddaaaa")
        else:
            color_map.append("#bbbbbb")
        for neighbor in node.neighbors:
            edges.add((get_node(nodes, node.genome_id), get_node(nodes, neighbor)))
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)

    pos = nx.spring_layout(graph)
    plt.figure()
    nx.draw(graph, pos=pos, labels={node: node.name for node in graph.nodes()},
            node_size=3000, node_color=color_map, edge_color="black", alpha=1, width=1)
    nx.draw_networkx_edge_labels(graph, pos=pos, edge_labels=distances)
    plt.title("Genome relations")
    plt.margins(x=0.1, y=0.1)
    plt.show()
