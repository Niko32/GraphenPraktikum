from typing import List, Tuple
import re
import networkx as nx
from matplotlib import pyplot as plt

from wp1 import draw_graph, _search_edges


def remove_no_transitions(G: nx.Graph) -> nx.Graph:
    '''
    Removes  the NO_TRANSITION edges from the network
    '''
    H = nx.Graph()
    H = G.copy()

    transition = nx.get_edge_attributes(H, "transition")
    for edge in H.edges:
        if transition[edge] == "TransitionType.NO_TRANSITION":
            H.remove_edge(*edge)
    return H


def remove_CO2(G: nx.Graph) -> nx.Graph:
    '''
    Removes all edges from or to CO2 from the network
    '''
    H = nx.Graph()
    H = G.copy()

    compounds = nx.get_node_attributes(H, "compound_name")
    for edge in H.edges:
        source, target = edge
        if compounds[source] == "CO2" or compounds[target] == "CO2":
            H.remove_edge(*edge)
    return H


def bfs_from_molecule(G: nx.Graph, molecule: str) -> nx.Graph:
    '''
    breadth first search starting at the given molecule 
    '''
    compounds = nx.get_node_attributes(G, "compound_name")
    start_nodes = [node for node, c in compounds.items() if c == molecule]

    H = _search_edges(G, start_nodes)

    return H


if __name__ == "__main__":
    G = nx.read_gml("data/sihumix/ecoli_adam/ecoli_adam_cleaned.gml")
    
    G_without_trans = remove_no_transitions(G)
    print(nx.number_connected_components(G_without_trans))
    components = sorted(nx.connected_components(G_without_trans), key=len, reverse=True)
    print(nx.density(G))
    G_without_CO2 = remove_CO2(G_without_trans)
    print(len(G.edges))
    print(len(G_without_trans.edges))
    print(len(G_without_CO2.edges))

    G_bfs_with_CO2 = bfs_from_molecule(G_without_trans, "D-glucose")


    #third_largest_cc = G.subgraph(components[2])
    draw_graph(G_bfs_with_CO2)
    



