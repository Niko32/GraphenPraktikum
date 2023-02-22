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
    compounds_reverse = {}
    for node in compounds:
        c = compounds[node]
        if c in compounds_reverse:
            compounds_reverse[c].append(node)
        else:
            compounds_reverse[c] = [node]

    start_nodes = compounds_reverse[molecule]

    H = _search_edges(G, start_nodes)

    remove_nodes = []
    for comp in compounds_reverse:
        reached_nodes = []
        for node in compounds_reverse[comp]:
            if node in H.nodes:
                reached_nodes.append(node)
        if not len(reached_nodes) == len(compounds_reverse[comp]):
            remove_nodes = remove_nodes + reached_nodes

    H2 = H.subgraph(list(set(H.nodes)-set(remove_nodes)))

    return H2


def rebuild_molecule_edges(G_full: nx.Graph, G_sub: nx.Graph) -> nx.Graph:
    '''
    Adds the molecule edges (NO_TRANSITION) again to the Graph G_sub
    '''
    add_edges = []
    for node in G_sub.nodes:
        edges_full = set(G_full.edges(node))
        edges_sub = set(G_sub.edges(node))
        for edge in edges_full - edges_sub:
            a, b = edge
            if a in G_sub.nodes and b in G_sub.nodes:
                add_edges.append((a, b))

    G_out = G_sub.copy()
    G_out.add_edges_from(add_edges)

    return G_out


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
    G_bfs_without_CO2 = bfs_from_molecule(G_without_CO2, "D-glucose")

    G_bfs_without_CO2_withMol = rebuild_molecule_edges(G, G_bfs_without_CO2)


    #third_largest_cc = G.subgraph(components[2])
    draw_graph(G_bfs_without_CO2_withMol)
    
    # Goals:
    # for every species medium combination
        # 1. Plot histogram of connected components by their size
        # 2. Calculate Density, Diameter
        # 3. Add them to the plot

    # Plot the avg. component sizes and other metrics for all combinations

    # Draw graph mti 3 kantenfarben

    # ATP oder so wegnehmen




