from typing import List, Tuple
import re
import networkx as nx
from matplotlib import pyplot as plt

from wp1 import draw_graph


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


if __name__ == "__main__":
    G = nx.read_gml("data/sihumix/ecoli_adam/ecoli_adam_cleaned.gml")
    
    G_without_trans = remove_no_transitions(G)
    print(nx.number_connected_components(G_without_trans))
    components = sorted(nx.connected_components(G_without_trans), key=len, reverse=True)
    print(nx.density(G))
    third_largest_cc = G.subgraph(components[2])
    draw_graph(third_largest_cc)
    
    # Goals:
    # 1. Plot histogram of connected components by their size
    # 2. Calculate Density, Diameter
    # 3. Add them to the plot

    # Draw graph mti 3 kantenfarben




