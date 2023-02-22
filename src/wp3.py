from typing import List, Tuple
import re
import networkx as nx
from matplotlib import pyplot as plt
import os
import statistics
import numpy as np

from wp1 import draw_graph, _search_edges
from constants import SEPCIES_MEDIUM_COMBINATIONS
from custom_types import SpeciesMediumCombination


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

def draw_graph(G: nx.Graph, output = ""):

    print("Drawing Graph...")

    node_element = nx.get_node_attributes(G, "element")
    node_compounds = nx.get_node_attributes(G, "compound_name")
    nodes, labels, color_map_nodes = len(G)*[None], {}, len(G)*[None]
    for i, node in enumerate(G.nodes):
        nodes[i] = node
        labels[node] = node_element[node]
        if node_compounds[node] == "D-Glucose":
            color_map_nodes[i] = "Red"
        else:
            color_map_nodes[i] = "White"

    edge_trans = nx.get_edge_attributes(G, "transition")
    print(edge_trans)
    edges, color_map_edges = len(G.edges)*[None], len(G.edges)*[None]
    for i, (u, v) in enumerate(G.edges):
        edges[i] = (u,v)
        trans_type = edge_trans[(u,v)]
        if trans_type == "TransitionType.REACTION":
            color_map_edges[i] = "Green"
        elif trans_type == "TransitionType.NO_TRANSITION":
            color_map_edges[i] = "Black"
        else:
            color_map_edges[i] = "Yellow"

    nx.draw(G, pos=nx.layout.kamada_kawai_layout(G), node_color=color_map_nodes, with_labels=True, 
            nodelist=nodes, edgelist=edges, font_size=5, labels=labels, edge_color=color_map_edges)
    
    # Save the figure
    if output:
        plt.savefig(output)
    # Or show it
    else:
        plt.show()

    # Clear the figure
    plt.clf()

def plot_no_of_components(no_connected_comp: dict[SpeciesMediumCombination, int]):
    plt.bar(no_connected_comp.keys(), no_connected_comp.values())
    plt.xlabel("species and media")
    plt.ylabel("number of connected components")
    plt.yticks(np.arange(0, 6, step=1))

    plt.savefig("output/plots/no_of_components_bar_plot.png")

def plot_component_size(connected_comp_sizes: dict[SpeciesMediumCombination, dict]):
    #min_values = [sizes["min"] for sizes in connected_comp_sizes.values()]
    max_values = [sizes["max"] for sizes in connected_comp_sizes.values()]
    avg_values = [sizes["avg"] for sizes in connected_comp_sizes.values()]
    median_values = [sizes["median"] for sizes in connected_comp_sizes.values()]

    x_axis = np.arange(len(connected_comp_sizes.keys()))
    #plt.bar(x_axis - 0.4, min_values, 0.2, label="min")
    plt.bar(x_axis - 0.3, max_values, 0.3, label="max")
    plt.bar(x_axis, avg_values, 0.3, label="average")
    plt.bar(x_axis + 0.3, median_values, 0.3, label="median")

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(0, 4001, step=500))
    plt.legend()

    plt.savefig("output/plots/component_sizes_bar_plot.png")

    plt.close()

    # plot just average and median
    plt.bar(x_axis - 0.2, avg_values, 0.4, label="average")
    plt.bar(x_axis + 0.2, median_values, 0.4, label="median")

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(100, 1501, step=200))
    plt.legend()

    plt.savefig("output/plots/component_sizes_avg_median_bar_plot.png")

def draw_graph(G: nx.Graph, output = ""):

    print("Drawing Graph...")

    node_element = nx.get_node_attributes(G, "element")
    node_compounds = nx.get_node_attributes(G, "compound_name")
    nodes, labels, color_map_nodes = len(G)*[None], {}, len(G)*[None]
    for i, node in enumerate(G.nodes):
        nodes[i] = node
        labels[node] = node_element[node]
        if node_compounds[node] == "D-Glucose":
            color_map_nodes[i] = "Red"
        else:
            color_map_nodes[i] = "White"

    edge_trans = nx.get_edge_attributes(G, "transition")
    print(edge_trans)
    edges, color_map_edges = len(G.edges)*[None], len(G.edges)*[None]
    for i, (u, v) in enumerate(G.edges):
        edges[i] = (u,v)
        trans_type = edge_trans[(u,v)]
        if trans_type == "TransitionType.REACTION":
            color_map_edges[i] = "Green"
        elif trans_type == "TransitionType.NO_TRANSITION":
            color_map_edges[i] = "Black"
        else:
            color_map_edges[i] = "Yellow"

    nx.draw(G, pos=nx.layout.kamada_kawai_layout(G), node_color=color_map_nodes, with_labels=True, 
            nodelist=nodes, edgelist=edges, font_size=5, labels=labels, edge_color=color_map_edges)
    
    # Save the figure
    if output:
        plt.savefig(output)
    # Or show it
    else:
        plt.show()

    # Clear the figure
    plt.clf()

def plot_no_of_components(no_connected_comp: dict[SpeciesMediumCombination, int]):
    plt.bar(no_connected_comp.keys(), no_connected_comp.values())
    plt.xlabel("species and media")
    plt.ylabel("number of connected components")
    plt.yticks(np.arange(0, 6, step=1))

    plt.savefig("output/plots/no_of_components_bar_plot.png")

def plot_component_size(connected_comp_sizes: dict[SpeciesMediumCombination, dict]):
    #min_values = [sizes["min"] for sizes in connected_comp_sizes.values()]
    max_values = [sizes["max"] for sizes in connected_comp_sizes.values()]
    avg_values = [sizes["avg"] for sizes in connected_comp_sizes.values()]
    median_values = [sizes["median"] for sizes in connected_comp_sizes.values()]

    x_axis = np.arange(len(connected_comp_sizes.keys()))
    #plt.bar(x_axis - 0.4, min_values, 0.2, label="min")
    plt.bar(x_axis - 0.3, max_values, 0.3, label="max")
    plt.bar(x_axis, avg_values, 0.3, label="average")
    plt.bar(x_axis + 0.3, median_values, 0.3, label="median")

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(0, 4001, step=500))
    plt.legend()

    plt.savefig("output/plots/component_sizes_bar_plot.png")

    plt.clf()

    # plot just average and median
    plt.bar(x_axis - 0.2, avg_values, 0.4, label="average")
    plt.bar(x_axis + 0.2, median_values, 0.4, label="median")

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(100, 1501, step=200))
    plt.legend()

    plt.savefig("output/plots/component_sizes_avg_median_bar_plot.png")

if __name__ == "__main__":

    no_connected_comp = {}
    connected_comp_sizes = {}

    for species_medium_combination in SEPCIES_MEDIUM_COMBINATIONS:
        path = f"data/sihumix/{species_medium_combination}/{species_medium_combination}_cleaned.gml"

        # Load graph from gml file
        if not os.path.isfile(path):
            continue

        G = nx.read_gml(path)

        remove_no_transitions(G)

        no_connected_comp[species_medium_combination] = nx.number_connected_components(G)
        component_sizes = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
        connected_comp_sizes[species_medium_combination] = {
            "min": min(component_sizes), 
            "max": max(component_sizes),
            "avg": sum(component_sizes)/len(component_sizes),
            "median": statistics.median(component_sizes)}
    
        G_without_trans = remove_no_transitions(G)
        print(nx.density(G))
        G_without_CO2 = remove_CO2(G_without_trans)
        print(len(G.edges))
        print(len(G_without_trans.edges))
        print(len(G_without_CO2.edges))

        G_bfs_with_CO2 = bfs_from_molecule(G_without_trans, "D-glucose")
        G_bfs_without_CO2 = bfs_from_molecule(G_without_CO2, "D-glucose")

        G_bfs_without_CO2_withMol = rebuild_molecule_edges(G, G_bfs_without_CO2)

    # Plot number and size of connected components
    plot_no_of_components(no_connected_comp)
    plot_component_size(connected_comp_sizes)


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




