from typing import List, Tuple
import re
import networkx as nx
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import os
import statistics
import numpy as np

from wp1 import draw_graph, _search_edges
from constants import SEPCIES_MEDIUM_COMBINATIONS, AMINO_ACIDS, LILA
from custom_types import SpeciesMediumCombination

combination = ["blongum_adam","btheta_adam","ecoli_adam","eramosum_adam"]


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
    reached_AS = []
    reached_compounds = set(compounds_reverse.keys())

    for comp in compounds_reverse:

        if comp in AMINO_ACIDS:
            reached_AS.append(comp)

        reached_nodes = []
        for node in compounds_reverse[comp]:
            if node in H.nodes:
                reached_nodes.append(node)

        if not len(reached_nodes) == len(compounds_reverse[comp]):
            remove_nodes = remove_nodes + reached_nodes
            reached_compounds.remove(comp)

    H2 = H.subgraph(list(set(H.nodes)-set(remove_nodes)))

    return H2, reached_AS, reached_compounds


def rebuild_molecule_edges(G_full: nx.Graph, G_sub: nx.Graph) -> nx.Graph:
    '''
    Adds the molecule edges (NO_TRANSITION) again to the Graph G_sub
    '''

    edges = G_full.edges(data=True)
    no_trans_edges = [(u,v) for u,v,data in edges if data["transition"]=="TransitionType.NO_TRANSITION"]

    G_out = G_sub.copy()
    for u,v in no_trans_edges:
        if u in G_sub.nodes and v in G_sub.nodes:
                G_out.add_edge(u,v,transition="TransitionType.NO_TRANSITION")

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
    trans_types = []
    for i, edge in enumerate(G.edges):
        edges[i] = edge
        trans_type = edge_trans[edge] if edge in edge_trans else None
        if trans_type == "TransitionType.REACTION":
            color_map_edges[i] = "Green"
        elif trans_type == "TransitionType.NO_TRANSITION":
            color_map_edges[i] = "Black"
        else:
            color_map_edges[i] = "Yellow"

        trans_types.append(trans_type)

    nx.draw(G, pos=nx.layout.kamada_kawai_layout(G), node_color=color_map_nodes, with_labels=True, 
            nodelist=nodes, edgelist=edges, font_size=5, labels=labels, edge_color=color_map_edges, node_size=30)
    
    # Save the figure
    if output:
        plt.savefig(output)
    # Or show it
    else:
        plt.show()

    # Clear the figure
    plt.clf()

def plot_no_of_components(no_connected_comp: dict[SpeciesMediumCombination, int]):
    plt.bar(no_connected_comp.keys(), no_connected_comp.values(), color=LILA[3])
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
    plt.bar(x_axis - 0.3, max_values, 0.3, label="max", color=LILA[3])
    plt.bar(x_axis, avg_values, 0.3, label="average", color=LILA[2])
    plt.bar(x_axis + 0.3, median_values, 0.3, label="median", color=LILA[1])

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(0, 4001, step=500))
    plt.legend()

    plt.savefig("output/plots/component_sizes_bar_plot.png")

    plt.clf()

    # plot just average and median
    plt.bar(x_axis - 0.2, avg_values, 0.4, label="average", color='#7B37D9')
    plt.bar(x_axis + 0.2, median_values, 0.4, label="median", color='#AC25C7')

    plt.xlabel("species and media")
    plt.ylabel("component size")
    plt.xticks(x_axis, connected_comp_sizes)
    plt.yticks(np.arange(100, 1501, step=200))
    plt.legend()

    plt.savefig("output/plots/component_sizes_avg_median_bar_plot.png")


def plot_reached_aa(reached_aa: dict[SpeciesMediumCombination: Tuple[list[str],list[str]]]):

    number_aa = []
    number_aa_noCO2 = []
    df = pd.DataFrame([([0] * len(AMINO_ACIDS)) for i in range(len(combination))], columns = AMINO_ACIDS, index = combination)
    for c in combination:
        aa, aa_noCO2 = reached_aa[c]
        number_aa.append(len(aa))
        number_aa_noCO2.append(len(aa_noCO2))
        for a in AMINO_ACIDS:
            if a in reached_aa[c]:
                df[a][c] = 1

    # heatmap that shows which amino acids are synthesized per combination
    plt.title("reached amino acids")
    sns.heatmap(df, cbar=False, cmap="BuPu")
    # TODO: Plot beschriftung schöner machen
    plt.savefig("output/plots/species_medium_aa_heatmap_atn.png")

    # plot just average and median
    x_axis = np.arange(len(reached_aa.keys()))
    plt.bar(x_axis - 0.2, number_aa, 0.4, label="with CO2")
    plt.bar(x_axis + 0.2, number_aa_noCO2, 0.4, label="without CO2")

    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.xticks(x_axis, number_aa)
    plt.yticks(np.arange(100, 1501, step=200))
    plt.legend()

    plt.savefig("output/plots/reached_aa_compare_bar_plot.png")



def generate_atn_graph() -> Tuple[dict[SpeciesMediumCombination: int], dict[SpeciesMediumCombination: dict[str: float]],
                            dict[SpeciesMediumCombination: Tuple[list[str],list[str]]], dict[SpeciesMediumCombination: Tuple[list[str],list[str]]]]:

    no_connected_comp = {}
    connected_comp_sizes = {}
    reached_AS = {}
    reached_compounds = {}

    for species_medium_combination in SEPCIES_MEDIUM_COMBINATIONS:
        path = f"data/sihumix/{species_medium_combination}/{species_medium_combination}_cleaned.gml"

        # Load graph from gml file
        if not os.path.isfile(path):
            continue

        G = nx.read_gml(path)
        
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

        G_bfs_with_CO2, AS, compounds = bfs_from_molecule(G_without_trans, "D-glucose")
        G_bfs_without_CO2, AS_noCO2, compounds_noCO2 = bfs_from_molecule(G_without_CO2, "D-glucose")
        reached_AS[species_medium_combination] = (AS, AS_noCO2)
        reached_compounds[species_medium_combination] = (compounds, compounds_noCO2)

        #G_bfs_without_CO2_withMol = rebuild_molecule_edges(G, G_bfs_without_CO2)

        #draw_graph(G_bfs_without_CO2_withMol, f"output/plots/atn_graphs/{species_medium_combination}.png")

    return no_connected_comp, connected_comp_sizes, reached_AS, reached_compounds

if __name__ == "__main__":

    no_connected_comp, connected_comp_sizes, reached_AS, reached_compounds = generate_atn_graph()
    
    # Plot number of reached amino acids and conected components
    plot_reached_aa(reached_AS)

    # Plot number and size of connected components
    plot_no_of_components(no_connected_comp)
    plot_component_size(connected_comp_sizes)
    
    

    # Goals:
    # for every species medium combination
        # 1. Plot histogram of connected components by their size
        # 2. Calculate Density, Diameter
        # 3. Add them to the plot