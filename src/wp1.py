from typing import List, Tuple
import networkx as nx
from networkx.classes.reportviews import OutEdgeDataView
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image
import os
import yaml

from custom_types import Reaction
from constants import AMINO_ACIDS, COFACTORS, SEPCIES_MEDIUM_COMBINATIONS


def seperate_blocks(file_path: str) -> List[List[str]]:
    """ 
    Takes in the .smiles_list file and outputs the reaction blocks 
    as a list of lists of four strings
    """

    extra_block = [
        "Bigg ID: R_LEUTA MetaNetXId: MNXR192685 Reversible: True",
        "ECs: 2.6.1.42;2.6.1.6;2.6.1.67",
        "L-leucine + 2-oxoglutarate = 4-methyl-2-oxopentanoate + L-glutamate",
        "CC(C)C[CH](C(=O)O)N.O=C([O-])CCC(=O)C(=O)[O-]>>CC(C)CC(=O)C(=O)[O-].[NH3+][CH](CCC(=O)[O-])C(=O)[O-]"
    ]

    reaction_blocks = [extra_block]
    reaction_block = []
    skip = False

    with open(file_path, 'r') as f:
        for l in f:
            if l == "\n":
                skip = False
                if not reaction_block == []:
                    reaction_blocks.append(reaction_block)
                    reaction_block = []
            elif "Growth" in l:
                skip = True
            else:
                if not skip:
                    reaction_block.append(l)
    reaction_blocks.append(reaction_block)
    return reaction_blocks
    

def extract_compounds(reaction_block: List[str]) -> Reaction:
    """
    Takes a block of four string representing one reaction and parses it into
    our reaction class
    """
    reac: Reaction = {}

    line1 = reaction_block[0].split()
    reac["bigg_id"] = line1[2]
    reac["metanetx_id"] = line1[4]
    reac["reversible"] = (line1[6] == "True")

    line3_educts = reaction_block[2].split('=')[0].strip()
    line3_products = reaction_block[2].split('=')[1].strip()
    reac["educts"] = line3_educts.split(" + ")
    reac["products"] = line3_products.split(" + ")

    line4_educts = reaction_block[3].split('>>')[0].strip()
    line4_products = reaction_block[3].split('>>')[1].strip()
    reac["smiles_educts"] = line4_educts.split(".")
    reac["smiles_products"] = line4_products.split(".")

    return reac


def construct_graph(reactions: List[Reaction]) -> nx.DiGraph:
    """ Takes a list of reactions to construct a network x graph from it """
    G = nx.DiGraph()
    for r in reactions:
        G.add_node(r["bigg_id"], reaction=True)
        for e, smiles_e in zip(r["educts"], r["smiles_educts"]):
            weight = r["educts"].count(e)
            G.add_node(e, smiles=smiles_e)
            G.add_edge(e, r["bigg_id"], weight=weight)

            if r["reversible"]:
                G.add_node(r["bigg_id"] + "_rev", reaction=True)
                G.add_edge(r["bigg_id"] + "_rev", e, weight=weight)

        for p, smiles_p in zip(r["products"], r["smiles_products"]):
            weight = r["products"].count(p)
            G.add_node(p, smiles=smiles_p)
            G.add_edge(r["bigg_id"], p, weight=weight)

            if r["reversible"]:
                G.add_node(r["bigg_id"] + "_rev", reaction=True)
                G.add_edge(p, r["bigg_id"] + "_rev", weight=weight)

    return G

def _visit_node(G: nx.DiGraph, v: str) -> Tuple[nx.DiGraph, OutEdgeDataView]:
    """ 
    Visits a node if it has not been visited yet and returns its outgoing edges 
    """

    if not G.nodes[v].get("visited"):
        print(f"Marking node {v} as visited")
        G.nodes[v]["visited"] = True
        return G, G.edges(v, data="visited")
    else:
        return G, []
    
def _inputs_complete(G: nx.DiGraph, v: str):
    """ Check if a node meets all requirements to be visited """

    inputs_complete = True

    # Check for educts if the node is a reaction
    if G.nodes[v].get("reaction"):
        for u, v in G.in_edges(v):
            if not G.nodes[u].get("visited"):
                inputs_complete = False
                break

    return inputs_complete

def _build_visited_subgraph(G: nx.DiGraph) -> nx.DiGraph:
    """ Builds a subgraph from all the nodes that have been visited """

    # Include all visited nodes
    visited_nodes = set()
    for node, visited in G.nodes(data="visited"):
        if visited:
            visited_nodes.add(node)
    H: nx.DiGraph = G.subgraph(visited_nodes)

    # Remove molecules that are not used in any reaction
    connected_nodes = set()
    for u, v in H.edges():
        connected_nodes.add(u)
        connected_nodes.add(v)
    return H.subgraph(connected_nodes)

def _create_gif(image_paths: List[str]):
    """ Creates a gif from a list of png image paths """
    print("Creating traversal gif...")
    images = [Image.open(path) for path in image_paths if path.endswith(".png")]
    images[0].save("output/plots/traversal.gif", save_all = True, 
                   append_images = images[1:], optimize = False, 
                   duration = len(image_paths)*10)
    
def _search_edges(G: nx.DiGraph, start_nodes: List[str], reverse = False, verbose = False, n: int = None):
    """
    Takes a graph and a set of starting nodes to perform bf search and rerturn the resulting subgraph 
    """
    
    outgoing_edges: List[Tuple[str, str, dict]] = []
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    # Unvisit all nodes
    nx.set_node_attributes(G, False, "visited")

    # Assure the starting nodes are in G
    start_nodes = set(start_nodes).intersection(G.nodes)

    # Mark the starting nodes as visited
    for node in start_nodes:
        G.nodes[node]["visited"] = True
        outgoing_edges.extend(G.edges(node, data="visited"))

    iteration = 0

    # Loop as long as there is a non empty queue
    while len(outgoing_edges) > 0:
        print(f"### Iteration {iteration} ###")

        # Save the graph figure
        if verbose:
            H = _build_visited_subgraph(G)
            traversal_output_dir = f"output/plots/bf_traversal"
            output_path = f"{traversal_output_dir}/{iteration:02}.png"

            # Use the layout as produced by the last graph
            with open("output/plots/bf_traversal/positions.yml", "r") as f:
                positions = yaml.safe_load(f)
            draw_graph(H, output_path, pos=positions)

        # Visit all edges in the current queue
        for u, v, visited in outgoing_edges:

            print(f"Traversing edge ({u},{v})")

            # Skip visiting the node if it does not have all required inputs
            if not reverse and not _inputs_complete(G, v):
                continue

            # Visit the node
            G, edges = _visit_node(G, v)
            new_outgoing_edges.extend(edges)

        # Continue with the next iteration
        iteration += 1
        outgoing_edges = new_outgoing_edges
        new_outgoing_edges = []

        # Stop after number of iterations is reached
        if n and n <= iteration:
            break

    H = _build_visited_subgraph(G)

    if verbose:
        # Compute node postitions to save them
        print("Computing final graph layout...")
        positions = nx.layout.kamada_kawai_layout(H)
        positions = {key: value.tolist() for key, value in positions.items()}
        with open("output/plots/bf_traversal/positions.yml", "w") as f:
            yaml.safe_dump(positions, f)

        # Save the last graph figure
        draw_graph(H, output_path, pos=positions)
        _create_gif([f"{traversal_output_dir}/{file}" for file in os.listdir(traversal_output_dir)])

    return H

def bf_traversal(G: nx.DiGraph, metabolites: List[str] = [], verbose = False, n: int = None, use_cofactors = True) -> nx.DiGraph:
    """ 
    Takes a set of metabolites to form a subgraph constructed from them
    """  

    print("##### bf_traversal #####")

    start_nodes = COFACTORS + metabolites if use_cofactors else metabolites  
    H = _search_edges(G, start_nodes, verbose=verbose, n=n)

    # Print out amino acids that have been reached
    reached_acids = set(H.nodes.keys()).intersection(AMINO_ACIDS)
    print(len(reached_acids), "/", len(AMINO_ACIDS), "acids have been reached")

    return H

def reverse_bf_traversal(G: nx.DiGraph) -> List[nx.DiGraph]:

    print("##### reverse_bf_traversal #####")
    
    # Reverse the edges and set new start nodes
    G = G.reverse()
    nx.set_node_attributes(G, False, "visited")

    # Exclude acids that are not in the graph
    acids = set(AMINO_ACIDS).intersection(G.nodes)

    subgraphs = []

    # Create one subgraph for every amino acid
    for acid in acids:

        print(f"Creating subgraph for {acid}")
        H = _search_edges(G, [acid], reverse=True)

        subgraphs.append(H)
            
    return subgraphs

def intersect_subgraph(s: nx.DiGraph, subgraphs: List[nx.DiGraph]) -> nx.DiGraph:
    """ Intersection of the given subgraphs """
    A = nx.DiGraph()
    for sub in subgraphs:
        sub = sub.reverse()
        A.add_nodes_from(sub.nodes(data=True))
        A.add_edges_from(sub.edges(data=True))

    G = A.copy()
    G.remove_nodes_from(n for n in A if n not in s)
    G.remove_edges_from(e for e in A.edges if e not in s.edges)
    return G

def draw_graph(G: nx.DiGraph, output = "", show_reactions = False, pos=None):

    print("Drawing Graph...")

    # Configure display settings based on the reaction flag
    nodes, sizes, labels, color_map = len(G)*[None], np.empty(len(G), dtype=int), {}, len(G)*[None]
    for i, (node, reaction_flag) in enumerate(G.nodes(data="reaction")):
        nodes[i] = node
        sizes[i] = 50 - 49 * (reaction_flag or 0)
        if show_reactions:
            labels[node] = node
        else:
            labels[node] = "" if reaction_flag else node
        if node in AMINO_ACIDS:
            color_map[i] = "Red"
        elif node in ["D-glucose", "biomass"]:
            color_map[i] = "Yellow"
        elif reaction_flag:
            color_map[i] = "Blue"
        else:
            color_map[i] = "Green"

    weights = len(G.edges)*[None]
    for i, (u, v) in enumerate(G.edges):
        weight = G[u][v].get("weight")
        weights[i] = weight if weight else 1
        weights[i] *= 0.25

    if not pos:
        pos = nx.layout.kamada_kawai_layout(G)

    nx.draw(G, pos=pos, node_color=color_map, with_labels=True, 
            nodelist=nodes, node_size=sizes, font_size=5, arrowsize=4, labels=labels, width=weights)
    
    # Save the figure
    if output:
        plt.savefig(output)
    # Or show it
    else:
        plt.show()

    # Clear the figure
    plt.clf()

def build_graph(file_path: str) -> nx.DiGraph:
    """ Builds a graph from a path to a smiles list file """
    reaction_block_list = seperate_blocks(file_path)
    reactions = [extract_compounds(rb) for rb in reaction_block_list]
    return construct_graph(reactions)

def build_subgraph(file_path: str, verbose = False) -> nx.DiGraph:
    """ 
    Builds a subgraph containing the paths from glucose to all amino acids 
    for one combination of molecule and medium 
    """
    G = build_graph(file_path)

    S = bf_traversal(G, ["D-glucose"], verbose)
    A = reverse_bf_traversal(S)

    subgraph = intersect_subgraph(G,A)

    # Remove COFACTORS
    # for cofactor in COFACTORS:
    #     if cofactor in subgraph:
    #         subgraph.remove_node(cofactor)

    return subgraph

if __name__ == "__main__":
    # file_paths = [f"data/sihumix/{c}/{c}.smiles_list" for c in SEPCIES_MEDIUM_COMBINATIONS]

    # subgraphs: List[nx.DiGraph] = []
    # for i, path in enumerate(file_paths):

    #     save_path = f"output/subgraphs/{SEPCIES_MEDIUM_COMBINATIONS[i]}"
        
    #     with open(save_path, "wb") as f:
    #         # Build and save the subgraph
    #         subgraph = build_subgraph(path)
    #         pickle.dump(subgraph, f)
    file_path = "data/sihumix/acacae_adam/acacae_adam.smiles_list"
    build_subgraph(file_path, verbose=True)
