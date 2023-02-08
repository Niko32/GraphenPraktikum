from typing import List, Tuple, Union
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image
import os

amino_acid_list = ["L-arginine", "L-valine", "L-methionine", "L-glutamate", "L-glutamine", "L-tyrosine", "L-tryptophan",
                   "L-proline", "L-cysteine", "L-histidine", "L-asparagine", "L-aspartate", "L-phenylalanine",
                   "L-threonine", "L-lysine", "L-serine", "L-isoleucine", "glycine", "L-alanine", "L-Leucine"]

amino_acid_dict = {
    "R": {"name": "L-arginine"},
    "V": {"name": "L-valine"},
    "M": {"name": "L-methionine"},
    "E": {"name": "L-glutamate"},
    "Q": {"name": "L-glutamine"},
    "Y": {"name": "L-tyrosine"},
    "W": {"name": "L-tryptophan"},
    "P": {"name": "L-proline"},
    "C": {"name": "L-cysteine"},
    "H": {"name": "L-histidine"},
    "N": {"name": "L-asparagine"},
    "D": {"name": "L-aspartate"},
    "F": {"name": "L-phenylalanine"},
    "T": {"name": "L-threonine"},
    "K": {"name": "L-lysine"},
    "S": {"name": "L-serine"},
    "I": {"name": "L-isoleucine"},
    "G": {"name": "glycine"},
    "A": {"name": "L-alanine"},
    "L": {"name": "L-Leucine"}
}

cofactors = ["AMP", "ADP", "ATP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "CTP", "CoA", "H2O", "NH4(+)", "hydrogen sulfide"]

class Reaction():
    def __init__(self, bigg_id: str, metanetx_id: str, reversible: bool, educts: List[str],
                 products: List[str], smiles_educts: List[str], smiles_products: List[str]):
        self.bigg_id = bigg_id
        self.metanetx_id = metanetx_id
        self.reversible = reversible
        self.products = products
        self.educts = educts
        self.smiles_educts = smiles_educts
        self.smiles_products = smiles_products


def seperate_blocks(file_path: str) -> List[List[str]]:
    """ 
    Takes in the .smiles_list file and outputs the reaction blocks 
    as a list of lists of four strings
    """

    reaction_blocks = []
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

    return reaction_blocks
    

def extract_compounds(reaction_block: List[str]) -> Reaction:
    """
    Takes a block of four string representing one reaction and parses it into
    our reaction class
    """
    line1 = reaction_block[0].split()
    bigg = line1[2]
    metanet = line1[4]
    reverse = (line1[6] == "True\n")

    line3_educts = reaction_block[2].split('=')[0].strip()
    line3_products = reaction_block[2].split('=')[1].strip()
    educts = line3_educts.split(" + ")
    products = line3_products.split(" + ")

    line4_educts = reaction_block[3].split('>>')[0].strip()
    line4_products = reaction_block[3].split('>>')[1].strip()
    sm_educts = line4_educts.split(".")
    sm_products = line4_products.split(".")

    reac = Reaction(bigg, metanet, reverse, educts, products, sm_educts, sm_products)

    return reac


def construct_graph(reactions: List[Reaction]) -> nx.DiGraph:
    """ Takes a list of reactions to construct a network x graph from it """
    G = nx.DiGraph()
    for r in reactions:
        G.add_node(r.bigg_id, reaction=True)
        for e, smiles_e in zip(r.educts, r.smiles_educts):
            weight = r.educts.count(e)
            G.add_node(e, smiles=smiles_e)
            G.add_edge(e, r.bigg_id, weight=weight)

            if r.reversible:
                G.add_edge(r.bigg_id, e, weight=weight)

        for p, smiles_p in zip(r.products, r.smiles_products):
            weight = r.products.count(p)
            G.add_node(p, smiles=smiles_p)
            G.add_edge(r.bigg_id, p, weight=weight)

            if r.reversible:
                G.add_edge(p, r.bigg_id, weight=weight)

    return G

def _visit_node(G: nx.DiGraph, v: str) -> Union[nx.DiGraph, List[Tuple[str, str, bool]]]:
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
    images = [Image.open(path) for path in image_paths]
    images[0].save("plots/traversal.gif", save_all = True, 
                   append_images = images[1:], optimize = False, 
                   duration = len(image_paths)*5)

def bf_traversal(G: nx.DiGraph, metabolite: str, verbose = False) -> nx.DiGraph:
    """ 
    Takes a metabolites to form a subgraph constructed from them
    """

    print("##### bf_traversal #####")
    
    # Assume there is enough energy, water and phosphor in the system
    for cofactor in cofactors:
        G.nodes[cofactor]["visited"] = True

    # Mark the starting node as visited
    G.nodes[metabolite]["visited"] = True

    # Define nodes to search
    outgoing_edges: List[Tuple[str, str, dict]] = G.edges(metabolite, data="visited")
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    iteration = 0

    # Search every outgoing edge
    while len(outgoing_edges) > 0:
        print(f"### Iteration {iteration} ###")

        # Save the graph figure
        if verbose:
            s = _build_visited_subgraph(G)
            draw_graph(s, f"plots/{iteration:02}.png")

        for u, v, visited in outgoing_edges:

            print(f"Traversing edge ({u},{v})")

            # Skip visiting the node if it does not have all required inputs
            if not _inputs_complete(G, v):
                continue

            # Visit the node
            G, edges = _visit_node(G, v)
            new_outgoing_edges.extend(edges)

        # Continue with the next iteration
        iteration += 1
        outgoing_edges = new_outgoing_edges
        new_outgoing_edges = []

    # Save the graph figure
    s = _build_visited_subgraph(G)

    if verbose:
        draw_graph(s, f"plots/{iteration:02}.png")
        _create_gif([f"plots/{file}" for file in os.listdir("plots")])

    # Print out amino acids that have been reached
    reached_acids = set(s.nodes.keys()).intersection(amino_acid_list)
    print(len(reached_acids), "/", len(amino_acid_list), "acids have been reached")

    return s

def reverse_bf_traversal(G: nx.DiGraph) -> List[nx.DiGraph]:

    print("##### reverse_bf_traversal #####")
    
    # Reverse the edges
    G = G.reverse()

    # Unvisit the nodes
    nx.set_node_attributes(G, False, "visited")

    subgraphs = []
    outgoing_edges: List[Tuple[str, str, bool]] = []
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    # Create one subgraph for every amino acid
    for acid in amino_acid_list:

        # Skip acids are not in the graph
        if not G.nodes.get(acid):
            continue

        # Mark the amino acid as visited
        G, outgoing_edges = _visit_node(G, acid)

        # Search every outgoing edge
        while len(outgoing_edges) > 0:
            for u, v, visited in outgoing_edges:

                # Visit the next node
                G, edges = _visit_node(G, v)
                new_outgoing_edges.extend(edges) 
            
            # Continue with the next iteration
            outgoing_edges = new_outgoing_edges
            new_outgoing_edges = []

            # Construct a subgraph from all nodes that have been visited
            s = _build_visited_subgraph(G)

            subgraphs.append(s)
            
    return subgraphs

def intersect_subgraph(s: nx.DiGraph, subgraphs: List[nx.DiGraph]) -> nx.DiGraph:
    """ Intersection of the given subgraphs """
    A = subgraphs[0]
    for sub in subgraphs:
        A = nx.union(A, sub)

    return nx.intersection(s,A)

def draw_graph(G: nx.DiGraph, output = ""):

    # Configure display settings based on the reaction flag
    nodes, sizes, labels, color_map = len(G)*[None], np.empty(len(G), dtype=int), {}, len(G)*[None]
    for i, (node, reaction_flag) in enumerate(G.nodes(data="reaction")):
        nodes[i] = node
        sizes[i] = 50 - 49 * (reaction_flag or 0)
        labels[node] = "" if reaction_flag else node
        color_map[i] = "Blue" if reaction_flag else "Green"

    nx.draw(G, pos=nx.layout.kamada_kawai_layout(G), node_color=color_map, with_labels=True, 
            nodelist=nodes, node_size=sizes, font_size=5, width=0.1, arrowsize=4, labels=labels)
    
    # Save the figure
    if output:
        plt.savefig(output)
    # Or show it
    else:
        plt.show()

    # Clear the figure
    plt.clf()

def build_subgraph(file_path: str) -> nx.DiGraph:
    """ 
    Builds a subgraph containing the paths from glucose to all amino acids 
    for one combination of molecule and medium 
    """
    reaction_block_list = seperate_blocks(file_path)
    reactions = [extract_compounds(rb) for rb in reaction_block_list]
    G = construct_graph(reactions)
    s = bf_traversal(G, "D-glucose", True)
    # Why is not a single amino aicid reached?
    draw_graph(s)
    A = reverse_bf_traversal(s)
    draw_graph(A[0])
    return intersect_subgraph(s,A)

if __name__ == "__main__":
    build_subgraph("sihumix/acacae_adam/acacae_adam.smiles_list")
