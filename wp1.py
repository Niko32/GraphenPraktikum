from typing import List, Tuple
import networkx as nx
from matplotlib import pyplot as plt
import logging

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


def bf_traversal(g: nx.DiGraph, metabolite: str) -> nx.DiGraph:
    """ 
    Takes a metabolites to form a subgraph constructed from them
    """
    
    # Assume there is enough energy, water and phosphor in the system
    for cofactor in cofactors:
        g.nodes[cofactor]["visited"] = True

    # Mark the starting node as visited
    g.nodes[metabolite]["visited"] = True

    # Define nodes to search
    outgoing_edges: List[Tuple[str, str]] = g.edges(metabolite, data="visited")
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    # Search every outgoing edge
    while len(outgoing_edges) > 0:
        for u, v, visited in outgoing_edges:

            logging.info(f"Traversing edge ({u},{v})")

            # Check if the node is a reaction or not
            if g.nodes[v].get("reaction"):

                # TODO: edge weights
                # Check if the node has all its required inputs
                ingoing_edges = g.in_edges(v)
                inputs_complete = True
                for u, v in ingoing_edges:
                    if not g.nodes[u].get("visited"):
                        inputs_complete = False

                # Visit the node if it has not been visited yet and has all inputs
                if not visited and inputs_complete:
                    g.nodes[v][visited] = True

                    # Remember the neighbourhood of this node for the next iteration
                    for edge in g.edges(v, data="visited"):
                        new_outgoing_edges.append(edge) 

            # If the next node is not a reaction
            else:
                # Remember the neighbourhood of this node for the next iteration
                for edge in g.edges(v, data="visited"):
                    new_outgoing_edges.append(edge) 

        # Continue with the next iteration
        outgoing_edges = new_outgoing_edges
        new_outgoing_edges = []

    # Construct a subgraph from all nodes that have been visited
    visited_nodes = []
    for node, visited in g.nodes(data="visited"):
        if visited:
            visited_nodes.append(node)
    subgraph = g.subgraph(visited_nodes)

    return subgraph

def reverse_bf_traversal(g: nx.DiGraph):
    
    # Reverse the edges
    g = g.reverse()

    subgraphs = []
    outgoing_edges: List[Tuple[str, str, bool]] = []
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    # Define amino acids as starting points
    for acid in amino_acid_list:
        outgoing_edges = g.edges(acid, data="visited")

        # Search every outgoing edge
        while len(outgoing_edges) > 0:
            for u, v, visited in outgoing_edges:

                # Visit the next node
                if not visited:
                    g.nodes[v][visited] = True

                    # Remember the neighbourhood of this node for the next iteration
                    for edge in g.edges(v).data("visited", default=False):
                        new_outgoing_edges.append(edge) 
            
            # Continue with the next iteration
            outgoing_edges = new_outgoing_edges
            new_outgoing_edges = []

            # Construct a subgraph from all nodes that have been visited
            visited_nodes = []
            for node, visited in g.nodes(data="visited"):
                if visited:
                    visited_nodes.append(node)
            subgraphs.append(g.subgraph(visited_nodes))
            
    return subgraphs

def intersect_subgraph(s: nx.DiGraph, subgraphs: List[nx.DiGraph]) -> nx.DiGraph:
    """ Intersection of the given subgraphs """
    A = subgraphs[0]
    for sub in subgraphs:
        A = nx.intersection(A, sub)

    return nx.intersection(s,A)

def draw_graph(g: nx.DiGraph):
    color_map = []

    for node in g.nodes(data=True):
        if node[1].get("reaction"):
            color_map.append("Blue")
        else:
            color_map.append("Green")
    nx.draw(g, node_color=color_map, with_labels=True)
    plt.show()

def build_subgraph(file_path: str) -> nx.DiGraph:
    """ 
    Builds a subgraph containing the paths from glucose to all amino acids 
    for one combination of molecule and medium 
    """
    reaction_block_list = seperate_blocks(file_path)
    reactions = [extract_compounds(rb) for rb in reaction_block_list]
    g = construct_graph(reactions)
    s = bf_traversal(g, "D-glucose")
    draw_graph(s)
    A = reverse_bf_traversal(s)
    return intersect_subgraph(s,A)
