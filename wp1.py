from typing import List, Tuple
import networkx as nx

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

cofactors = ["AMP", "ADP", "ATP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "CTP", "CoA", "H2O", "NH4(+)", "hydrogensulfide"]

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
    # Nicola
    pass

def extract_compounds(reaction_block: List[str]) -> Reaction:
    """
    Takes a block of four string representing one reaction and parses it into
    our reaction class
    """
    # Nicola
    pass

def construct_graph(reactions: List[Reaction]) -> nx.DiGraph:
    """ Takes a list of reations to construct a network x graph from it """
    # Remember to take the reversible flag into account
    # The reaction should be named after the bigg_id
    # Alwine
    pass

def bf_traversal(g: nx.DiGraph, metabolite: str) -> nx.DiGraph:
    """ 
    Takes a metabolites to form a subgraph constructed from them
    """
    
    # Assume there is enough energy in the system

    # Mark the starting node as visited
    g.nodes[metabolite]["visited"] = True

    # Define nodes to search
    outgoing_edges: List[Tuple[str, str, bool]] = g.edges(metabolite).data("visited", default=False)
    new_outgoing_edges: List[Tuple[str, str, bool]] = []

    # Search every outgoing edge
    while len(outgoing_edges) > 0:
        for u, v, visited in outgoing_edges:

            # Check if the node has all its required inputs

            # Visit the node if it has not been visited yet
            if not visited:
                g.nodes[v][visited] = True

                # Remember the neighbourhood of this node for the next iteration
                for edge in g.edges(v).data("visited", default=False):
                    new_outgoing_edges.append(edge) 

        # Continue with the next iteration
        outgoing_edges = new_outgoing_edges
        new_outgoing_edges = []


    # Mark it as visited
    g.nodes[metabolite]["visited"] = True

    # Iterate over the OUTGOING neighborhood of the starting element
    for u,v in g.edges(metabolite):

        # Mark every searched node as visited
        g.nodes[v]["visited"] = True

        # Continue search if all input nodes have already been visited




    

    pass

def reverse_bf_traversal(g: nx.DiGraph):
    return bf_traversal(g.reverse())

def intersect_subgraph(s: nx.DiGraph, subgraphs: List[nx.DiGraph]) -> nx.DiGraph:
    """ Intersection of the given subgraphs """
    # Alwine
    pass

if __name__ == "__main__":
    # 1. Done
    # 2. Parse the file and construct the graph

