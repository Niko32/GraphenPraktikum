from typing import List
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
    # Niko
    pass

def reverse_bf_traversal():
    pass
    # Niko

def intersect_subgraph(s: nx.DiGraph, subgraphs: List[nx.DiGraph]) -> nx.DiGraph:
    """ Intersection of the given subgraphs """
    # Alwine
    pass