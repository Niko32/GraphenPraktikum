from typing import List 
import networkx as nx
import numpy as np

from custom_types import AminoAcid, Protein
from constants import proteins

def find_protein_composition(G: nx.DiGraph, protein: Protein) -> dict[AminoAcid, float]:
    """ 
    Get the composition of amino acids of a give protein. 
    Returns a dict containing amino acids as keys and their ratio in the protein as values.
    """

    composition = {}

    raise NotImplementedError

    # Assure the ratios add to one
    assert np.sum(composition.values()) == 1.

    return composition

def parse_proteome(file_path: str) -> list[list[str]]:
    """ Takes in a file_path of the proteome file to output a block of string representing one protein """

def parse_block(protein_block: list[str]) -> dict[str, str]:
    """ Takes a block of strings to output a dict containing the name of the protein and its amnino acid chain """

def get_ratios(proteins: list[dict[str, str]]) -> dict[AminoAcid, float]:
    """ Takes a list of proteins to calculate the ratio of amino acids in the entire proteome """