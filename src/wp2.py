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