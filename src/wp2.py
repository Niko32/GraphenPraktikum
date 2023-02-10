from typing import List 
import networkx as nx
import numpy as np
import pulp
import pickle

from custom_types import AminoAcid, Protein
from constants import PROTEINS, SEPCIES_MEDIUM_COMBINATIONS


def parse_fasta(file_path: str) -> dict[str, str]:
    """ Takes in a fasta file and outputs a dict containing the name of the protein and its amino acid chain """

def get_ratios(proteins: dict[str, str]) -> dict[AminoAcid, float]:
    """ Takes a list of proteins to calculate the ratio of amino acids in the entire proteome """

def add_biomass_reaction(G: nx.DiGraph, ratios: dict[AminoAcid, float]) -> nx.DiGraph:
    """ Takes the subgraph containing all amino acids and adds a biomass reaction and node to it """

def load_graph():
    save_path = f"output/subgraphs/{SEPCIES_MEDIUM_COMBINATIONS[0]}"
    with open(save_path, "rb") as f:
        G: nx.DiGraph = pickle.load(f)
        return G
    
def get_variables(G: nx.DiGraph) -> dict[str, pulp.LpVariable]:
    """ Get the pulp variables representing the reactions in a dict with the reaction name """

def add_constraints(model: pulp.LpProblem, G: nx.DiGraph) -> pulp.LpProblem:
    """ Add the constraints to the model  including input constraints and a glucose constraint """


if __name__ == "__main__":
    proteins = parse_fasta("data/proteomes/Anaerostipes_caccae.faa")
    ratios = get_ratios(proteins)
    G = load_graph()
    model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)
    variables = get_variables(G)
    model += variables["R_biomass"], "Profit"
    model = add_constraints(model, G)
    model.solve()
    pulp.LpStatus[model.status]