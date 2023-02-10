from typing import List 
import networkx as nx
import numpy as np
import pulp
import pickle

from custom_types import AminoAcid, Protein
from constants import COFACTORS, SEPCIES_MEDIUM_COMBINATIONS, AMINO_ACIDS, AMINO_ACID_DICT



def parse_fasta(file_path: str) -> dict[str, str]:
    """ Takes in a fasta file and outputs a dict containing the name of the protein and its amino acid chain """

    protein_dict = {}

    with open(file_path, 'r') as f:
        for l in f:
            if l[0] == '>':
                # save the curent protein name and create a dict entry
                curent_protein = l[1:].rstrip('\n')
                protein_dict[curent_protein] = ""
            else:
                # append the sequence line to the sequence of the curently read in protein
                protein_dict[curent_protein] = protein_dict[curent_protein] + l.rstrip('\n')

    return protein_dict


def get_ratios(proteins: dict[str, str]) -> dict[AminoAcid, float]:
    """ Takes a list of proteins to calculate the ratio of amino acids in the entire proteome """

    # first set number of each amino acid to 0
    aminoacid_ratio = {}
    for aa in AMINO_ACIDS:
        aminoacid_ratio[aa] = 0

    # iterate over all proteins sum up the length of each protein and count the amino acids
    proteom_len = 0
    for p in proteins:
        proteom_len += len(proteins[p])
        for aa in proteins[p]:
            aminoacid_ratio[AMINO_ACID_DICT[aa]] += 1

    # divide the number of amino acids through the length of the proteom
    for aa in AMINO_ACIDS:
        aminoacid_ratio[aa] = aminoacid_ratio[aa] / proteom_len

    return aminoacid_ratio

def add_biomass_reaction(G: nx.DiGraph, ratios: dict[AminoAcid, float]) -> nx.DiGraph:
    """ Takes the subgraph containing all amino acids and adds a biomass reaction and node to it """

    # Add the biomass reaction and output node
    G.add_node("R_biomass", reaction=True)
    G.add_node("biomass")
    G.add_edge("R_biomass", "biomass")

    # Add edges from the amino acids to the biomass node
    amino_acids_in_the_graph = set(AMINO_ACIDS).intersection(G.nodes)
    for amino_acid in amino_acids_in_the_graph:
        G.add_edge(amino_acid, "R_biomass", weight=ratios[amino_acid])

    return G

def add_input_reactions(G: nx.DiGraph) -> nx.DiGraph:
    """ Take the graph and adds input and output reactions to it """

    # Add a node for every cofactor and connect it to the cofactor
    for cofactor in [*COFACTORS, "D-glucose"]:
        input_reaction = f"R_{cofactor}"
        G.add_node(input_reaction)
        G.add_edge(input_reaction, cofactor, weight=1)

    return G

def load_graph():
    save_path = f"output/subgraphs/{SEPCIES_MEDIUM_COMBINATIONS[0]}"
    with open(save_path, "rb") as f:
        G: nx.DiGraph = pickle.load(f)
        return G
    
def get_variables(G: nx.DiGraph) -> dict[str, pulp.LpVariable]:
    """ Get the pulp variables representing the reactions in a dict with the reaction name """

    variables = {}

    # Add a reaction variable for every reaction in the graph
    for node in G.nodes:
        if node.startswith("R_"):
            var = pulp.LpVariable(node, lowBound=0, upBound=50, cat='Integer')
            variables[node] = var

    return variables

def add_constraints(model: pulp.LpProblem, V: dict[str: pulp.LpVariable], G: nx.DiGraph) -> pulp.LpProblem:
    """ Add the constraints to the model  including input constraints and a glucose constraint """
    # Add glucose constraint
    model += V["R_D-glucose"] >= 10

    # Add constraints for each cofactors
    for cofactor in COFACTORS:
        model += V["R_" + cofactor] >= -1000

    # Add constraints for all inner nodes
    for v in V.keys():
        if not v in COFACTORS:
            model += V[v] >= 0

    # Add reaction equations iterating over all compound nodes
    nodes = G.nodes(data="reaction")
    for node in nodes:
        if not node.get("reaction"):
            predessecors = [G.edges[node, p]["weight"] * V[p] for p in G.predecessors(node)]
            successors = [G.edges[node, p]["weight"] * V[p] for p in G.successors(node)]

            model += pulp.lpSum(predessecors) == pulp.lpSum(successors)

    return model

if __name__ == "__main__":
    proteins = parse_fasta("data/proteomes/Anaerostipes_caccae.faa")
    ratios = get_ratios(proteins)
    G = load_graph()
    G = add_biomass_reaction(G, ratios)
    G = add_input_reactions(G)
    model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)
    variables = get_variables(G)
    model += variables["R_biomass"], "Profit"
    model = add_constraints(model, variables, G)
    model.solve()
    pulp.LpStatus[model.status]