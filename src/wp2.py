from typing import List 
import networkx as nx
import numpy as np
import pulp
import pickle
import pandas as pd
from matplotlib import pyplot as plt

from custom_types import AminoAcid, Protein
from constants import COFACTORS, SEPCIES_MEDIUM_COMBINATIONS, AMINO_ACIDS, AMINO_ACID_DICT, SPECIES_DICT
from wp1 import bf_traversal, draw_graph

START_NODE = "D-glucose"


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
    for protein in proteins.values():
        protein = protein.replace("U", "")
        proteom_len += len(protein)
        for aa in protein:
            aminoacid_ratio[AMINO_ACID_DICT[aa]] += 1

    # divide the number of amino acids through the length of the proteom
    for aa in AMINO_ACIDS:
        aminoacid_ratio[aa] = aminoacid_ratio[aa] / proteom_len

    return aminoacid_ratio

def add_biomass_reactions(G: nx.DiGraph, ratios: dict[AminoAcid, float]) -> nx.DiGraph:
    """ Takes the subgraph containing all amino acids and adds a biomass reaction and node to it """

    # Add the biomass reaction and output node
    G.add_node("biomass")

    # Add edges from the amino acids to the biomass node
    amino_acids_in_the_graph = set(AMINO_ACIDS).intersection(G.nodes)
    for amino_acid in amino_acids_in_the_graph:
        G.add_node(f"R_{amino_acid}", reaction=True)
        G.add_edge(amino_acid, f"R_{amino_acid}", weight=1)
        G.add_edge(f"R_{amino_acid}", "biomass", weight=ratios[amino_acid])

    return G

def add_input_reactions(G: nx.DiGraph) -> nx.DiGraph:
    """ Take the graph and adds input and output reactions to it """

    cofactors = set(COFACTORS).intersection(G.nodes)

    # Add a node for every cofactor and connect it to the cofactor
    for cofactor in [*cofactors, START_NODE]:
        input_reaction = f"R_in_{cofactor}"
        G.add_node(input_reaction, reaction=True)
        G.add_edge(input_reaction, cofactor, weight=1)

    return G

def add_output_reactions(G: nx.DiGraph):
    compounds = []
    for compound, is_reaction in G.nodes(data="reaction"):

        if is_reaction:
            continue

        compounds.append(compound)

    for compound in compounds:
        G.add_node(f"R_out_{compound}", reaction=True)
        G.add_edge(compound, f"R_out_{compound}", weight=1)

    return G

def load_graph(species_medium_combination: str):
    save_path = f"output/subgraphs/{species_medium_combination}"
    with open(save_path, "rb") as f:
        G: nx.DiGraph = pickle.load(f)
        return G
    
def get_variables(G: nx.DiGraph) -> dict[str, pulp.LpVariable]:
    """ Get the pulp variables representing the reactions in a dict with the reaction name """

    variables = {}

    # Add a reaction variable for every reaction in the graph
    for node in G.nodes:
        if node.startswith("R_"):
            var = pulp.LpVariable(node, cat='Integer')
            variables[node] = var

    return variables

def add_constraints(model: pulp.LpProblem, V: dict[str: pulp.LpVariable], G: nx.DiGraph) -> pulp.LpProblem:
    """ Add the constraints to the model  including input constraints and a glucose constraint """
    # Add glucose constraint
    model += V[f"R_in_{START_NODE}"] >= 10
    model += V[f"R_in_{START_NODE}"] <= 1000

    # Add constraints for each cofactors
    cofactors = set(COFACTORS).intersection(G.nodes)
    for cofactor in cofactors:
        model += V["R_in_" + cofactor] >= -1000
        model += V["R_in_" + cofactor] <= 1000

    # Add constraints for all inner nodes
    for v in V.keys():
        if not v in [f"R_in_{cofactor}" for cofactor in cofactors]:
            model += V[v] >= 0
            model += V[v] <= 1000

    # Add reaction equations iterating over all compound nodes
    for compound, is_reaction in G.nodes(data="reaction"):

        # Skip reactions
        if is_reaction:
            continue

        predessecors = [data["weight"] * V[u] for u, v, data in G.in_edges(compound, data=True)]
        successors = [data["weight"] * V[v] for u, v, data in G.out_edges(compound, data=True)]

        c = pulp.lpSum(predessecors) >= pulp.lpSum(successors)
        model += c

    return model

def create_models():
    proteomes = {key: parse_fasta(f"data/proteomes/{value}.faa") for key, value in SPECIES_DICT.items()}
    for species_medium_combination in SEPCIES_MEDIUM_COMBINATIONS:
        species = "_".join(species_medium_combination.split("_")[:-1])
        proteome = proteomes[species]
        ratios = get_ratios(proteome)
        G = load_graph(species_medium_combination)
        G = add_biomass_reactions(G, ratios)
        G = add_input_reactions(G)
        G = add_output_reactions(G)
        #draw_graph(G, show_reactions=True)
        model = pulp.LpProblem(species_medium_combination, pulp.LpMaximize)
        variables = get_variables(G)
        model += variables["R_out_biomass"], "Profit"
        model = add_constraints(model, variables, G)
        model.solve()

        # Print results
        print(f"Solved problem for {species_medium_combination}")
        print("Model Objective:", pulp.value(model.objective))

        for v in variables.values():
            if v.varValue > 0:
                print(v.name,":", v.varValue)

        # Save model
        with open(f"output/models/{species_medium_combination}.pkl", "wb") as f:
            pickle.dump(model, f)

def plot_biomasses(df: pd.DataFrame):
    
    # Fixing random state for reproducibility
    np.random.seed(19680801)

    plt.rcdefaults()
    fig, ax = plt.subplots()
    plt.tight_layout()

    # Example data
    y_pos = np.arange(len(df))

    ax.barh(y_pos, df["Objective Function"], align='center')
    ax.set_yticks(y_pos, labels=list(df["Species Medium Combination"]))
    ax.invert_yaxis()
    ax.set_xlabel('Biomass')
    ax.set_title('Which medium produces optimal amino acid distributions?')

    plt.savefig("output/plots/biomasses.png")

if __name__ == "__main__":
    create_models(),

    # Create an empty dataframe
    models_df = pd.DataFrame([], columns=["Species Medium Combination", "Objective Function"])

    # Load the model
    for species_medium_combination in SEPCIES_MEDIUM_COMBINATIONS:
        with open(f"output/models/{species_medium_combination}.pkl", "rb") as f:
            model: pulp.LpProblem = pickle.load(f)

        # Add an entry for the species medium combination
        new_entry = pd.DataFrame({"Species Medium Combination": species_medium_combination, "Objective Function": pulp.value(model.objective)}, index=[0])
        models_df = pd.concat([models_df, new_entry], ignore_index=True)

    plot_biomasses(models_df)

