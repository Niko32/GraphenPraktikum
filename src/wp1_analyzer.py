import networkx as nx
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import yaml
from typing import List
from matplotlib.colors import LinearSegmentedColormap

import wp1
from constants import SEPCIES_MEDIUM_COMBINATIONS, AMINO_ACIDS, LILA
from custom_types import SpeciesMediumCombination, AminoAcid


def generate_subgraphs(load=True) -> dict[SpeciesMediumCombination, nx.DiGraph]:
    '''
    Generate all subgraphs for all organism and medium combinations and output them in a dictionary
    '''
    
    file_paths = [f"data/sihumix/{c}/{c}.smiles_list" for c in SEPCIES_MEDIUM_COMBINATIONS]

    subgraphs: List[nx.DiGraph] = []
    for i, path in enumerate(file_paths):

        save_path = f"output/subgraphs/{SEPCIES_MEDIUM_COMBINATIONS[i]}"
        
        if load:
            with open(save_path, "rb") as f:
                # Load the subgraphs instead of building them
                subgraph: nx.DiGraph = pickle.load(f)
        else:
            with open(save_path, "wb") as f:
                # Build and save the subgraph
                subgraph = wp1.build_subgraph(path)
                pickle.dump(subgraph, f)
            
        subgraphs.append(subgraph)

    return dict(zip(SEPCIES_MEDIUM_COMBINATIONS, subgraphs))

def generate_pathways() -> dict[SpeciesMediumCombination, dict[AminoAcid, list[str]]]:
    """ Returns the shortest path for every amino acid in every combination of medium and species """
    subgraphs = generate_subgraphs()
    pathways = {}
    print("Finding shortest paths...")
    for combination, s in subgraphs.items():
        amino_acids = set(AMINO_ACIDS).intersection(s.nodes)
        amino_acid_paths = {}
        for a in amino_acids:
            # print(f"Finding shortest simple path for {a} in {combination}")
            shortest_path = [node for node in list(nx.shortest_path(s, "D-glucose", a)) if node.startswith("R_")]
            amino_acid_paths[a] = shortest_path
        pathways[combination] = amino_acid_paths

    # Save the pathways
    with open("pathways.yml", "w") as f:
        yaml.safe_dump(pathways, f)

    return pathways

def compare_nr_amino_acids(pathways: dict[SpeciesMediumCombination, dict[AminoAcid, list[str]]]):
    '''
    Visualization of the number of all amino acids that can be synthesized 
    '''

    number_amino_acids = []
    df = pd.DataFrame([([0] * len(AMINO_ACIDS)) for i in range(len(SEPCIES_MEDIUM_COMBINATIONS))], columns = AMINO_ACIDS, index = SEPCIES_MEDIUM_COMBINATIONS)
    for c in SEPCIES_MEDIUM_COMBINATIONS:
        number_amino_acids.append(len(pathways[c]))
        for aa in AMINO_ACIDS:
            if aa in pathways[c].keys():
                df[aa][c] = 1

    print("Number of amino acids reached for each species medium combination: ")
    print(list(zip(SEPCIES_MEDIUM_COMBINATIONS, number_amino_acids)))

    # barplot for overall number of amino acids that are synthesized per species and medium
    plt.bar(SEPCIES_MEDIUM_COMBINATIONS, number_amino_acids, color=LILA[3])
    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.yticks(np.arange(0, 20, step=2))
    plt.xticks(rotation = 90)
    plt.tight_layout()
    # TODO: Plot beschriftung sch??ner machen
    plt.savefig("output/plots/species_medium_aa_bar_plot.png")

    # heatmap that shows which amino acids are synthesized per combination
    plt.title("synthezised amino acids")
    sns.heatmap(df.transpose(), cbar=False, cmap=["#FFFFFF", LILA[3]])
    # TODO: Plot beschriftung sch??ner machen
    plt.savefig("output/plots/species_medium_aa_heatmap.png", bbox_inches = "tight")

def compare_rec_based_on_organism(pathways: dict[SpeciesMediumCombination, dict[AminoAcid, list[str]]]):
    '''
    For the same organism, are there differences in the reconstruction pathways based on the cultivation media?
    Visualization as a heatmap
    '''

    df = pd.DataFrame(columns = AMINO_ACIDS)

    # iterate over the species
    for i in range(0,len(SEPCIES_MEDIUM_COMBINATIONS),2):
        species = SEPCIES_MEDIUM_COMBINATIONS[i].split('_')[0]
        new_row = pd.DataFrame(np.zeros((1,len(AMINO_ACIDS))), index=[species], columns = AMINO_ACIDS)
        df = pd.concat([df, new_row])
        # get the synthesized amino acid set for each media in the given species
        media1_amino_acids = pathways[SEPCIES_MEDIUM_COMBINATIONS[i]].keys()
        media2_amino_acids = pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]].keys()
        # iterate over the amino acids in media1
        for aa in media1_amino_acids:
            path1 = pathways[SEPCIES_MEDIUM_COMBINATIONS[i]][aa]
            # if the amino acid is also synthesized in media2 count the different reactions in the pathways
            if aa in media2_amino_acids:
                path2 = pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]][aa]
                df[aa][species] = len(list(set(path1).intersection(path2)))
            # if the amino acid is not synthesized in media2 all reactions in the pathway of media1 are different
            else:
                df[aa][species] = len(path1)
        # all amino acids that are synthezised in media2 but not in media 1
        for aa in set(media2_amino_acids).difference(media1_amino_acids):
            path2 = pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]][aa]
            df[aa][species] = len(pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]][aa])
        # all amino acids which are neither synthesized in media1 nor in media2
        for aa in set(AMINO_ACIDS).difference(media1_amino_acids & media2_amino_acids):
            df[aa][species] = np.nan

    # create a heatmap (draw another heatmap, with a transparent color and with only values where the original dataframe is NaN)
    annot_df = df.applymap(lambda f: f'{f:.1f}')
    fig, ax = plt.subplots(squeeze=False)
    sns.heatmap(
        np.where(df.isna(), 0, np.nan),
        ax=ax[0, 0],
        cbar=False,
        annot=np.full_like(df, "NA", dtype=object),
        fmt="",
        annot_kws={"size": 10, "va": "center_baseline", "color": "black"},
        cmap=ListedColormap(['none']),
        linewidth=0)
    sns.heatmap(
        df,
        ax=ax[0, 0],
        cbar=False,
        annot=annot_df,
        fmt="",
        annot_kws={"size": 10, "va": "center_baseline"},
        cmap="coolwarm",
        linewidth=0.5,
        linecolor="black",
        vmin=-1,
        xticklabels=True,
        yticklabels=True)

    plt.tight_layout()
    plt.show()

def path_exists(pathways: list[list[str]], path: list[str]):
    '''
    Checking if path is already in the list pathways
    '''

    for p in pathways:
        if p == path:
            return True
    return False

def compare_rec_based_on_medium(pathways: dict[SpeciesMediumCombination, dict[AminoAcid, list[str]]]):
    '''
    For the same medium, are there notable differences in the reconstructed pathways between species?
    Visualization with a heatmap that shows how many different pathways exist over all species per media and amino acid
    '''

    # data frame with medias in the rows and amino acids in the columns
    df = pd.DataFrame(np.zeros((2, len(AMINO_ACIDS))), columns = AMINO_ACIDS, index = ["adam", "cimIV"])

    for aa in AMINO_ACIDS:
        paths_media1_aa = []
        paths_media2_aa = []
        for i in range(0,16,2):
            if aa in pathways[SEPCIES_MEDIUM_COMBINATIONS[i]].keys():
                if not path_exists(paths_media1_aa, pathways[SEPCIES_MEDIUM_COMBINATIONS[i]][aa]):
                    paths_media1_aa.append(pathways[SEPCIES_MEDIUM_COMBINATIONS[i]][aa])
                if not path_exists(paths_media2_aa, pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]][aa]):
                    paths_media2_aa.append(pathways[SEPCIES_MEDIUM_COMBINATIONS[i+1]][aa])
        df[aa]["adam"] = len(paths_media1_aa)
        df[aa]["cimIV"] = len(paths_media2_aa)

    plt.title("no. of diff. pathways per media in diff. species")
    sns.heatmap(df)
    plt.tight_layout()
    plt.savefig("output/plots/reaction_medium.png")

def alternative_react_paths(pathways: dict[SpeciesMediumCombination, dict[AminoAcid, list[str]]]):
    '''
    Are there/How many alternative reaction paths exist to synthesize each amino acid?
    '''
    # barplot showing the number of alternative paths for each amino acid
    number_of_paths = {}
    for a in AMINO_ACIDS:
        all_paths = []
        for path in pathways.values():
            if a in path.keys():
                all_paths.append(path[a])
            else:
                all_paths.append([])
        all_paths = [path for path in all_paths if path != []]
        number_of_paths[a] = len(np.unique(all_paths))

    plt.bar(AMINO_ACIDS, number_of_paths.values())
    plt.title("no. of alternative reaction paths")
    plt.xlabel("amino acid")
    plt.ylabel("number of amino acids")
    plt.tight_layout()
    plt.show()

def compare_subgraph_sizes(load=False):
    """ Print edges and node set sizes for the original graph and after forward and back traversal """

    # Generate lists of the lengths of the original graph and the graph after bf traversal
    original_lengths, glucose_subgraph_lengths = {}, {}
    for species_medium_combination in SEPCIES_MEDIUM_COMBINATIONS:

        # Build the original graph from .smiles_list files
        G = wp1.build_graph(f"data/sihumix/{species_medium_combination}/{species_medium_combination}.smiles_list")
        original_lengths[species_medium_combination] = len(G)

        # Get the subgraphs after forward traversal
        glucose_subgraph_path = f"output/glucose_subgraphs/{species_medium_combination}.gml"
        if load:
            S = nx.read_gml(glucose_subgraph_path)
        else:
            S = wp1.bf_traversal(G, ["D-glucose"])
            nx.write_gml(S, glucose_subgraph_path)
        glucose_subgraph_lengths[species_medium_combination] = len(S)

    # Generate subgraphs
    subgraphs = generate_subgraphs(load=True)

    # Put them into a dict
    sizes_list = []
    for species_medium_combination, subgraph in subgraphs.items():
        sizes_list.append((
            species_medium_combination, 
            original_lengths[species_medium_combination], 
            glucose_subgraph_lengths[species_medium_combination],
            len(subgraph)
        ))

    return pd.DataFrame(sizes_list, columns=["Species Medium Combination", "Original Graph", "Glucose Graph", "Amino Acids Graph"])

def plot_sizes(sizes_df: pd.DataFrame):
    """ Takes in the sizes_df from compare_subgraph_sizes to create a plot for it """
    sizes_df.index = SEPCIES_MEDIUM_COMBINATIONS
    cmap = LinearSegmentedColormap.from_list("lila", colors=list(LILA.values())[1:]).reversed()
    sizes_df.plot.bar(cmap=cmap)
    plt.savefig("output/plots/subgraph_sizes.png", bbox_inches="tight")

if __name__ == "__main__":
    # pathways = generate_pathways()
    # compare_nr_amino_acids(pathways)
    # compare_rec_based_on_organism(pathways)
    # compare_rec_based_on_medium(pathways)
    # alternative_react_paths(pathways)
    sizes_df = compare_subgraph_sizes(load=True)
    plot_sizes(sizes_df)
    