import networkx as nx
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import yaml
from typing import List

import wp1
from constants import SEPCIES_MEDIUM_COMBINATIONS, AMINO_ACIDS
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
    for combination, s in subgraphs.items():
        amino_acids = set(AMINO_ACIDS).intersection(s.nodes)
        amino_acid_paths = {}
        for a in amino_acids:
            print(f"Finding shortest simple path for {a} in {combination}")
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

    # barplot for overall number of amino acids that are synthesized per species and medium
    plt.bar(SEPCIES_MEDIUM_COMBINATIONS, number_amino_acids)
    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.yticks(np.arange(0, 20, step=2))
    plt.xticks(rotation = 90)
    plt.show()

    # heatmap that shows which amino acids are synthesized per combination
    plt.title("synthezised amino acids")
    sns.heatmap(df, cbar=False)
    plt.show()

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
    plt.show()

def path_exists(pathways: list[list[str]], path: list[str]):
    '''
    Checking if path is already in the list pathways
    '''

    for p in pathways:
        if p == path:
            return True
    return False

def compare_rec_based_on_medium(pathways: dict[SpeciesMediumCombination, dict[AminoAcid, list[nx.DiGraph]]]):
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
    plt.show()

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
    plt.show()

if __name__ == "__main__":
    pathways = generate_pathways()
    compare_nr_amino_acids(pathways)
    compare_rec_based_on_organism(pathways)
    compare_rec_based_on_medium(pathways)
    alternative_react_paths(pathways)
    