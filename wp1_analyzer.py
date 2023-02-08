import networkx as nx
import wp1
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap

combinations = [
        "acacae_adam",
        "acacae_cimIV",
        "blongum_adam",
        "blongum_cimIV",
        "bproducta_adam",
        "bproducta_cimIV",
        "btheta_adam",
        "btheta_cimIV",
        "cbuty_adam",
        "cbuty_cimIV",
        "ecoli_adam",
        "ecoli_cimIV",
        "eramosum_adam",
        "eramosum_cimIV",
        "lplantarum_adam",
        "lplantarum_cimIV"
    ]

def generate_subgraphs() -> dict[str, nx.DiGraph]:
    '''
    Generate all subgraphs for all organism and medium combinations and output them in a dictionary
    '''
    
    file_paths = ["sihumix/" + c + "/" + c + ".smiles_list" for c in combinations]

    subgraphs = []
    for i, path in enumerate(file_paths):
        subgraph = wp1.build_subgraph(path)

        # Save the subgraph
        with open(f"subgraphs/{combinations[i]}", "wb") as f:
            pickle.dump(subgraph, f)
            
        subgraphs.append(subgraph)

    return dict(zip(combinations, subgraphs))

def generate_pathways() -> dict[str, dict[str, list[str]]]:
    subgraphs = generate_subgraphs()
    amino_acid_list = wp1.amino_acid_list
    p = {}
    for c, s in subgraphs.items():
        for a in amino_acid_list:
            p[c] = dict(zip(a, list(nx.shortest_simple_paths(s, "D-glucose", a))[0]))

    return dict(zip(subgraphs.keys, p))

def compare_nr_amino_acids(pathways: dict[str, dict[str, list[str]]]):
    '''
    Visualization of the number of all amino acids that can be synthesized 
    '''

    number_amino_acids = []
    df = pd.DataFrame([([0] * len(wp1.amino_acid_list)) for i in len(combinations)], columns = wp1.amino_acid_list, index = combinations)
    for c in combinations:
        number_amino_acids.append(len(pathways[c]))
        for aa in wp1.amino_acid_list:
            if aa in pathways[c].keys():
                df[aa][c] = 'X'

    # barplot for overall number of amino acids that are synthesized per species and medium
    plt.bar(combinations, number_amino_acids)
    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.show()

    # table that shows which amino acids are synthesized per combination
    df.head()

def compare_rec_based_on_organism(pathways: dict[str, dict[str, list[str]]]):
    '''
    For the same organism, are there differences in the reconstruction pathways based on the cultivation media?
    Visualization as a heatmap
    '''

    df = pd.DataFrame(np.zeros((8, len(wp1.amino_acid_list))), columns = wp1.amino_acid_list)

    # iterate over the species
    for i in range(0,len(combinations),2):
        species = combinations[i].split('_')[0]
        df.rename(index = {i/2 : species})
        # get the synthesized amino acid set for each media in the given species
        media1_amino_acids = pathways[combinations[i]].keys()
        media2_amino_acids = pathways[combinations[i+1]].keys()
        # iterate over the amino acids in media1
        for aa in media1_amino_acids:
            if aa in media2_amino_acids:
                # iterate over the paths to the amino acid in the different medias and check if they exist is both medias
                for path1 in pathways[combinations[i]][aa]:
                    path_found = False
                    for path2 in pathways[combinations[i+1]][aa]:
                        if nx.utils.graphs_equal(path1, path2):
                            path_found = True
                    if not path_found:
                        df[aa][species] += 1
                for path2 in pathways[combinations[i]][aa]:
                    path_found = False
                    for path1 in pathways[combinations[i+1]][aa]:
                        if nx.utils.graphs_equal(path1, path2):
                            path_found = True
                    if not path_found:
                        df[aa][species] += 1
            # if the amino acid is not synthesized in media 2 all pathways are different
            else:
                df[aa][species] = len(pathways[combinations[i]][aa])
        for aa in media2_amino_acids:
            # if the amino acid is not synthesized in media 1 all pathways are different
            if not aa in media1_amino_acids:
                df[aa][species] = len(pathways[combinations[i+1]][aa])
        for aa in set(wp1.amino_acid_list).difference(media1_amino_acids & media2_amino_acids):
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
        vmax=1,
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

def compare_rec_based_on_medium(pathways: dict[str, dict[str, list[nx.DiGraph]]]) -> list[dict[str, list]]:
    '''
    For the same medium, are there notable differences in the reconstructed pathways between species?
    Save the number of pairwise different pathways for each amino acid and medium in a matrix which coulms and rows represent the different species
    '''

    # dictionaries with a cubic matrix for every amino acid
    df = pd.DataFrame(np.zeros((2, len(wp1.amino_acid_list))), columns = wp1.amino_acid_list, index = ["adam", "cimIV"])

    for aa in wp1.amino_acid_list:
        paths_media1_aa = []
        paths_media2_aa = []
        for i in range(0,16,2):
            if aa in pathways[combinations[i]].keys():
                if not path_exists(paths_media1_aa, pathways[combinations[i]][aa]):
                    paths_media1_aa.append(pathways[combinations[i]][aa])
                if not path_exists(paths_media2_aa, pathways[combinations[i+1]][aa]):
                    paths_media2_aa.append(pathways[combinations[i+1]][aa])
        df[aa]["adam"] = len(paths_media1_aa)
        df[aa]["cimIV"] = len(paths_media2_aa)

    plt.title("no. of diff. pathways per media in diff. species")
    sns.heatmap(df)
    plt.show()

def alternative_react_paths():
    '''
    Are there/How many alternative reaction paths exist to synthesize each amino acid?
    '''
    # barplot showing the number of alternative paths for each amino acid
    number_of_paths = {}
    for a in wp1.amino_acid_list:
        all_paths = []
        pathway_dict = generate_pathways()
        for p in pathway_dict.values():
            all_paths.append(p[a])
        number_of_paths[a] = len(np.unique(all_paths))

    plt.bar(wp1.amino_acid_list, number_of_paths.values())
    plt.title("no. of alternative reaction paths")
    plt.xlabel("amino acid")
    plt.ylabel("number of amino acids")
    plt.show()

if __name__ == "__main__":
    pathways = generate_pathways()
    compare_nr_amino_acids(pathways)
    compare_rec_based_on_organism(pathways)
    media1, media2 = compare_rec_based_on_medium(pathways)
    