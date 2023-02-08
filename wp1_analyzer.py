import networkx as nx
import wp1
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
        "Iplantarum_adam",
        "Iplantarum_cimIV"
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
    subgraphs = generate_pathways()
    amino_acid_list = wp1.amino_acid_list
    p = {}
    for c, s in subgraphs.items():
        for a in amino_acid_list:
            p[c] = dict(zip(a, nx.shortest_simple_paths(s, "D-glucose", a)[0]))

    return dict(zip(subgraphs.keys, p))

def compare_nr_amino_acids(pathways: dict[str, dict[str, list[nx.DiGraph]]]):
    '''
    Visualization of the number of all amino acids that can be synthesized 
    '''

    number_amino_acids = []
    df = pd.DataFrame(np.zeros((len(combinations), len(wp1.amino_acid_list))), columns = wp1.amino_acid_list, index = combinations)
    for c in combinations:
        number_amino_acids.append(len(pathways[c]))
        for aa in wp1.amino_acid_list:
            if aa in pathways[c].keys():
                df[aa][c] = len(pathways[c][aa])

    # barplot for overall number of amino acids that are synthesized per species and medium
    plt.bar(combinations, number_amino_acids)
    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.show()

    # barplot that shows the number of pathways per amino acid and combination
    plt.title("no. of pathways")
    sns.heatmap(df)
    plt.show()

    pass

def compare_rec_based_on_organism(pathways: dict[str, dict[str, list[nx.DiGraph]]]):
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

    plt.title("no. of diff. pathways per species in diff. medias")
    sns.heatmap(df)
    plt.show()

    pass

def compare_pathways_species(media: dict[str, list[list]], species1: int, species2: int, pathways: dict[str, dict[str, list[nx.DiGraph]]]):
    '''
    Counting different pathways for each amino acid in the species1 and species and saving in the dict media
    '''
    amino_acids1 = set(pathways[combinations[species1]].keys())
    amino_acids2 = set(pathways[combinations[species2]].keys())
    # iterate over all amino acids that are synthesized in species i and species j
    for aa in amino_acids1 & amino_acids2:
        # iterate over the paths to the amino acid in the different medias and check if they exist is both medias
        for path1 in pathways[combinations[species1]][aa]:
            path_found = False
            for path2 in pathways[combinations[species2]][aa]:
                if nx.utils.graphs_equal(path1, path2):
                    path_found = True
            if not path_found:
                media[aa][species1][species2] += 1
        for path2 in pathways[combinations[species1]][aa]:
            path_found = False
            for path1 in pathways[combinations[species2]][aa]:
                if nx.utils.graphs_equal(path1, path2):
                    path_found = True
            if not path_found:
                media[aa][species1][species2] += 1
    # iterate over all amino acids that are only synthesized in species i
    for aa in amino_acids1.difference(amino_acids2):
        media[aa][species1][species2] = len(pathways[combinations[species1]][aa])
    # iterate over all amino acids that are only synthesized in species j
    for aa in amino_acids2.difference(amino_acids1):
        media[aa][species1][species2] = len(pathways[combinations[species2]][aa])
    # iterate over all amino acids that are not synthesized neither in i nor in j
    for aa in set(wp1.amino_acid_list).difference(amino_acids1 & amino_acids2):
        media[aa][species1][species2] = 'X'

def compare_rec_based_on_medium(pathways: dict[str, dict[str, list[nx.DiGraph]]]) -> list[dict[str, list]]:
    '''
    For the same medium, are there notable differences in the reconstructed pathways between species?
    Save the number of pairwise different pathways for each amino acid and medium in a matrix which coulms and rows represent the different species
    '''

    # dictionaries with a cubic matrix for every amino acid
    media1 = {k:np.zeros((8, 8)) for k in wp1.amino_acid_list}
    media2 = {k:np.zeros((8, 8)) for k in wp1.amino_acid_list}

    for i in range(0,15,2):
        for j in range(i+2,16,2):

            compare_pathways_species(media1, i, j, pathways)
            compare_pathways_species(media2, i+1, j+1, pathways)

    return [media1, media2]

def alternative_react_paths():
    pass

