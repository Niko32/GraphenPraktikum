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

def get_all_paths(g: nx.DiGraph, amino_acid: str) -> list[nx.DiGraph]:
    paths = nx.all_simple_paths(g, "D-glucose", amino_acid)

    return [g.edge_subgraph(list(p)) for p in map(nx.utils.pairwise, paths)]

def generate_pathways(subgraphs: dict[str, nx.DiGraph]) -> dict[str, dict[str, list[nx.DiGraph]]]:
    amino_acid_list = wp1.amino_acid_list
    p = {}
    for c, s in subgraphs.items():
        for a in amino_acid_list:
            p[c] = dict(zip(a, get_all_paths(a, s)))

    return dict(zip(subgraphs.keys, p))

def compare_nr_amino_acids(pathways: dict[str, dict[str, list[nx.DiGraph]]]):
    '''
    Visualization of the number of all amino acids that can be synthesized? 
    '''
    number_amino_acids = []
    for c in combinations:
        number_amino_acids.append(len(pathways[c]))

    plt.bar(combinations, number_amino_acids)
    plt.xlabel("species and media")
    plt.ylabel("number of amino acids")
    plt.show()

    pass

def compare_rec_based_on_organism(pathways: dict[str, dict[str, list[nx.DiGraph]]]):
    '''
    For the same organism, are there differences in the reconstruction pathways based on the cultivation media?
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
                        if nx.graphs_equal(path1, path2):
                            path_found = True
                    if not path_found:
                        df[aa][species] += 1
                for path2 in pathways[combinations[i]][aa]:
                    path_found = False
                    for path1 in pathways[combinations[i+1]][aa]:
                        if nx.graphs_equal(path1, path2):
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

def compare_rec_based_on_medium():
    '''
    For the same medium, are there notable differences in the reconstructed pathways between species?
    '''
    pass

def alternative_react_paths():
    pass

