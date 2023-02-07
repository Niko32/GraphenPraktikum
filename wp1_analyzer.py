import networkx as nx
import wp1
import pickle

def generate_subgraphs() -> dict[str, nx.DiGraph]:
    '''
    Generate all subgraphs for all organism and medium combinations and output them in a dictionary
    '''
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

def get_amino_acids():
    '''
    How many/can all amino acids be synthesized?
    '''
    pass

def compare_rec_based_on_organism():
    '''
    For the same organism, are there differences in the reconstruction based on the cultivation media?
    '''
    pass

def compare_rec_based_on_medium():
    '''
    For the same medium, are there notable differences in the reconstructed pathways between species?
    '''
    pass

def alternative_react_paths():
    pass

