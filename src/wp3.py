from typing import List, Tuple
import re
import networkx as nx

from custom_types import Atom

def seperate_blocks(file_path: str) -> Tuple[List[List[str]], List[List[str]]]:
    '''
    Seperate component blocks from the gml file
    '''
    
    # first list contains all node blocks and second list contains all edge blocks
    component_blocks = [[],[]]

    with open(file_path, 'r') as f:
        # avoid reading first and last line
        for l in f.readlines()[1:len(f.readlines())-2]:
            if "node [" in l:
                component = []
                component_type = 0
            elif "edge [" in l:
                component = []
                component_type = 1
            elif "]" in l:
                component_blocks[component_type].append(component)
            else:
                component.append(l.strip())

    return component_blocks


def extract_atom(node_block: List[str]) -> Atom:
    """
    Takes a block of four string representing one node and parses it into
    our node class
    """
    
    new_atom: Atom = {}
    new_atom["id"] = int(node_block[0].strip().split(" "))
    new_atom["label"] = re.findall(r'"([^"]*)"', reaction_block[1])[0]
    new_atom["charge"] = int(node_block[2].strip().split(" "))
    new_atom["hcount"] = int(node_block[3].strip().split(" "))
    new_atom["aromatic"] = int(node_block[4].strip().split(" "))
    new_atom["element"] = re.findall(r'"([^"]*)"', reaction_block[5])[0]
    new_atom["atom_class"] = int(node_block[6].strip().split(" "))
    new_atom["compound_id"] = int(node_block[7].strip().split(" "))
    new_atom["compound_name"] = re.findall(r'"([^"]*)"', reaction_block[8])[0]

    return new_atom


def construct_graph(components: Tuple[List[List[str]], List[List[str]]]) -> nx.DiGraph:
    """ 
    Takes a list of components to construct a network x graph from it
    """
    G = nx.DiGraph()

    atoms = {}
    # creates the atom nodes 
    for atom_block in components[0]:
        new_atom = extract_atom(atom_block)
        atoms[new_atom["id"]] = new_atom["label"]
        G.add_node(new_atom["label"], reaction=False, atom=new_atom)

    # creates the reaction nodes and the edges
    for edge_block in components[1]:
        source = int(edge_block[0].strip().split(" "))
        target = int(edge_block[1].strip().split(" "))
        # checks if the block represents a reaction
        if "transition" in edge_block[2]:
            reaction_id = "R_" + str(re.findall(r'"([^"]*)"', edge_block[3])[0])
            # checks if the reaction already exists and only more educts and products have to be added
            if not G.has_node(reaction_id):
                G.add_node(reaction_id, reaction=True)
            else:
                G.add_edge(atoms[source]["label"], reaction_id)
                G.add_edge(reaction_id, atoms[target]["label"])
        # create an edge in the compound
        else:
            G.add_edge(atoms[source]["label"], atoms[target]["label"], order = int(edge_block[2].strip().split(" ")))

    return G
