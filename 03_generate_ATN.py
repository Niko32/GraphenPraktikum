#!/usr/bin/env python3

import re
import os
import sys
import enum

import logging

from custom_pysmiles import read_smiles
import networkx as nx
from networkx.algorithms import isomorphism as nxisomorphism

from pyvis.network import Network

logging.basicConfig(level=logging.INFO)

mappedsmiles = sys.argv[1]
outputgml = sys.argv[2] 

## ========== DEF BLOCK

@enum.unique
class TransitionType(enum.IntEnum):
    """Possible SMILES token types"""
    NO_TRANSITION = 0
    SYMMETRY = 1
    REACTION = 2
   
def findIsomorphATNStructure(ATN, mol):

    em = nxisomorphism.categorical_edge_match(['order'],[0])
    nm = nxisomorphism.categorical_node_match(['element', 'isotope', 'hcount', 'charge'],['', 0, 0, 0])
    GM = nxisomorphism.GraphMatcher(ATN, mol, node_match=nm, edge_match=em)   
    result_iso = GM.subgraph_is_monomorphic()
    mapping_ATN_to_mol = GM.mapping   

    return (result_iso, mapping_ATN_to_mol)

def addAutomorphisms(mol, limit_to_orbits=True):
    em = nxisomorphism.categorical_edge_match(['order'],[0])
    nm = nxisomorphism.categorical_node_match(['element', 'isotope', 'hcount', 'charge'],['', 0, 0, 0])
    GM = nxisomorphism.GraphMatcher(mol, mol, node_match=nm, edge_match=em)   

    if limit_to_orbits:
        blockset = set()
        for node in mol.nodes():
            if node in blockset:
                continue
            for isomorphism in GM.isomorphisms_iter():
                mol.add_edge(node, isomorphism[node])
                mol.edges[node, isomorphism[node]]['transistion'] = TransitionType.SYMMETRY
                blockset.add(isomorphism[node])
    else:
        for isomorphism in GM.isomorphisms_iter():
            for i in isomorphism:
                if i != isomorphism[i] and not mol.has_edge(i, isomorphism[i]):
                    mol.add_edge(i, isomorphism[i])
                    mol.edges[i, isomorphism[i]]['transistion'] = TransitionType.SYMMETRY
        
def parseXDuct(name, smiles, hydro, compound_to_subgraph, compoundId_to_compound, ATN, mapped_atoms):

    logging.debug("Parse " + name + " : " + smiles)

    mol = read_smiles(smiles) # read smile
    hydromol = read_smiles(hydro)
    
    for n, h in zip(mol.nodes(), hydromol.nodes()):
       if 'hcount' in hydromol.nodes[h]:
          mol.nodes[n]['hcount'] = hydromol.nodes[h]['hcount']

    if name in compound_to_subgraph:

        logging.debug("Use Existing " + name)

        # we already added this compound to the network
        compound_subgraph = ATN.subgraph(compound_to_subgraph[name])
        has_isomorph_subgraph, mapping_ATN_to_mol = findIsomorphATNStructure(compound_subgraph, mol)
        if has_isomorph_subgraph:
            for atn_node in mapping_ATN_to_mol:
                if 'class' in mol.nodes[mapping_ATN_to_mol[atn_node]]:
                    mapped_atoms[ mol.nodes[mapping_ATN_to_mol[atn_node]]['class'] ] = atn_node
        else:
            logging.error("Compound Naming Error: " + name + " : " + smiles)
    else:

        logging.debug("Add New Compound "+ name)

        # this is definitely new
        nextCId = len(compound_to_subgraph)
        compoundId_to_compound[nextCId] = name

        rename = {node : str(nextCId)+'_'+str(node) for node in mol.nodes()} # rename all nodes so that we cannot have collisions in the ATN
        nx.relabel_nodes(mol, rename, copy=False)
  
        for e in mol.edges():
            mol.edges[e]['transistion'] = TransitionType.NO_TRANSITION

        addAutomorphisms(mol)

        ATN.add_nodes_from(mol.nodes(data=True))
        ATN.add_edges_from(mol.edges(data=True))
        for mol_node, data in mol.nodes(data=True):
            if 'class' in data:
                mapped_atoms[data['class']] = mol_node
            compound_to_subgraph.setdefault(name, []).append(mol_node)
        
# ======== MAIN

# first an intermediary with bonding and transition edges that is trimmed afterwards
ATN=nx.Graph()
compound_to_subgraph = {}
compoundId_to_compound = {}

reactions = []

# load all SMILES from txt files
with open( mappedsmiles , 'r') as smiles_file:
  while True:    
    name_line = smiles_file.readline().strip()
    hydrogen_smiles_line = smiles_file.readline().strip()
    mapped_smiles_line = smiles_file.readline().strip()
    if not name_line or not  mapped_smiles_line:
       break;

    logging.info("Next Reaction ==============")
    logging.info(name_line.strip())
    logging.info(mapped_smiles_line.strip())

    smiles_str = mapped_smiles_line.strip().replace("@", '').replace("/", '')
    hydrogen_str = hydrogen_smiles_line.strip()
    names_str = name_line.strip()
    #filename = ''.join(letter for letter in name_line if letter.isalnum())

    reactions.append( (names_str, smiles_str) )

    # create left and right list with the smiles
    smiles_sides = smiles_str.split('>>')
    smiles_left = smiles_sides[0].split('.')
    smiles_right = smiles_sides[1].split('.')
    
     # create left and right list with the hydrogen smiles
    hydrogen_sides = hydrogen_str.split('>>')
    hydrogen_left = hydrogen_sides[0].split('.')
    hydrogen_right = hydrogen_sides[1].split('.')

    # create left and right list with the metabolite names
    names_sides = names_str.split('=')
    names_left = map(str.strip, names_sides[0].split(' + '))
    names_right = map(str.strip, names_sides[1].split(' + '))

    logging.debug("Parse Educts")
    mapped_educt_atoms = {}
    for name, smiles, hydro in zip(names_left, smiles_left, hydrogen_left):
        parseXDuct(name, smiles, hydro, compound_to_subgraph, compoundId_to_compound, ATN, mapped_educt_atoms)

    logging.debug("Parse Products")
    mapped_product_atoms = {}
    for name, smiles, hydro in zip(names_right, smiles_right, hydrogen_right):
        parseXDuct(name, smiles, hydro, compound_to_subgraph, compoundId_to_compound, ATN, mapped_product_atoms)

    for c in mapped_educt_atoms:
        n1 = mapped_educt_atoms[c]
        n2 = mapped_product_atoms[c]
        ATN.add_edge(n1, n2)
        ATN.edges[n1, n2]['transistion'] = TransitionType.REACTION
        ATN.edges[n1, n2].setdefault('reaction_id', []).append(len(reactions)-1)

nx.write_gml(ATN, outputgml)

# DRAWING
"""
draw = nx.Graph()
draw.add_nodes_from(ATN.nodes())
draw.add_edges_from(ATN.edges())

for n in ATN.nodes():
    draw.nodes[n]['label'] = ATN.nodes[n]['element']+str(ATN.nodes[n]['hcount'])+'H'
    draw.nodes[n]['color'] = "blue"

for e in ATN.edges():
    if ATN.edges[e]['transistion'] == TransitionType.SYMMETRY:
        draw.edges[e]['color'] = "green"
    elif ATN.edges[e]['transistion'] == TransitionType.REACTION:
        draw.edges[e]['color'] = "red"
    else:
        draw.edges[e]['label'] = str(ATN.edges[e]['order'])
        draw.edges[e]['color'] = "blue"

nt = Network('1000px', '1000px')
nt.show_buttons()
nt.from_nx(draw)
nt.show('nx2.html')
"""
