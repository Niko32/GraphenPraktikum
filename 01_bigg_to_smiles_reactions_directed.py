#!/usr/bin/env python3

import sys
import csv
from bs4 import BeautifulSoup
import os
import re

from rdkit import Chem
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

db_file_path = os.path.realpath(os.path.dirname(__file__))+'/metanetx'


reaction_xml = sys.argv[1]
outputsmiles = sys.argv[2] 
path_database = sys.argv[3]

# Load Database for missing SMILES
metaNetXError = {}
idList_MissMeta = []
with open(path_database, mode='r') as db:
    for j in db:
        if re.search('MetaNetXError:',j):   
           invalid = j.split(';')
           id_name = invalid[0].split()
           idList_MissMeta.append(id_name[5])
           print(invalid)
           metaNetXError[id_name[5]] = { 'Networkname' : invalid[1], 'SMILES' : invalid[2]}  
#print(metaNetXError)

comp_deprecated = {}
with open(db_file_path+'/chem_depr.tsv', mode='r') as tsv:
   for line in tsv:
      if line.startswith("#"):
         continue
       
      id_old, id_new, version = line.strip().split('\t')
      version = int(version.replace('*', '0').replace('.',''))
     
      if id_old in comp_deprecated and version < comp_deprecated[id_old][1]:
         continue

      comp_deprecated[id_old] = (id_new, version)


compounds = {}
#ID	name	reference	formula	charge	mass	InChI	InChIKey	SMILES
with open(db_file_path+'/chem_prop.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      meta_id = splits[0]
      name = splits[1]
      smiles = splits[8]
      inChI = splits[6]    
      compounds[meta_id] = (name, smiles)

reactions = {}
#ID	mnx_equation	reference	classifs	is_balanced	is_transport
with open(db_file_path+'/reac_prop.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      meta_id = splits[0]
      mnx_equation = splits[1]
      ecs = [ec for ec in splits[3].split(";") if not ec == ''] 

      reactions[meta_id] = (mnx_equation, ecs)

checkref = {}
with open(db_file_path+'/reac_xref.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      if ":" not in splits[0]:
         continue

      source, ref_id = splits[0].split(":")
      meta_id = splits[1]

      if source == "bigg.reaction":
         checkref[ref_id] = meta_id


print("# Compounds in DB", len(compounds), file=sys.stderr)
print("# Reaction in DB", len(reactions), file=sys.stderr)

def get_compound_info(compound_id):
   while compound_id in comp_deprecated:
      compound_id = comp_deprecated[comp_deprecated]

   return compounds[compound_id]

def parse_and_print_reaction(reaction_id):
   # 1 MNXM10@MNXD1 + 1 MNXM1312@MNXD1 + 2 MNXM1@MNXD1 = 1 MNXM1895@MNXD1 + 1 MNXM8@MNXD1 + 1 WATER@MNXD1

   name_formula = ""
   smiles_formula = ""

   reaction_string = reactions[reaction_id][0]
   if reaction_string == " = ":
      print("MetaNetXError: Empty reaction in MetaNetX", reaction_id, file=sys.stderr)
      return (False, "", "")

   left_reaction, right_reaction = reaction_string.split('=')
   left_reac_comps = left_reaction.split('+')
   right_reac_comps = right_reaction.split('+')

   first = True 
   for lr in left_reac_comps:
      count_string, comp_id_str = lr.strip().split()
      count = int(count_string.strip())
      comp_id = comp_id_str.strip().split('@')[0]
      name, smiles = get_compound_info(comp_id)

      if not smiles:
         if comp_id in metaNetXError:
            print('COORECTION POSSIBILE' + comp_id)
            entry =  metaNetXError[comp_id]
            smiles = entry['SMILES']
         else:    
             print("MetaNetXError: No Smiles for Compound", comp_id, name, "InChI", inChI, file=sys.stderr)
             return (False, "", "")

      for i in range(count):
         if not first:
            name_formula += " + "  
         name_formula += name
      
         if not first:
            smiles_formula += "."
            
         first = False
         smiles_formula += smiles
      
   name_formula += " = "   
   smiles_formula += ">>"
   
   first = True 
   for rr in right_reac_comps:
      count_string, comp_id_str = rr.strip().split()
      count = int(count_string.strip())
      comp_id = comp_id_str.strip().split('@')[0]
      name, smiles = get_compound_info(comp_id)
      #name, smiles, inChI = get_compound_info(comp_id)
      
      if not smiles:
         if comp_id in metaNetXError:
            print('COORECTION POSSIBILE: ' + comp_id)
            entry =  metaNetXError[comp_id]
            smiles = entry['SMILES']
         else:   
            print("MetaNetXError: No Smiles for Compound", comp_id, name, "InChI", inChI, file=sys.stderr)
            return (False, "", "")

      for i in range(count):
         if not first:
            name_formula += " + "  
         name_formula += name
      
         if not first:
            smiles_formula += "."
            
         first = False
         smiles_formula += smiles
      
   return (True, name_formula, smiles_formula)

def count_atoms(mol, counter):
   molWithH = Chem.AddHs(mol)
   
   for atom in molWithH.GetAtoms():
       sym = atom.GetSymbol()
       counter[sym] = counter.get(sym, 0) + 1

def is_balanced(smiles_formula):
    smiles_sides = smiles_formula.split('>>')

    left_mols = [ Chem.MolFromSmiles(s) for s in smiles_sides[0].split('.')  ]
    right_mols = [ Chem.MolFromSmiles(s) for s in smiles_sides[1].split('.')  ]       

    counter_left = {}
    for m in left_mols:
        count_atoms(m, counter_left)

    counter_right = {}
    for m in right_mols:
        count_atoms(m, counter_right)
    
    diff = {}
    if not counter_left == counter_right:
        diff = {k: counter_left.get(k,0) - counter_right.get(k,0) for k in counter_left.keys() | counter_right.keys() if counter_left.get(k,0) - counter_right.get(k,0) != 0}
       
    return (counter_left == counter_right, diff)
    
with open(reaction_xml, 'r') as f:
   reaction_data = f.read()

invalid_reaction_count = 0;
big_missing_count = 0
reaction_count = 0
unbalanced_count = 0

with open(outputsmiles, 'w') as omf:
  reaction_parser = BeautifulSoup(reaction_data, "xml")
  for xml_react in reaction_parser.find_all('reaction'):
     reaction_count += 1
     bigg_id = xml_react.get('id')
     reverse = xml_react.get('reversible')
     if bigg_id in checkref:
        meta_reaction_id = checkref[bigg_id]   
        valid, name_formula, smiles_formula = parse_and_print_reaction(meta_reaction_id)
        
        #Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(smiles_formula)),allHsExplicit=False, canonical=False)
        
        if valid:
           print(file=omf)
           print("Bigg ID:", bigg_id, "MetaNetXId:", meta_reaction_id, file=omf)
           print("ECs:", ";".join(reactions[meta_reaction_id][1]), file=omf)
           print(name_formula, file=omf)   
           print(smiles_formula, file=omf)
           print("reversibile:",reverse, file=omf)
        else:
           invalid_reaction_count += 1
     else:
     
        invalid_reaction_count += 1

print("Unable to parse", big_missing_count, "and unable to convert", invalid_reaction_count, "reactions of ", reaction_count, "total", file=sys.stderr)      
print("# of Unbalanced Reactions:", unbalanced_count)





