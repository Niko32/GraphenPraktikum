#!/usr/bin/env python3

import sys
import csv
from bs4 import BeautifulSoup
import os

db_file_path = os.path.realpath(os.path.dirname(__file__))+'/metanetx'


reaction_xml = sys.argv[1]
outputsmiles = sys.argv[2] 

metaNetXError = {}
with open(db_file_path+'/db_corrections', mode='r') as db:
   for l in db:
      if l.startswith('MetaNetXError:'):   
           invalid = l.strip().split(';')
           id_name = invalid[0].split()
           metaNetXError[id_name[5]] = { 'name' : invalid[1], 'SMILES' : invalid[2]} 

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
      compound_id, _ = comp_deprecated[compound_id]

   name, smiles = compounds[compound_id]
   
   if smiles == "" and compound_id in metaNetXError:
         smiles = metaNetXError[compound_id]["SMILES"]
         name = metaNetXError[compound_id]["name"]
   
   if smiles.count('.') > 0:
      newname = name + " Subcomponent 1"
      for i in range(smiles.count('.')):
         newname += " + " + name + " Subcomponent "+str(i+2)
      name = newname

   return name, smiles


def parse_and_print_reaction(reaction_id):
   # 1 MNXM10@MNXD1 + 1 MNXM1312@MNXD1 + 2 MNXM1@MNXD1 = 1 MNXM1895@MNXD1 + 1 MNXM8@MNXD1 + 1 WATER@MNXD1

   name_formula = ""
   smiles_formula = ""

   reaction_string = reactions[reaction_id][0]
   if reaction_string == " = ":
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

      if not smiles:
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


invalid_reaction_count = 0;
reaction_count = 0

with open(outputsmiles, 'w') as omf:
 with open(reaction_xml, 'r') as f:
  for line in f:
     meta_reaction_id = line.strip()
     if meta_reaction_id == "":
         continue
   
     valid, name_formula, smiles_formula = parse_and_print_reaction(meta_reaction_id)
     if valid:
        print(file=omf)
        print("Bigg ID:", "-", "MetaNetXId:", meta_reaction_id, file=omf)
        print("ECs:", ";".join(reactions[meta_reaction_id][1]), file=omf)
        print(name_formula, file=omf)   
        print(smiles_formula, file=omf)
     else:
           invalid_reaction_count += 1
     
print("Unable to convert", invalid_reaction_count, "reactions of ", reaction_count, "total", file=sys.stderr)      




