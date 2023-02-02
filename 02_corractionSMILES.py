#!/usr/bin/env python3

import sys
import re 
import itertools

inputFileName = sys.argv[1]
path_database = sys.argv[2]
path_error = sys.argv[3]
corrMissingBiggId = sys.argv[4]
transportRXN = sys.argv[5]
outputFileName = sys.argv[6]

def isa_group_separator(line):
    return line=='\n'

def transporter(reaction):
    try:
        reaction = reaction.split('=')
        left_side = reaction[0]
        left_side = left_side.strip()
        right_side = reaction[1]
        right_side = right_side.strip()
        if left_side == right_side:
            return True
        else:
            return False
    except:
        print('-- Not a adequate Reaction --') 

# Read Database    
invalidBiGG = {}
unbalanced = {}
with open(path_database, mode='r') as db:
    for j in db:
        if re.search('InvalidBiGG:',j):   
           invalid = j.split(';')
           id_invalid = invalid[0].split()

           invalidBiGG[id_invalid[1]] =  {'Reaction:' : invalid[1], 'Metabolites' : invalid[3], 'corrSMILES' : invalid[4], 'reversible' : invalid[5]  } 
        elif re.search('Unbalanced:',j):   
           invalid = j.split(';')
           id_invalid = invalid[0].split()
           #print(invalid)
           unbalanced[id_invalid[1]] = { 'MetanetXID': id_invalid[2], 'CorrSMILES' : invalid[1]}  

# Read Error List
idList_invalidBiGG = []
idList_unbalanced = []

with open(path_error, mode='r') as er:
    for j in er:
        if re.search('InvalidBiGG:',j):   
           invalid = j.split()
           idList_invalidBiGG.append(invalid[1])
        elif re.search('Unbalanced:',j):   
           invalid = j.split()
           idList_unbalanced.append(invalid[1])    
           
           

# correcting unbalancing Reaktion and add lost Reactions
with open(inputFileName, mode='r') as inp:
    with open(outputFileName, mode='w') as out: 
      for key,group in itertools.groupby(inp,isa_group_separator):
          if key:
             continue
          x = list(group)
          ids = x[0]
          ec = x[1]
          name_reaction = x[2]
          smiles_reaction = x[3]
          reverse = x[4]
          reverse = reverse.replace('\n','')

          # Delete Transport Reactions
          if transportRXN == 'False':
            print(name_reaction)
            if transporter(name_reaction):
                continue
          # correced unbalanced reaktion
          bigg_id = ids.split()
          if bigg_id[1] in idList_unbalanced:
            smiles_reaction = unbalanced['CorrSMILES']
          #print reaction
          out.write(ids)
          out.write(ec)
          out.write(name_reaction)
          out.write(smiles_reaction)
          out.write(reverse +'\n')
      # add Missing Reaction:             
      for i in idList_invalidBiGG:
        entry =  invalidBiGG[i]
        print(entry['Metabolites'])
        if transportRXN == 'False':
            if transporter(entry['Metabolites']):
                continue
        out.write('Bigg ID: '+ i +'\n')
        out.write('ECs: \n')
        out.write(entry['Metabolites'] + '\n')
        out.write(entry['corrSMILES'] + '\n') 
        out.write(entry['reversible'] + '\n') 
            
# correcting unbalancing Reaktion and add lost Reactions
with open(inputFileName, mode='r') as inp:
    with open(outputFileName + '_aam', mode='w') as out: 
      for key,group in itertools.groupby(inp,isa_group_separator):
          if key:
             continue
          x = list(group)
          name_reaction = x[2]
          smiles_reaction = x[3]
          reverse = x[4]
          #reverse = reverse.replace('\n','')
          #print reaction
          out.write(name_reaction)
          out.write(smiles_reaction)
          out.write(reverse)

           




