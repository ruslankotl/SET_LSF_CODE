'''
This module is used to process the data to make it the standardized csv file,
ready for loading into the GNN module.
'''
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem

def init_args():
     '''
     Setting up the input arguments.
     '''
     parser = argparse.ArgumentParser()
     parser.add_argument('--data_path', '-d', type=str)
     parser.add_argument('--save_path', '-s', type=str)
     return parser.parse_args()

def parent_centres(reactant_smiles:str, product_smiles:str)->np.array(list):
    reactant = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    centres = list()

    match = product.GetSubstructMatch(reactant)

    # get the reactive atoms indices on product: in substructure, but with neighbours outside MCS
    
    reaction_sites = [atom.GetIdx() for atom in product.GetAtoms() \
         if (atom.GetIdx() in match) and \
         not all((n.GetIdx() in match) for n in atom.GetNeighbors())]

    # noting that reactant atoms have a through numbering 
    for reactant_index, prod_index in enumerate(match):
          if prod_index in reaction_sites:
               centres.append(reactant_index)

    return np.array(centres)

def main():
    args = init_args()
    data_path = args.data_path
    save_path = args.save_path
    df = pd.read_excel(data_path)

    reaction_sites = []

    for i in range(len(df)):
        # Canonicalizing the SMILES strings - must be done for consistent mapping!
        sm = Chem.CanonSmiles(df.loc[i, 'reactant'])
        if pd.isna(df.loc[i, 'product']) == True or df.loc[i, 'product'] == 'NO STRUCTURE':
            p = Chem.CanonSmiles(df.loc[i, 'reactant'])
        else:
            p = Chem.CanonSmiles(df.loc[i, 'product'])
        reaction_site = parent_centres(sm, p)
        reaction_sites.append(reaction_site)

    df['parent_centres'] = reaction_sites
    df.to_pickle(save_path)

if __name__ == '__main__':
    main()