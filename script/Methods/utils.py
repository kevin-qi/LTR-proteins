import pandas as pd
import numpy as np
from Methods.scales import stScales

def AASeq(mutant_name, wildtype_sequence):
    """
        ARGS: rluc8_dataset
        RETURNS: standardized AA sequence relevant AA in wildtype rluc8
    """
    wildtype_sequence = wildtype_sequence.upper()
    
    seq = wildtype_sequence
    
    for mutation in mutant_name.split("/")[1:]:
        index = int(mutation[1:-1])
        original = mutation[0]
        new = mutation[-1]
        
        if(wildtype_sequence[index-1] != original):
            print("Alignment Error! Wild: {}, Original: {}".format(wildtype_sequence[index], original) )
            return
        
        seq = seq[:index-1] + new + seq[index:]
        if(seq[index-1] != new):
            print("Alignment Error: ")
            return
    return seq

def stEncode(seq):
    encoded = []
    for AA in seq:
        encoded += list(stScales(AA))
    return np.array(encoded)

def extractRelevantResidues(dataset):
    """
        Args: rluc8_dataset
        Returns: 1) Array of all relevant residue addresses without duplicates, ordered small to large.
                 2) Dictionary mapping of residue address to index
    """
    variants = np.array(dataset["mutant"])
    residue_set = set()
    
    # Extract set of residue addresses into a set
    for mutant in variants[1:]:
        mutations = mutant.split("/")[1:]
        for mutation in mutations:
            residue = ''.join(filter(lambda x: x.isdigit(), mutation))
            residue_set.add(residue)
    
    # Order small to large
    residue_array = list(residue_set)
    residue_array.sort(key = lambda x: int(x))
    
    # Create index mapping
    index_map = {}
    inverse_map = {}
    i = 0
    for residue in residue_array:
        inverse_map[i] = str(residue)
        index_map[str(residue)] = i
        i+=1
    return residue_array, index_map, inverse_map

def extractWildtypeSequence(dataset):
    """
        ARGS: rluc8_dataset
        RETURNS: standardized AA sequence relevant AA in wildtype rluc8
    """
    residue_array, address2index, index2address = extractRelevantResidues(dataset)
    wildtype_sequence = "Z"*50
    
    variants = np.array(dataset["mutant"])
    for mutant in variants:
        mutations = mutant.split("/")[1:]
        for mutation in mutations:
            address = ''.join(filter(lambda x: x.isdigit(), mutation)) # eg) 223 from I223C
            index = address2index[address] # eg) 0 from 11
            
            wildtype_sequence = wildtype_sequence[:index] + mutation[0] + wildtype_sequence[index+1:]
    return wildtype_sequence

def relevantAASeq(mutant_name, dataset):
    """
        ARGS: 1) name of mutant, eg) RLuc8/A123S/D162E/I163V/F262W
              2) rluc8_dataset
        RETURNS: standardized sequence of mutant
    """
    residue_array, address2index, index2address = extractRelevantResidues(dataset)
    wildtype_seq = extractWildtypeSequence(dataset)
    mutant_sequence = wildtype_seq
    for mutation in mutant_name.split("/")[1:]:
        address = int(''.join(filter(lambda x: x.isdigit(), mutation))) # eg) 223 from I223C
        index = address2index[str(address)] # eg) 0 from 11
        mutant_sequence = mutant_sequence[:index] + str(mutation[-1]) + mutant_sequence[index+1:]
    return mutant_sequence

def encode(df, wildtype_sequence):
    df = df.copy()
    df["X"] = df["mutant"].apply(lambda x: AASeq(x, wildtype_sequence))
    df["X"] = df["X"].apply(lambda x: stEncode(x))
    return df

def generatePairs(df):
    df = df.merge(df, on=df.assign(key_col=1)['key_col'], suffixes=('_1', '_2')).reset_index(drop=True)
    df = df.drop(["key_0"], axis=1)
    df = df.loc[df["mutant_1"] != df["mutant_2"]]
    df = df.loc[df["label_1"] > df["label_2"]]
    df["label"] = (df["label_1"] > df["label_2"])*1.
    return df

def contactMap(threshold):
    import Bio.PDB as pdb
    pdb_code = "2pse"
    pdb_filename = "Data/2pse.pdb"
    
    def residue_dist(res1, res2):
        diff = res1["CA"].coord - res2["CA"].coord
        return np.sqrt(np.sum(diff*diff))
    
    rluc = pdb.PDBParser().get_structure(pdb_code, pdb_filename)[0]
    residue_array, address2index, index2address = extractRelevantResidues(pd.read_csv("Data/data.txt"))

    contacts = np.zeros((305,305))
    residues = []
    for row, res in enumerate(rluc["A"]):
        residues += [(row,res)]
    for res1 in residues[:305]:
        for res2 in residues[:305]:
            contacts[res1[0]][res2[0]] = residue_dist(res1[1], res2[1])
    return contacts<threshold

def makeMask(contacts):
    mask = []
    for row in contacts:
        dup = list(row)*8
        mask += list(np.array(dup).reshape(8,len(row)))
    return np.array(mask)