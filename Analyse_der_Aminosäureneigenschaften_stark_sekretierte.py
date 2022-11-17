import numpy as np 
import random
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#stark sekretierte: alles über 60% von WT

df = pd.DataFrame ({'N-Terminus': ['MNICVNSLYRLSIPQFHSLYTEEVS', 'MTNITLSTQHYRIHRSDVEPVKEKT', 'MFNITNIQSTARHQSISNEASTEVP', 'MLISNVGINPAAYLNNHSVENSSQT', 'MPTGIKPIFINNMMSTYGLSHPHDS', 'MNISSSGINISTIPTQVKKSVETIR', 'MLKPICHSGSIKVPEYLETDKEKNA', 'MPFTFQIGNHSCQISERYLRDIIDN'],
                    'ID': ['SseL', 'SopE2', 'SlrP', 'SipC', 'GtgA', 'SboH', 'SopF', 'SteC'],
                    'description': ['pJA068', 'pJA049', 'pJA043', 'pJA041', 'pJA025', 'pJA035', 'pJA050', 'pJA073']})



AA2NA = {
    "A": list("nonpolar; ".split(",")), 
    "R": list("positive; ".split(",")), 
    "N": list("polar; ".split(",")), 
    "D": list("negative; ".split(",")), 
    "C": list("polar; ".split(",")), 
    "Q": list("polar; ".split(",")), 
    "E": list("negative; ".split(",")), 
    "G": list("polar; ".split(",")), 
    "H": list("positive; ".split(",")), 
    "I": list("nonpolar; ".split(",")), 
    "L": list("nonpolar; ".split(",")), 
    "K": list("positive; ".split(",")), 
    "M": list("nonpolar; ".split(",")), 
    "F": list("nonpolar; ".split(",")), 
    "P": list("nonpolar; ".split(",")), 
    "S": list("polar; ".split(",")), 
    "T": list("polar; ".split(",")), 
    "W": list("nonpolar; ".split(",")), 
    "Y": list("polar; ".split(",")), 
    "V": list("nonpolar; ".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["------"])) for c in seq]
    return "".join(na_seq)
      
i=0
n = len(df)
all_as_charges = []
all_ids = []
all_desc = []

while i < n:
    as_charge = aa2na(df['N-Terminus'][i])
    all_ids = np.append (all_ids, df ['ID'][i])
    all_desc = np.append (all_desc, df ['description'] [i])
    all_as_charges = np.append(all_as_charges, as_charge)
    i=i+1

df_all_ids = pd.DataFrame (all_ids, columns = ['ID'])
df_charges = pd.DataFrame(all_as_charges, columns = ['charge'])
df_charges_id = pd.merge(df_all_ids, df_charges, left_index=True, right_index=True)
df_charges = df_charges['charge']

df_charges_per_residue = df_charges.str.split(';', expand=True)
df_charges_per_residue.columns = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25']
print(df_charges_per_residue['0'].value_counts())
print(df_charges_per_residue['1'].value_counts())
print(df_charges_per_residue['2'].value_counts())
print(df_charges_per_residue['3'].value_counts())
print(df_charges_per_residue['4'].value_counts())
print(df_charges_per_residue['5'].value_counts())
print(df_charges_per_residue['6'].value_counts())
print(df_charges_per_residue['7'].value_counts())
print(df_charges_per_residue['8'].value_counts())
print(df_charges_per_residue['9'].value_counts())
print(df_charges_per_residue['10'].value_counts())
print(df_charges_per_residue['11'].value_counts())
print(df_charges_per_residue['12'].value_counts())
print(df_charges_per_residue['13'].value_counts())
print(df_charges_per_residue['14'].value_counts())
print(df_charges_per_residue['15'].value_counts())
print(df_charges_per_residue['16'].value_counts())
print(df_charges_per_residue['17'].value_counts())
print(df_charges_per_residue['18'].value_counts())
print(df_charges_per_residue['19'].value_counts())
print(df_charges_per_residue['20'].value_counts())
print(df_charges_per_residue['21'].value_counts())
print(df_charges_per_residue['22'].value_counts())
print(df_charges_per_residue['23'].value_counts())
print(df_charges_per_residue['24'].value_counts())

print (df_charges_per_residue)

AA2NA = {
    "A": list("aliphatic; ".split(",")), 
    "R": list("basic; ".split(",")), 
    "N": list("amidated; ".split(",")), 
    "D": list("acidic; ".split(",")), 
    "C": list("sulphurous; ".split(",")), 
    "Q": list("amidated; ".split(",")), 
    "E": list("acidic; ".split(",")), 
    "G": list("aliphatic; ".split(",")), 
    "H": list("basic; ".split(",")), 
    "I": list("aliphatic; ".split(",")), 
    "L": list("aliphatic; ".split(",")), 
    "K": list("basic; ".split(",")), 
    "M": list("aliphatic_sulphurous; ".split(",")), 
    "F": list("aromatic; ".split(",")), 
    "P": list("aliphatic; ".split(",")), 
    "S": list("hydroxylated; ".split(",")), 
    "T": list("hydroxylated; ".split(",")), 
    "W": list("aromatic; ".split(",")), 
    "Y": list("aromatic_hydroxylated; ".split(",")), 
    "V": list("aliphatic; ".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["------"])) for c in seq]
    return "".join(na_seq)
      
i=0
n = len(df)
all_as_properties = []
all_ids = []
all_desc = []

while i < n:
    as_properties = aa2na(df['N-Terminus'][i])
    all_ids = np.append (all_ids, df['ID'][i])
    all_desc = np.append (all_desc, df['description'] [i])
    all_as_properties = np.append(all_as_properties, as_properties)
    i=i+1

df_all_ids = pd.DataFrame (all_ids, columns = ['ID'])
df_properties = pd.DataFrame(all_as_properties, columns = ['property'])
df_properties_id = pd.merge(df_all_ids, df_properties, left_index=True, right_index=True)
df_properties = df_properties['property']

df_properties_per_residue = df_properties.str.split(';', expand=True)
df_properties_per_residue.columns = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25']
print(df_properties_per_residue['0'].value_counts())
print(df_properties_per_residue['1'].value_counts())
print(df_properties_per_residue['2'].value_counts())
print(df_properties_per_residue['3'].value_counts())
print(df_properties_per_residue['4'].value_counts())
print(df_properties_per_residue['5'].value_counts())
print(df_properties_per_residue['6'].value_counts())
print(df_properties_per_residue['7'].value_counts())
print(df_properties_per_residue['8'].value_counts())
print(df_properties_per_residue['9'].value_counts())
print(df_properties_per_residue['10'].value_counts())
print(df_properties_per_residue['11'].value_counts())
print(df_properties_per_residue['12'].value_counts())
print(df_properties_per_residue['13'].value_counts())
print(df_properties_per_residue['14'].value_counts())
print(df_properties_per_residue['15'].value_counts())
print(df_properties_per_residue['16'].value_counts())
print(df_properties_per_residue['17'].value_counts())
print(df_properties_per_residue['18'].value_counts())
print(df_properties_per_residue['19'].value_counts())
print(df_properties_per_residue['20'].value_counts())
print(df_properties_per_residue['21'].value_counts())
print(df_properties_per_residue['22'].value_counts())
print(df_properties_per_residue['23'].value_counts())
print(df_properties_per_residue['24'].value_counts())

print (df_properties_per_residue)