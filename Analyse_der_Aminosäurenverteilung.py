import numpy as np 
import random
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

with open("backbone.txt", "r", encoding='utf-8-sig') as plasmid:
    plasmid_str = plasmid.read()

full_dataset = '202204_T3SS_secreted_proteins_all_validated.fasta'

full_dataset_id = []
full_dataset_name = []
full_dataset_seq = []
full_dataset_desc = []

for record in SeqIO.parse(full_dataset, 'fasta'):
    full_dataset_id.append(record.id)
    full_dataset_seq.append(str(record.seq))
    full_dataset_desc.append(record.description)
    full_dataset_name.append(record.name)

full_dataset_dict = {'ID':full_dataset_id,'sequence':full_dataset_seq,'description':full_dataset_desc,'name':full_dataset_name}


df_full_dataset = pd.DataFrame.from_dict(full_dataset_dict)
df_salmonella_full_dataset = df_full_dataset[df_full_dataset['description'].str.contains("Salmonella")]

df_short_dataset= df_full_dataset.drop_duplicates(subset=['sequence'], keep='last',ignore_index=True)
df_salmonella_short_dataset = df_short_dataset[df_short_dataset['description'].str.contains("Salmonella")]
df_salmonella_short_dataset = df_salmonella_short_dataset.reset_index(drop=True)
df_salmonella_short_dataset.columns = ['ID', 'sequence', 'description', 'name']
df_salmonella_seq = df_salmonella_short_dataset['sequence']
df_salmonella_seq = df_salmonella_seq.values
df_description0 = df_salmonella_short_dataset['description']
df_id1 = df_salmonella_short_dataset ['ID']
df_id_trim = df_id1.str.split("|", expand=True)
df_id_trim.columns = ['DB', 'ID', 'ID2']
df_id=df_id_trim['ID']
df_id = df_id.values

df_desc = df_description0.str.split("|", expand=True)
df_desc.columns = ['DB', 'id', 'name', 'function', 'target', 'evidence', 'localization', 'host', 'PMID']
df_description2 = df_desc['name'].str.split("OS=", expand=True)
df_description2.columns = ['name', 'description']
df_description = df_description2 ['name']
df_species = df_description2['description']
df_species.columns = ['description']

def trim_ntermini (dna_seq):
    trim_seq = dna_seq[:25]
    return "".join(trim_seq) 

i=0
n = len(df_salmonella_seq)
all_ntermini_seq = []

while i < n:
    trimmed_ntermini = trim_ntermini(df_salmonella_seq[i])
    all_ntermini_seq = np.append (all_ntermini_seq, trimmed_ntermini)
    i=i+1
    
df_ntermini = pd.DataFrame(all_ntermini_seq, columns = ['ntermini'])
df_ntermini_id_description= pd.merge(df_ntermini, df_salmonella_short_dataset, left_index=True, right_index=True)

salmonella_ntermini_list = []

for i in range(len(df_ntermini)):
    entry = SeqRecord(Seq(df_ntermini['ntermini'].iloc[i]),
    id=df_salmonella_short_dataset['ID'].iloc[i],
    description=df_description.iloc[i])
    salmonella_ntermini_list.append(entry)

AA2NA = {
    "A": list("A; ".split(",")), 
    "R": list("R; ".split(",")), 
    "N": list("N; ".split(",")), 
    "D": list("D; ".split(",")), 
    "C": list("C; ".split(",")), 
    "Q": list("Q; ".split(",")), 
    "E": list("E; ".split(",")), 
    "G": list("G; ".split(",")), 
    "H": list("H; ".split(",")), 
    "I": list("I; ".split(",")), 
    "L": list("L; ".split(",")), 
    "K": list("K; ".split(",")), 
    "M": list("M; ".split(",")), 
    "F": list("F; ".split(",")), 
    "P": list("P; ".split(",")), 
    "S": list("S; ".split(",")), 
    "T": list("T; ".split(",")), 
    "W": list("W; ".split(",")), 
    "Y": list("Y; ".split(",")), 
    "V": list("V; ".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["------"])) for c in seq]
    return "".join(na_seq)
      
i=0
n = len(df_ntermini_id_description)
all_as_properties = []
all_ids = []
all_desc = []


while i < n:
    as_properties = aa2na(df_ntermini_id_description['ntermini'][i])
    all_ids = np.append (all_ids, df_ntermini_id_description ['ID'][i])
    all_desc = np.append (all_desc, df_ntermini_id_description ['description'] [i])
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

print(df_properties_per_residue [df_properties_per_residue == " A"].count())
print(df_properties_per_residue [df_properties_per_residue == " R"].count())
print(df_properties_per_residue [df_properties_per_residue == " N"].count())
print(df_properties_per_residue [df_properties_per_residue == " D"].count())
print(df_properties_per_residue [df_properties_per_residue == " C"].count())
print(df_properties_per_residue [df_properties_per_residue == " Q"].count())
print(df_properties_per_residue [df_properties_per_residue == " E"].count())
print(df_properties_per_residue [df_properties_per_residue == " G"].count())
print(df_properties_per_residue [df_properties_per_residue == " H"].count())
print(df_properties_per_residue [df_properties_per_residue == " I"].count())
print(df_properties_per_residue [df_properties_per_residue == " L"].count())
print(df_properties_per_residue [df_properties_per_residue == " K"].count())
print(df_properties_per_residue [df_properties_per_residue == "M"].count())
print(df_properties_per_residue [df_properties_per_residue == " M"].count())
print(df_properties_per_residue [df_properties_per_residue == " F"].count())
print(df_properties_per_residue [df_properties_per_residue == " P"].count())
print(df_properties_per_residue [df_properties_per_residue == " S"].count())
print(df_properties_per_residue [df_properties_per_residue == " T"].count())
print(df_properties_per_residue [df_properties_per_residue == " W"].count())
print(df_properties_per_residue [df_properties_per_residue == " Y"].count())
print(df_properties_per_residue [df_properties_per_residue == " V"].count())

