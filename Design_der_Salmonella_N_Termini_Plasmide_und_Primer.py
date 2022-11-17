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
    
print (all_ntermini_seq) 

df_ntermini = pd.DataFrame(all_ntermini_seq, columns = ['ntermini'])
df_ntermini_id_description= pd.merge(df_ntermini, df_salmonella_short_dataset, left_index=True, right_index=True)

salmonella_ntermini_list = []

for i in range(len(df_ntermini)):
    entry = SeqRecord(Seq(df_ntermini['ntermini'].iloc[i]),
    id=df_salmonella_short_dataset['ID'].iloc[i],
    description=df_description.iloc[i])
    salmonella_ntermini_list.append(entry)


AA2NA = {
    "A": list("GCT,GCC,GCA,GCG".split(",")),
    "R": list("CGT,CGC,CGA,CGG,AGA,AGG".split(",")),
    "N": list("AAT,AAC".split(",")),
    "D": list("GAT,GAC".split(",")),
    "C": list("TGT,TGC".split(",")),
    "Q": list("CAA,CAG".split(",")),
    "E": list("GAA,GAG".split(",")),
    "G": list("GGT,GGC,GGA,GGG".split(",")),
    "H": list("CAT,CAC".split(",")),
    "I": list("ATT,ATC,ATA".split(",")),
    "L": list("TTA,TTG,CTT,CTC,CTA,CTG".split(",")),
    "K": list("AAA,AAG".split(",")),
    "M": list("ATG".split(",")),
    "F": list("TTT,TTC".split(",")),
    "P": list("CCT,CCC,CCA,CCG".split(",")),
    "S": list("TCT,TCC,TCA,TCG,AGT,AGC".split(",")),
    "T": list("ACT,ACC,ACA,ACG".split(",")),
    "W": list("TGG".split(",")),
    "Y": list("TAT,TAC".split(",")),
    "V": list("GTT,GTC,GTA,GTG".split(",")),
    "*": list("TAA,TGA,TAG".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)
      
def complement (dna_seq):
    my_seq = Seq(dna_seq)
    my_com = my_seq.complement()
    return "".join(my_com)
        
def rev_complement (dna_seq):
    my_seq = Seq(dna_seq)
    my_rev_com = my_seq.reverse_complement()
    return "".join(my_rev_com)    
        
def plasmid_building (plasmid_seq, dna_seq):
    Backbone = Seq(plasmid_seq)
    Insert = Seq(dna_seq)
    Plasmid = Insert + Backbone
    return "".join(Plasmid)    

i=0
n = len(df_ntermini_id_description)
all_dna_seq = []
all_com_seq =[]
all_rev_com_seq = []
all_plasmids = []
all_ids = []

while i < n:
    dna_seq = aa2na(df_ntermini_id_description['ntermini'][i])
    com_dna_seq = complement(dna_seq)
    rev_com_dna_seq = rev_complement (dna_seq)
    insert_backbone_joined = plasmid_building (plasmid_str, dna_seq) 
    all_ids = np.append (all_ids, df_id [i])
    all_plasmids = np.append (all_plasmids, insert_backbone_joined)
    all_rev_com_seq = np.append(all_rev_com_seq, rev_com_dna_seq)
    all_com_seq = np.append(all_com_seq, com_dna_seq)
    all_dna_seq = np.append(all_dna_seq, dna_seq)
    i=i+1
print(all_dna_seq) 
print(all_com_seq) 
print (all_rev_com_seq)
print (all_plasmids)

df_all_ids = pd.DataFrame (all_ids, columns = ['ID'])
df_plasmids = pd.DataFrame(all_plasmids, columns = ['sequence'])
df_plasmids_id = pd.merge(df_all_ids, df_plasmids, left_index=True, right_index=True)
df_plasmids_id_desc = pd.merge(df_plasmids_id, df_description, left_index=True, right_index=True)

df_plasmids_id_desc.to_csv('salmonella_plasmids_ids_desc_full_dataset_shortend.csv')

plasmid_list = []

for i in range(len(df_plasmids_id_desc)):
    entry = SeqRecord(Seq(df_plasmids_id_desc['sequence'].iloc[i]),
    id=df_plasmids_id_desc['ID'].iloc[i])
    plasmid_list.append(entry)

def trim_backbone_primer_fw_right (dna_seq):
    trim_seq_right = dna_seq[:-11844]
    return (trim_seq_right) 
    
def trim_backbone_primer_fw_left (dna_seq):
    trim_seq_left = dna_seq[28:]
    return "".join(trim_seq_left) 

i=0
n = len(df_plasmids)
all_trim_seq = []
all_fw_primer = []

while i < n:
    trimmed_fw_primer_right = trim_backbone_primer_fw_right(all_plasmids[i])
    all_trim_seq = np.append (all_trim_seq, trimmed_fw_primer_right)
    trimmed_fw_primer_left = trim_backbone_primer_fw_left(all_trim_seq[i])
    all_fw_primer = np.append (all_fw_primer, trimmed_fw_primer_left)
    i=i+1

df_all_fw_primer= pd.DataFrame (all_fw_primer, columns = ['fw primer'])
print (df_all_fw_primer)

def trim_backbone_primer_rev (dna_seq):
    trim_seq_right2 = dna_seq[11915:]+dna_seq[:44]
    return (trim_seq_right2) 
    
i=0
n = len(df_plasmids)
all_rev_primer = []

while i < n:
    trimmed_plasmid = trim_backbone_primer_rev(all_plasmids [i])
    rev_com_primers = rev_complement (trimmed_plasmid)
    all_rev_primer = np.append (all_rev_primer, rev_com_primers)
    i=i+1
    
df_all_rev_primer= pd.DataFrame (all_rev_primer, columns = ['rev primer'])
print (df_all_rev_primer)

df_fw_and_rev_primers = pd.merge(df_all_fw_primer, df_all_rev_primer, left_index=True, right_index=True)
df_all_primers_ids = pd.merge(df_all_ids, df_fw_and_rev_primers, left_index=True, right_index=True)
df_fw_primers_ids=pd.merge(df_all_ids, df_all_fw_primer, left_index=True, right_index=True)
df_rev_primers_ids=pd.merge(df_all_ids, df_all_rev_primer, left_index=True, right_index=True)
df_fw_primers_ids_description = pd.merge (df_fw_primers_ids, df_description, left_index=True, right_index=True)
df_rev_primers_ids_description = pd.merge (df_rev_primers_ids, df_description, left_index=True, right_index=True)
df_all_primers_ids_description = pd.merge (df_all_rev_primer, df_fw_primers_ids_description, left_index=True, right_index=True)

df_all_primers_ids.to_csv('salmonella_primers_and_ids_full_dataset_shortend.csv')
df_all_primers_ids_description.to_csv('salmonella_primers_ids_descriptions_full_dataset_shortend.csv')