import numpy as np 
import random
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

##gering sekretierte unter 20% vom WT

df = pd.DataFrame ({'N-Terminus': ['MIPGTIPTSYLVPTADTEATGVVSL','MPITNASPENILRYLHAAGTGTKEA','MERSLDSLAGMAKSAFGAGTSAAMR','MATPWSGYLDDVSAKFDTGVDNLQT','MPITIGNGFLKSEILTNSPRNTKEA','MVTSVRTQPPVIMPGMQTEIKTQAT','MPVTLSFGNRHNYEINHSRLARLMS','MRVSGSASSQDIISRINSKNINNND','MSVVPVSTQSYVKSSAEPSQEQINF','MSSGNILWGSQNPIVFKNSFGVSNA','MPFHIGSGCLPAIISNRRIYRIAWS','MPLSVGQGYFTSSISSEKFNAIKES','MIPPLNRYVPALSKNELVKTVTNRD','MARFNAAFTRIKIMFSRIRGLISCQ','MFSRVRGFLSCQNYSHTATPAITLP','MPFHIGSGCLPATISNRRIYRIAWS'],
                    'ID': ['OrgC', 'PipB', 'PipB2', 'PrgI', 'SifA', 'SipA', 'SopD2', 'SpvD', 'SsaL', 'SseB', 'SseI', 'SseJ', 'SseK1', 'SseK2', 'SseK3', 'SspH2'],
                    'description': ['pJA027', 'pJA029', 'pJA030', 'pJA031', 'pJA037', 'pJA039', 'pJA047', 'pJA053', 'pJA057', 'pJA058', 'pJA063', 'pJA064', 'pJA065', 'pJA066', 'pJA067', 'pJA070']})


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
n = len(df)
all_as_properties = []
all_ids = []
all_desc = []

while i < n:
    as_properties = aa2na(df['N-Terminus'][i])
    all_ids = np.append (all_ids, df['ID'][i])
    all_desc = np.append (all_desc, df ['description'] [i])
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
print(df_properties_per_residue [df_properties_per_residue == " M"].count())
print(df_properties_per_residue [df_properties_per_residue == "M"].count())
print(df_properties_per_residue [df_properties_per_residue == " F"].count())
print(df_properties_per_residue [df_properties_per_residue == " P"].count())
print(df_properties_per_residue [df_properties_per_residue == " S"].count())
print(df_properties_per_residue [df_properties_per_residue == " T"].count())
print(df_properties_per_residue [df_properties_per_residue == " W"].count())
print(df_properties_per_residue [df_properties_per_residue == " Y"].count())
print(df_properties_per_residue [df_properties_per_residue == " V"].count())

