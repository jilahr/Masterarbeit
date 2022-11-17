import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import stats
import scipy.stats as sp 

#######
##Assay 1: 13.09.22
#######

file1 = '20220913_nl_plate_01_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pja067
file2 = '20220913_rfl_plate_01_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen + pja067
file3 = '20220913_nl_plate_02_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pja074
file4 = '20220913_rfl_plate_02_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pja074
file5 = '20220913_nl_plate_03_001.csv' # nanoluc luciferase plate 3 pJA076, 035, 036
file6 = '20220913_rfl_plate_03_001.csv' # red firefly luciferase plate 3 pJA076, 035, 036
file7 = '20220913_nl_plate_04_001.csv' # nanoluc luciferase plate 4 pJA037, 038, 039
file8 = '20220913_rfl_plate_04_001.csv' # red firefly luciferase plate 4 pJA037, 038, 039

# nanoluc luciferase plate 1 positiv Kontrollen + pJA067

df_p1 = pd.read_csv(file1, skiprows=1)
df_p1 = df_p1.iloc[4:12,:12]
df_p1.columns = ['L1','L2','L3']*4
df_p1 = df_p1.reset_index(drop=True)

s1 = df_p1.iloc[:8,:3]
s2 = df_p1.iloc[:8,4:7]
s3 = df_p1.iloc[:8,8:11]

df_pc1_nl = pd.concat([s1])
df_pc1_nl['number'] = np.arange(0,8,1)
df_pc1_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.1'
df_pc1_nl['Konstrukt'] = 'Positivkontrolle.1'
df_pc1_nl['Luciferase'] = 'NL'
df_pc1_nl['N-Terminus'] = 'PK'
df_pc1_nl['Normalisierung']='PC1'

df_pc2_nl = pd.concat([s2])
df_pc2_nl['number'] = np.arange(0,8,1)
df_pc2_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.2'
df_pc2_nl['Konstrukt'] = 'Positivkontrolle.2'
df_pc2_nl['Luciferase'] = 'NL'
df_pc2_nl['N-Terminus'] = 'PC2 SptP'
df_pc2_nl['Normalisierung']='PC2'

df_67_nl = pd.concat([s3])
df_67_nl['number'] = np.arange(0,8,1)
df_67_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA067'
df_67_nl['Konstrukt'] = 'P0DUJ7' #protein ID
df_67_nl['Luciferase'] = 'NL'
df_67_nl['N-Terminus'] = 'SseK3'
df_67_nl['Normalisierung']='PC1'

# red firefly luciferase plate 1 positiv Kontrollen + pJA067

df_p2 = pd.read_csv(file2, skiprows=1)
df_p2 = df_p2.iloc[4:12,:12]
df_p2.columns = ['L1','L2','L3']*4
df_p2 = df_p2.reset_index(drop=True)

s4 = df_p2.iloc[:8,:3]
s5 = df_p2.iloc[:8,4:7]
s6 = df_p2.iloc[:8,8:11]

df_pc1_rfl = pd.concat([s4])
df_pc1_rfl['number'] = np.arange(0,8,1)
df_pc1_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.1'
df_pc1_rfl['Konstrukt'] = 'Positivkontrolle.1'
df_pc1_rfl['Luciferase'] = 'RFL'
df_pc1_rfl['N-Terminus'] = 'PK'
df_pc1_rfl['Normalisierung']='PC1'

df_pc2_rfl = pd.concat([s5])
df_pc2_rfl['number'] = np.arange(0,8,1)
df_pc2_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.2'
df_pc2_rfl['Konstrukt'] = 'Positivkontrolle.2'
df_pc2_rfl['Luciferase'] = 'RFL'
df_pc2_rfl['N-Terminus'] = 'PC2 SptP'
df_pc2_rfl['Normalisierung']='PC2'

df_67_rfl = pd.concat([s6])
df_67_rfl['number'] = np.arange(0,8,1)
df_67_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA067'
df_67_rfl['Konstrukt'] = 'P0DUJ7' # protein ID
df_67_rfl['Luciferase'] = 'RFL'
df_67_rfl['N-Terminus'] = 'SseK3'
df_67_rfl['Normalisierung']='PC1'

# nanoluc luciferase plate 2 negativ Kontrollen + pJA074

df_p3 = pd.read_csv(file3, skiprows=1)
df_p3 = df_p3.iloc[4:12,:12]
df_p3.columns = ['L1','L2','L3']*4
df_p3 = df_p3.reset_index(drop=True)

s7 = df_p3.iloc[:8,:3]
s8 = df_p3.iloc[:8,4:7]
s9 = df_p3.iloc[:8,8:11]

df_pc_leer_nl = pd.concat([s7])
df_pc_leer_nl['number'] = np.arange(0,8,1)
df_pc_leer_nl['Linie'] = 'SB905 dSipA dSptP pHilA'
df_pc_leer_nl['Konstrukt'] = 'Leerkontrolle1'
df_pc_leer_nl['Luciferase'] = 'NL'
df_pc_leer_nl['N-Terminus'] = 'leer1'
df_pc_leer_nl['Normalisierung']='PC1'

df_nc_nl = pd.concat([s8])
df_nc_nl['number'] = np.arange(0,8,1)
df_nc_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049'
df_nc_nl['Konstrukt'] = 'Negativkontrolle.1'
df_nc_nl['Luciferase'] = 'NL'
df_nc_nl['N-Terminus'] = 'NC1 SptP'
df_nc_nl['Normalisierung']='PC1'

df_74_nl = pd.concat([s9])
df_74_nl['number'] = np.arange(0,8,1)
df_74_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA074'
df_74_nl['Konstrukt'] = 'Q8ZNP2' # protein ID
df_74_nl['Luciferase'] = 'NL'
df_74_nl['N-Terminus'] = 'SteD'
df_74_nl['Normalisierung']='PC1'

# red firefly luciferase plate 2 negativ Kontrollen + pJA074

df_p4 = pd.read_csv(file4, skiprows=1)
df_p4 = df_p4.iloc[4:12,:12]
df_p4.columns = ['L1','L2','L3']*4
df_p4 = df_p4.reset_index(drop=True)

s10 = df_p4.iloc[:8,:3]
s11 = df_p4.iloc[:8,4:7]
s12 = df_p4.iloc[:8,8:11]
  
df_pc_leer_rfl = pd.concat([s10])
df_pc_leer_rfl['number'] = np.arange(0,8,1)
df_pc_leer_rfl['Linie'] = 'SB905 dSipA dSptP pHilA'
df_pc_leer_rfl['Konstrukt'] = 'Leerkontrolle1'
df_pc_leer_rfl['Luciferase'] = 'RFL'
df_pc_leer_rfl['N-Terminus'] = 'leer1'
df_pc_leer_rfl['Normalisierung']='PC1'

df_nc_rfl = pd.concat([s11])
df_nc_rfl['number'] = np.arange(0,8,1)
df_nc_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049'
df_nc_rfl['Konstrukt'] = 'Negativkontrolle.1'
df_nc_rfl['Luciferase'] = 'RFL'
df_nc_rfl['N-Terminus'] = 'NC1 SptP'
df_nc_rfl['Normalisierung']='PC1'

df_74_rfl = pd.concat([s12])
df_74_rfl['number'] = np.arange(0,8,1)
df_74_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA074'
df_74_rfl['Konstrukt'] = 'Q8ZNP2' # protein ID
df_74_rfl['Luciferase'] = 'RFL'
df_74_nl['N-Terminus'] = 'SteD'
df_74_rfl['Normalisierung']='PC1'

# nanoluc luciferase plate 3 pJA076, 035, 036

df_p5 = pd.read_csv(file5, skiprows=1)
df_p5 = df_p5.iloc[4:12,:12]
df_p5.columns = ['L1','L2','L3']*4
df_p5 = df_p5.reset_index(drop=True)

s13 = df_p5.iloc[:8,:3]
s14 = df_p5.iloc[:8,4:7]
s15 = df_p5.iloc[:8,8:11]

df_76_nl = pd.concat([s13])
df_76_nl['number'] = np.arange(0,8,1)
df_76_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA076'
df_76_nl['Konstrukt'] = 'Q8Z7T2' # protein ID
df_76_nl['Luciferase'] = 'NL'
df_76_nl['N-Terminus'] = 'StoD'
df_76_nl['Normalisierung']='PC1'

df_35_nl = pd.concat([s14])
df_35_nl['number'] = np.arange(0,8,1)
df_35_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA035'
df_35_nl['Konstrukt'] = 'A0A0K0HC32' # protein ID
df_35_nl['Luciferase'] = 'NL'
df_35_nl['N-Terminus'] = 'SboH'
df_35_nl['Normalisierung']='PC2'

df_36_nl = pd.concat([s15])
df_36_nl['number'] = np.arange(0,8,1)
df_36_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA036'
df_36_nl['Konstrukt'] = 'A0A0K0H9B7' # protein ID
df_36_nl['Luciferase'] = 'NL'
df_36_nl['N-Terminus'] = 'SboI'
df_36_nl['Normalisierung']='PC2'

# red firefly luciferase plate 3 pJA076, 035, 036

df_p6 = pd.read_csv(file6, skiprows=1)
df_p6 = df_p6.iloc[4:12,:12]
df_p6.columns = ['L1','L2','L3']*4
df_p6 = df_p6.reset_index(drop=True)

s16 = df_p6.iloc[:8,:3]
s17 = df_p6.iloc[:8,4:7]
s18 = df_p6.iloc[:8,8:11]

df_76_rfl = pd.concat([s16])
df_76_rfl['number'] = np.arange(0,8,1)
df_76_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA076'
df_76_rfl['Konstrukt'] = 'Q8Z7T2' # protein ID
df_76_rfl['Luciferase'] = 'RFL'
df_76_rfl['N-Terminus'] = 'StoD'
df_76_rfl['Normalisierung']='PC1'

df_35_rfl = pd.concat([s17])
df_35_rfl['number'] = np.arange(0,8,1)
df_35_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA035'
df_35_rfl['Konstrukt'] = 'A0A0K0HC32' # protein ID
df_35_rfl['Luciferase'] = 'RFL'
df_35_rfl['N-Terminus'] = 'SboH'
df_35_rfl['Normalisierung']='PC2'

df_36_rfl = pd.concat([s18])
df_36_rfl['number'] = np.arange(0,8,1)
df_36_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA036'
df_36_rfl['Konstrukt'] = 'A0A0K0H9B7' # protein ID
df_36_rfl['Luciferase'] = 'RFL'
df_36_rfl['N-Terminus'] = 'SboI'
df_36_rfl['Normalisierung']='PC2'

# nanoluc luciferase plate 4 pJA037, 038, 039

df_p7 = pd.read_csv(file7, skiprows=1)
df_p7 = df_p7.iloc[4:12,:12]
df_p7.columns = ['L1','L2','L3']*4
df_p7 = df_p7.reset_index(drop=True)

s19 = df_p7.iloc[:8,:3]
s20 = df_p7.iloc[:8,4:7]
s21 = df_p7.iloc[:8,8:11]

df_37_nl = pd.concat([s19])
df_37_nl['number'] = np.arange(0,8,1)
df_37_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA037'
df_37_nl['Konstrukt'] = 'Q56061' # protein ID
df_37_nl['Luciferase'] = 'NL'
df_37_nl['N-Terminus'] = 'SifA'
df_37_nl['Normalisierung']='PC2'

df_38_nl = pd.concat([s20])
df_38_nl['number'] = np.arange(0,8,1)
df_38_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA038'
df_38_nl['Konstrukt'] = 'Q9KIB9' # protein ID
df_38_nl['Luciferase'] = 'NL'
df_38_nl['N-Terminus'] = 'SifB'
df_38_nl['Normalisierung']='PC2'

df_39_nl = pd.concat([s21])
df_39_nl['number'] = np.arange(0,8,1)
df_39_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA039'
df_39_nl['Konstrukt'] = 'P0CL52' # protein ID
df_39_nl['Luciferase'] = 'NL'
df_39_nl['N-Terminus'] = 'SipA'
df_39_nl['Normalisierung']='PC2'

# red firefly luciferase plate 4 pJA037, 038, 039

df_p8 = pd.read_csv(file8, skiprows=1)
df_p8 = df_p8.iloc[4:12,:12]
df_p8.columns = ['L1','L2','L3']*4
df_p8 = df_p8.reset_index(drop=True)

s22 = df_p8.iloc[:8,:3]
s23 = df_p8.iloc[:8,4:7]
s24 = df_p8.iloc[:8,8:11]

df_37_rfl = pd.concat([s22])
df_37_rfl['number'] = np.arange(0,8,1)
df_37_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA037'
df_37_rfl['Konstrukt'] = 'Q56061' # protein ID
df_37_rfl['Luciferase'] = 'RFL'
df_37_rfl['N-Terminus'] = 'SifA'
df_37_rfl['Normalisierung']='PC2'

df_38_rfl = pd.concat([s23])
df_38_rfl['number'] = np.arange(0,8,1)
df_38_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA038'
df_38_rfl['Konstrukt'] = 'Q9KIB9' # protein ID
df_38_rfl['Luciferase'] = 'RFL'
df_38_rfl['N-Terminus'] = 'SifB'
df_38_rfl['Normalisierung']='PC2'

df_39_rfl = pd.concat([s24])
df_39_rfl['number'] = np.arange(0,8,1)
df_39_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA039'
df_39_rfl['Konstrukt'] = 'P0CL52' # protein ID
df_39_rfl['Luciferase'] = 'RFL'
df_39_rfl['N-Terminus'] = 'SipA'
df_39_rfl['Normalisierung']='PC2'


###########
#Assay 2: 06.09.22
###########

file9 = '20220906_plate_01_NL_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pja031
file10 = '20220906_plate_01_RFL_10min_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen + pja031
file11 = '20220906_plate_02_NL_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pja032
file12 = '20220906_plate_02_RFL_10min_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pja032
file13 = '20220906_plate_03_NL_001.csv' # nanoluc luciferase plate 3 pJA033, 034, 068
file14 = '20220906_plate_03_RFL_10min_001.csv' # red firefly luciferase plate 3 pJA033, 034, 068
file15 = '20220906_plate_04_NL_001.csv' # nanoluc luciferase plate 4 pJA070, 072, 073
file16 = '20220906_plate_04_RFL_10min_001.csv' # red firefly luciferase plate 4 pJA070,072, 073

# nanoluc luciferase plate 1 positiv Kontrollen + pJA031

df_p9 = pd.read_csv(file9, skiprows=1)
df_p9 = df_p9.iloc[4:12,:12]
df_p9.columns = ['L1','L2','L3']*4
df_p9 = df_p9.reset_index(drop=True)

s25 = df_p9.iloc[:8,:3]
s26 = df_p9.iloc[:8,4:7]
s27 = df_p9.iloc[:8,8:11]

df_pc3_nl = pd.concat([s25])
df_pc3_nl['number'] = np.arange(0,8,1)
df_pc3_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.3'
df_pc3_nl['Konstrukt'] = 'Positivkontrolle.3'
df_pc3_nl['Luciferase'] = 'NL'
df_pc3_nl['N-Terminus'] = 'PC3 SptP'
df_pc3_nl['Normalisierung']='PC3'

df_pc4_nl = pd.concat([s26])
df_pc4_nl['number'] = np.arange(0,8,1)
df_pc4_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.4'
df_pc4_nl['Konstrukt'] = 'Positivkontrolle.4'
df_pc4_nl['Luciferase'] = 'NL'
df_pc4_nl['N-Terminus'] = 'PC4 SptP'
df_pc4_nl['Normalisierung']='PC4'

df_31_nl = pd.concat([s27])
df_31_nl['number'] = np.arange(0,8,1)
df_31_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA031'
df_31_nl['Konstrukt'] = 'P41784' # protein ID
df_31_nl['Luciferase'] = 'NL'
df_31_nl['N-Terminus'] = 'PrgI'
df_31_nl['Normalisierung']='PC3'

# red firefly luciferase plate 1 positiv Kontrollen + pJA031

df_p10 = pd.read_csv(file10, skiprows=1)
df_p10 = df_p10.iloc[4:12,:12]
df_p10.columns = ['L1','L2','L3']*4
df_p10 = df_p10.reset_index(drop=True)

s28 = df_p10.iloc[:8,:3]
s29 = df_p10.iloc[:8,4:7]
s30 = df_p10.iloc[:8,8:11]

df_pc3_rfl = pd.concat([s28])
df_pc3_rfl['number'] = np.arange(0,8,1)
df_pc3_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.3'
df_pc3_rfl['Konstrukt'] = 'Positivkontrolle.3'
df_pc3_rfl['Luciferase'] = 'RFL'
df_pc3_rfl['N-Terminus'] = 'PC3 SptP'
df_pc3_rfl['Normalisierung']='PC3'

df_pc4_rfl = pd.concat([s29])
df_pc4_rfl['number'] = np.arange(0,8,1)
df_pc4_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.4'
df_pc4_rfl['Konstrukt'] = 'Positivkontrolle.4'
df_pc4_rfl['Luciferase'] = 'RFL'
df_pc4_rfl['N-Terminus'] = 'PC4 SptP'
df_pc4_rfl['Normalisierung']='PC4'

df_31_rfl = pd.concat([s30])
df_31_rfl['number'] = np.arange(0,8,1)
df_31_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA031'
df_31_rfl['Konstrukt'] = 'P41784' # protein ID
df_31_rfl['Luciferase'] = 'RFL'
df_31_rfl['N-Terminus'] = 'PrgI'
df_31_rfl['Normalisierung']='PC3'

# nanoluc luciferase plate 2 negativ Kontrollen + pja032

df_p11 = pd.read_csv(file11, skiprows=1)
df_p11 = df_p11.iloc[4:12,:12]
df_p11.columns = ['L1','L2','L3']*4
df_p11 = df_p11.reset_index(drop=True)

s31 = df_p11.iloc[:8,:3]
s32 = df_p11.iloc[:8,4:7]
s33 = df_p11.iloc[:8,8:11]

df_pc_leer2_nl = pd.concat([s32])
df_pc_leer2_nl['number'] = np.arange(0,8,1)
df_pc_leer2_nl['Linie'] = 'SB905 dSipA dSptP pHilA.2'
df_pc_leer2_nl['Konstrukt'] = 'Leerkontrolle2'
df_pc_leer2_nl['Luciferase'] = 'NL'
df_pc_leer2_nl['N-Terminus'] = 'LK'
df_pc_leer2_nl['Normalisierung']='PC3'

df_nc2_nl = pd.concat([s31])
df_nc2_nl['number'] = np.arange(0,8,1)
df_nc2_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.2'
df_nc2_nl['Konstrukt'] = 'Negativkontrolle.2'
df_nc2_nl['Luciferase'] = 'NL'
df_nc2_nl['N-Terminus'] = 'NK'
df_nc2_nl['Normalisierung']='PC3'

df_32_nl = pd.concat([s33])
df_32_nl['number'] = np.arange(0,8,1)
df_32_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA032'
df_32_nl['Konstrukt'] = 'P41785' # protein ID
df_32_nl['Luciferase'] = 'NL'
df_32_nl['N-Terminus'] = 'PrgJ'
df_32_nl['Normalisierung']='PC3'

# red firefly luciferase plate 2 negativ Kontrollen + pJA032

df_p12 = pd.read_csv(file12, skiprows=1)
df_p12 = df_p12.iloc[4:12,:12]
df_p12.columns = ['L1','L2','L3']*4
df_p12 = df_p12.reset_index(drop=True)

s34 = df_p12.iloc[:8,:3]
s35 = df_p12.iloc[:8,4:7]
s36 = df_p12.iloc[:8,8:11]
  
df_pc_leer2_rfl = pd.concat([s35])
df_pc_leer2_rfl['number'] = np.arange(0,8,1)
df_pc_leer2_rfl['Linie'] = 'SB905 dSipA dSptP pHilA.2'
df_pc_leer2_rfl['Konstrukt'] = 'Leerkontrolle2'
df_pc_leer2_rfl['Luciferase'] = 'RFL'
df_pc_leer2_rfl['N-Terminus'] = 'LK'
df_pc_leer2_rfl['Normalisierung']='PC3'

df_nc2_rfl = pd.concat([s34])
df_nc2_rfl['number'] = np.arange(0,8,1)
df_nc2_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.2'
df_nc2_rfl['Konstrukt'] = 'Negativkontrolle.2'
df_nc2_rfl['Luciferase'] = 'RFL'
df_nc2_rfl['N-Terminus'] = 'NK'
df_nc2_rfl['Normalisierung']='PC3'

df_32_rfl = pd.concat([s36])
df_32_rfl['number'] = np.arange(0,8,1)
df_32_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA032'
df_32_rfl['Konstrukt'] = 'P41785' # protein ID
df_32_rfl['Luciferase'] = 'RFL'
df_32_rfl['N-Terminus'] = 'PrgJ'
df_32_rfl['Normalisierung']='PC3'

# nanoluc luciferase plate 3 pJA033, 034, 068

df_p13 = pd.read_csv(file13, skiprows=1)
df_p13 = df_p13.iloc[4:12,:12]
df_p13.columns = ['L1','L2','L3']*4
df_p13 = df_p13.reset_index(drop=True)

s37 = df_p13.iloc[:8,:3]
s38 = df_p13.iloc[:8,4:7]
s39 = df_p13.iloc[:8,8:11]

df_33_nl = pd.concat([s37])
df_33_nl['number'] = np.arange(0,8,1)
df_33_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA033'
df_33_nl['Konstrukt'] = 'A0A0K0H8V0' # protein ID
df_33_nl['Luciferase'] = 'NL'
df_33_nl['N-Terminus'] = 'SboA'
df_33_nl['Normalisierung']='PC3'

df_34_nl = pd.concat([s38])
df_34_nl['number'] = np.arange(0,8,1)
df_34_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA034'
df_34_nl['Konstrukt'] = 'A0A0K0HD42' # protein ID
df_34_nl['Luciferase'] = 'NL'
df_34_nl['N-Terminus'] = 'SboC'
df_34_nl['Normalisierung']='PC4'

df_68_nl = pd.concat([s39])
df_68_nl['number'] = np.arange(0,8,1)
df_68_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA068'
df_68_nl['Konstrukt'] = 'Q8ZNG2' # protein ID
df_68_nl['Luciferase'] = 'NL'
df_68_nl['N-Terminus'] = 'SseL'
df_68_nl['Normalisierung']='PC4'

# red firefly luciferase plate 3 pJA033, 034, 068

df_p14 = pd.read_csv(file14, skiprows=1)
df_p14 = df_p14.iloc[4:12,:12]
df_p14.columns = ['L1','L2','L3']*4
df_p14 = df_p14.reset_index(drop=True)

s40 = df_p14.iloc[:8,:3]
s41 = df_p14.iloc[:8,4:7]
s42 = df_p14.iloc[:8,8:11]

df_33_rfl = pd.concat([s40])
df_33_rfl['number'] = np.arange(0,8,1)
df_33_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA033'
df_33_rfl['Konstrukt'] = 'A0A0K0H8V0' # protein ID
df_33_rfl['Luciferase'] = 'RFL'
df_33_rfl['N-Terminus'] = 'SboA'
df_33_rfl['Normalisierung']='PC3'

df_34_rfl = pd.concat([s41])
df_34_rfl['number'] = np.arange(0,8,1)
df_34_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA034'
df_34_rfl['Konstrukt'] = 'A0A0K0HD42' # protein ID
df_34_rfl['Luciferase'] = 'RFL'
df_34_rfl['N-Terminus'] = 'SboC'
df_34_rfl['Normalisierung']='PC4'

df_68_rfl = pd.concat([s42])
df_68_rfl['number'] = np.arange(0,8,1)
df_68_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA068'
df_68_rfl['Konstrukt'] = 'Q8ZNG2' # protein ID
df_68_rfl['Luciferase'] = 'RFL'
df_68_rfl['N-Terminus'] = 'SseL'
df_68_rfl['Normalisierung']='PC4'

# nanoluc luciferase plate 4 pJA070,072, 073

df_p15 = pd.read_csv(file15, skiprows=1)
df_p15 = df_p15.iloc[4:12,:12]
df_p15.columns = ['L1','L2','L3']*4
df_p15 = df_p15.reset_index(drop=True)

s43 = df_p15.iloc[:8,:3]
s44 = df_p15.iloc[:8,4:7]
s45 = df_p15.iloc[:8,8:11]

df_70_nl = pd.concat([s43])
df_70_nl['number'] = np.arange(0,8,1)
df_70_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA070'
df_70_nl['Konstrukt'] = 'P0CE12' # protein ID
df_70_nl['Luciferase'] = 'NL'
df_70_nl['N-Terminus'] = 'SspH2'
df_70_nl['Normalisierung']='PC4'

df_72_nl = pd.concat([s44])
df_72_nl['number'] = np.arange(0,8,1)
df_72_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA072'
df_72_nl['Konstrukt'] = 'Q8ZPA6' # protein ID
df_72_nl['Luciferase'] = 'NL'
df_72_nl['N-Terminus'] = 'SteB'
df_72_nl['Normalisierung']='PC4'

df_73_nl = pd.concat([s45])
df_73_nl['number'] = np.arange(0,8,1)
df_73_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA073'
df_73_nl['Konstrukt'] = 'Q8ZP57' # protein ID
df_73_nl['Luciferase'] = 'NL'
df_73_nl['N-Terminus'] = 'SteC'
df_73_nl['Normalisierung']='PC4'

# red firefly luciferase plate 4 pJA070,072, 073

df_p16 = pd.read_csv(file16, skiprows=1)
df_p16 = df_p16.iloc[4:12,:12]
df_p16.columns = ['L1','L2','L3']*4
df_p16 = df_p16.reset_index(drop=True)

s46 = df_p16.iloc[:8,:3]
s47 = df_p16.iloc[:8,4:7]
s48 = df_p16.iloc[:8,8:11]

df_70_rfl = pd.concat([s46])
df_70_rfl['number'] = np.arange(0,8,1)
df_70_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA070'
df_70_rfl['Konstrukt'] = 'P0CE12' # protein ID
df_70_rfl['Luciferase'] = 'RFL'
df_70_rfl['N-Terminus'] = 'SspH2'
df_70_rfl['Normalisierung']='PC4'

df_72_rfl = pd.concat([s47])
df_72_rfl['number'] = np.arange(0,8,1)
df_72_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA072'
df_72_rfl['Konstrukt'] = 'Q8ZPA6' # protein ID
df_72_rfl['Luciferase'] = 'RFL'
df_72_rfl['N-Terminus'] = 'SteB'
df_72_rfl['Normalisierung']='PC4'

df_73_rfl = pd.concat([s48])
df_73_rfl['number'] = np.arange(0,8,1)
df_73_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA073'
df_73_rfl['Konstrukt'] = 'Q8ZP57' # protein ID
df_73_rfl['Luciferase'] = 'RFL'
df_73_rfl['N-Terminus'] = 'SteC'
df_73_rfl['Normalisierung']='PC4'


###########
#Assay 3: 01.09.22
###########

file17 = '20220901 nl 01_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pja029
file18 = '20220901 rfl 01_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen + pja029
file19 = '20220901 nl 02_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pja040
file20 = '20220901 rfl 02_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pja040
file21 = '20220901 nl 03_001.csv' # nanoluc luciferase plate 3 pJA043, 044, 048
file22 = '20220901 rfl 03_001.csv' # red firefly luciferase plate 3 pJA043, 044, 048
file23 = '20220901 nl 04_001.csv' # nanoluc luciferase plate 4 pJA058, 065, 066
file24 = '20220901 rfl 04_001.csv' # red firefly luciferase plate 4 pJA058,065, 066

# nanoluc luciferase plate 1 positiv Kontrollen + pJA029

df_p17 = pd.read_csv(file17, skiprows=1)
df_p17 = df_p17.iloc[4:12,:12]
df_p17.columns = ['L1','L2','L3']*4
df_p17 = df_p17.reset_index(drop=True)

s49 = df_p17.iloc[:8,:3]
s50 = df_p17.iloc[:8,4:7]
s51 = df_p17.iloc[:8,8:11]

df_pc5_nl = pd.concat([s49])
df_pc5_nl['number'] = np.arange(0,8,1)
df_pc5_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.5'
df_pc5_nl['Konstrukt'] = 'Positivkontrolle.5'
df_pc5_nl['Luciferase'] = 'NL'
df_pc5_nl['N-Terminus'] = 'PC5 SptP'
df_pc5_nl['Normalisierung']='PC5'

df_pc6_nl = pd.concat([s50])
df_pc6_nl['number'] = np.arange(0,8,1)
df_pc6_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.6'
df_pc6_nl['Konstrukt'] = 'Positivkontrolle.6'
df_pc6_nl['Luciferase'] = 'NL'
df_pc6_nl['N-Terminus'] = 'PC6 SptP'
df_pc6_nl['Normalisierung']='PC6'

df_29_nl = pd.concat([s51])
df_29_nl['number'] = np.arange(0,8,1)
df_29_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA029'
df_29_nl['Konstrukt'] = 'Q8ZQ59' # protein ID
df_29_nl['Luciferase'] = 'NL'
df_29_nl['N-Terminus'] = 'PipB'
df_29_nl['Normalisierung']='PC5'

# red firefly luciferase plate 1 positiv Kontrollen + pJA029

df_p18 = pd.read_csv(file18, skiprows=1)
df_p18 = df_p18.iloc[4:12,:12]
df_p18.columns = ['L1','L2','L3']*4
df_p18 = df_p18.reset_index(drop=True)

s52 = df_p18.iloc[:8,:3]
s53 = df_p18.iloc[:8,4:7]
s54 = df_p18.iloc[:8,8:11]
  
df_pc5_rfl = pd.concat([s52])
df_pc5_rfl['number'] = np.arange(0,8,1)
df_pc5_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.5'
df_pc5_rfl['Konstrukt'] = 'Positivkontrolle.5'
df_pc5_rfl['Luciferase'] = 'RFL'
df_pc5_rfl['N-Terminus'] = 'PC5 SptP'
df_pc5_rfl['Normalisierung']='PC5'

df_pc6_rfl = pd.concat([s53])
df_pc6_rfl['number'] = np.arange(0,8,1)
df_pc6_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.6'
df_pc6_rfl['Konstrukt'] = 'Positivkontrolle.6'
df_pc6_rfl['Luciferase'] = 'RFL'
df_pc6_rfl['N-Terminus'] = 'PC6 SptP'
df_pc6_rfl['Normalisierung']='PC6'

df_29_rfl = pd.concat([s54])
df_29_rfl['number'] = np.arange(0,8,1)
df_29_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA029'
df_29_rfl['Konstrukt'] = 'Q8ZQ59' # protein ID
df_29_rfl['Luciferase'] = 'RFL'
df_29_rfl['N-Terminus'] = 'PipB'
df_29_rfl['Normalisierung']='PC5'

# nanoluc luciferase plate 2 negativ Kontrollen + pJA040

df_p19 = pd.read_csv(file19, skiprows=1)
df_p19 = df_p19.iloc[4:12,:12]
df_p19.columns = ['L1','L2','L3']*4
df_p19 = df_p19.reset_index(drop=True)

s55 = df_p19.iloc[:8,:3]
s56 = df_p19.iloc[:8,4:7]
s57 = df_p19.iloc[:8,8:11]

df_pc_leer3_nl = pd.concat([s55])
df_pc_leer3_nl['number'] = np.arange(0,8,1)
df_pc_leer3_nl['Linie'] = 'SB905 dSipA dSptP pHilA.3'
df_pc_leer3_nl['Konstrukt'] = 'Leerkontrolle3'
df_pc_leer3_nl['Luciferase'] = 'NL'
df_pc_leer3_nl['N-Terminus'] = 'leer3'
df_pc_leer3_nl['Normalisierung']='PC5'

df_nc3_nl = pd.concat([s56])
df_nc3_nl['number'] = np.arange(0,8,1)
df_nc3_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.3'
df_nc3_nl['Konstrukt'] = 'Negativkontrolle.3'
df_nc3_nl['Luciferase'] = 'NL'
df_nc3_nl['N-Terminus'] = 'NC3 SptP'
df_nc3_nl['Normalisierung']='PC5'

df_40_nl = pd.concat([s57])
df_40_nl['number'] = np.arange(0,8,1)
df_40_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA040'
df_40_nl['Konstrukt'] = 'Q56019' # protein ID
df_40_nl['Luciferase'] = 'NL'
df_40_nl['N-Terminus'] = 'SipB'
df_40_nl['Normalisierung']='PC5'

# red firefly luciferase plate 2 negativ Kontrollen + pJA040

df_p20 = pd.read_csv(file20, skiprows=1)
df_p20 = df_p20.iloc[4:12,:12]
df_p20.columns = ['L1','L2','L3']*4
df_p20 = df_p20.reset_index(drop=True)

s58 = df_p20.iloc[:8,:3]
s59 = df_p20.iloc[:8,4:7]
s60 = df_p20.iloc[:8,8:11]
  
df_pc_leer3_rfl = pd.concat([s58])
df_pc_leer3_rfl['number'] = np.arange(0,8,1)
df_pc_leer3_rfl['Linie'] = 'SB905 dSipA dSptP pHilA.3'
df_pc_leer3_rfl['Konstrukt'] = 'Leerkontrolle3'
df_pc_leer3_rfl['Luciferase'] = 'RFL'
df_pc_leer3_rfl['N-Terminus'] = 'leer3'
df_pc_leer3_rfl['Normalisierung']='PC5'

df_nc3_rfl = pd.concat([s59])
df_nc3_rfl['number'] = np.arange(0,8,1)
df_nc3_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.3'
df_nc3_rfl['Konstrukt'] = 'Negativkontrolle.3'
df_nc3_rfl['Luciferase'] = 'RFL'
df_nc3_rfl['N-Terminus'] = 'NC3 SptP'
df_nc3_rfl['Normalisierung']='PC5'

df_40_rfl = pd.concat([s60])
df_40_rfl['time'] = 4
df_40_rfl['number'] = np.arange(0,8,1)
df_40_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA040'
df_40_rfl['Konstrukt'] = 'Q56019' # protein ID
df_40_rfl['Luciferase'] = 'RFL'
df_40_rfl['N-Terminus'] = 'SipB'
df_40_rfl['Normalisierung']='PC5'

# nanoluc luciferase plate 3 pJA043, 044, 048

df_p21 = pd.read_csv(file21, skiprows=1)
df_p21 = df_p21.iloc[4:12,:12]
df_p21.columns = ['L1','L2','L3']*4
df_p21 = df_p21.reset_index(drop=True)

s61 = df_p21.iloc[:8,:3]
s62 = df_p21.iloc[:8,4:7]
s63 = df_p21.iloc[:8,8:11]

df_43_nl = pd.concat([s61])
df_43_nl['number'] = np.arange(0,8,1)
df_43_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA043'
df_43_nl['Konstrukt'] = 'Q8ZQQ2' # protein ID
df_43_nl['Luciferase'] = 'NL'
df_43_nl['N-Terminus'] = 'SlrP'
df_43_nl['Normalisierung']='PC5'

df_44_nl = pd.concat([s62])
df_44_nl['number'] = np.arange(0,8,1)
df_44_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA044'
df_44_nl['Konstrukt'] = 'Q8ZNR3' # protein ID
df_44_nl['Luciferase'] = 'NL'
df_44_nl['N-Terminus'] = 'SopA'
df_44_nl['Normalisierung']='PC6'

df_48_nl = pd.concat([s63])
df_48_nl['number'] = np.arange(0,8,1)
df_48_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA048'
df_48_nl['Konstrukt'] = 'O52623' # protein ID
df_48_nl['Luciferase'] = 'NL'
df_48_nl['N-Terminus'] = 'SopE'
df_48_nl['Normalisierung']='PC6'

# red firefly luciferase plate 3 pJA043, 044, 048

df_p22 = pd.read_csv(file22, skiprows=1)
df_p22 = df_p22.iloc[4:12,:12]
df_p22.columns = ['L1','L2','L3']*4
df_p22 = df_p22.reset_index(drop=True)

s64 = df_p22.iloc[:8,:3]
s65 = df_p22.iloc[:8,4:7]
s66 = df_p22.iloc[:8,8:11]

df_43_rfl = pd.concat([s64])
df_43_rfl['number'] = np.arange(0,8,1)
df_43_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA043'
df_43_rfl['Konstrukt'] = 'Map' # protein ID
df_43_rfl['Luciferase'] = 'Q8ZQQ2'
df_43_rfl['N-Terminus'] = 'SlrP'
df_43_rfl['Normalisierung']='PC5'

df_44_rfl = pd.concat([s65])
df_44_rfl['number'] = np.arange(0,8,1)
df_44_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA044'
df_44_rfl['Konstrukt'] = 'Q8ZNR3' # protein ID
df_44_rfl['Luciferase'] = 'RFL'
df_44_rfl['N-Terminus'] = 'SopA'
df_44_rfl['Normalisierung']='PC6'

df_48_rfl = pd.concat([s66])
df_48_rfl['number'] = np.arange(0,8,1)
df_48_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA048'
df_48_rfl['Konstrukt'] = 'O52623' # protein ID
df_48_rfl['Luciferase'] = 'RFL'
df_48_rfl['N-Terminus'] = 'SopE'
df_48_rfl['Normalisierung']='PC6'

# nanoluc luciferase plate 4 pJA058,065, 066

df_p23 = pd.read_csv(file23, skiprows=1)
df_p23 = df_p23.iloc[4:12,:12]
df_p23.columns = ['L1','L2','L3']*4
df_p23 = df_p23.reset_index(drop=True)

s67 = df_p23.iloc[:8,:3]
s68 = df_p23.iloc[:8,4:7]
s69 = df_p23.iloc[:8,8:11]

df_58_nl = pd.concat([s67])
df_58_nl['number'] = np.arange(0,8,1)
df_58_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA058'
df_58_nl['Konstrukt'] = 'Q7BVH7' # protein ID
df_58_nl['Luciferase'] = 'NL'
df_58_nl['N-Terminus'] = 'SseB'
df_58_nl['Normalisierung']='PC6'

df_65_nl = pd.concat([s68])
df_65_nl['number'] = np.arange(0,8,1)
df_65_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA065'
df_65_nl['Konstrukt'] = 'Q9L9J3' # protein ID
df_65_nl['Luciferase'] = 'NL'
df_65_nl['N-Terminus'] = 'SseK1'
df_65_nl['Normalisierung']='PC6'

df_66_nl = pd.concat([s69])
df_66_nl['number'] = np.arange(0,8,1)
df_66_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA066'
df_66_nl['Konstrukt'] = 'Q8ZNP4' # protein ID
df_66_nl['Luciferase'] = 'NL'
df_66_nl['N-Terminus'] = 'SseK2'
df_66_nl['Normalisierung']='PC6'

# red firefly luciferase plate 4 pJA058,065, 066

df_p24 = pd.read_csv(file24, skiprows=1)
df_p24 = df_p24.iloc[4:12,:12]
df_p24.columns = ['L1','L2','L3']*4
df_p24 = df_p24.reset_index(drop=True)

s70 = df_p24.iloc[:8,:3]
s71 = df_p24.iloc[:8,4:7]
s72 = df_p24.iloc[:8,8:11]

df_58_rfl = pd.concat([s70])
df_58_rfl['number'] = np.arange(0,8,1)
df_58_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA058'
df_58_rfl['Konstrukt'] = 'Q7BVH7' # protein ID
df_58_rfl['Luciferase'] = 'RFL'
df_58_rfl['N-Terminus'] = 'SseB'
df_58_rfl['Normalisierung']='PC6'

df_65_rfl = pd.concat([s71])
df_65_rfl['number'] = np.arange(0,8,1)
df_65_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA065'
df_65_rfl['Konstrukt'] = 'Q9L9J3' # protein ID
df_65_rfl['Luciferase'] = 'RFL'
df_65_rfl['N-Terminus'] = 'SseK1'
df_65_rfl['Normalisierung']='PC6'

df_66_rfl = pd.concat([s72])
df_66_rfl['number'] = np.arange(0,8,1)
df_66_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA066'
df_66_rfl['Konstrukt'] = 'Q8ZNP4' # protein ID
df_66_rfl['Luciferase'] = 'RFL'
df_66_rfl['N-Terminus'] = 'SseK2'
df_66_rfl['Normalisierung']='PC6'


###########
#Assay 5.10.22
###########

file25 = '221005_plate01_nl_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pJA023
file26 = '221005_plate01_rfl_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pJA023
file27 = '221005_plate02_nl_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pJA024
file28 = '221005_plate02_rfl_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pJA024
file29 = '221005_plate03_nl_001.csv' # nanoluc luciferase plate 3 pJA025, 027, 030
file30 = '221005_plate03_rfl_001.csv' # red firefly luciferase plate 3 pJA025, 027, 030

###Nanluc luciferase plate 1: positivkontrollen + pja023

df_p25 = pd.read_csv(file25, skiprows=1)
df_p25 = df_p25.iloc[4:12,:12]
df_p25.columns = ['L1','L2','L3']*4
df_p25 = df_p25.reset_index(drop=True)

s73 = df_p25.iloc[:8,:3]
s74 = df_p25.iloc[:8,4:7]
s75 = df_p25.iloc[:8,8:11]

df_pc7_nl = pd.concat([s73])
df_pc7_nl['number'] = np.arange(0,8,1)
df_pc7_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.7'
df_pc7_nl['Konstrukt'] = 'Positivkontrolle.7'
df_pc7_nl['Luciferase'] = 'NL'
df_pc7_nl['N-Terminus'] = 'PC7 SptP'
df_pc7_nl['Normalisierung'] = 'PC7'

df_pc8_nl = pd.concat([s74])
df_pc8_nl['number'] = np.arange(0,8,1)
df_pc8_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.8'
df_pc8_nl['Konstrukt'] = 'Positivkontrolle.8'
df_pc8_nl['Luciferase'] = 'NL'
df_pc8_nl['N-Terminus'] = 'PC8 SptP'
df_pc8_nl['Normalisierung'] = 'PC8'

df_23_nl = pd.concat([s75])
df_23_nl['number'] = np.arange(0,8,1)
df_23_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA023'
df_23_nl['Konstrukt'] = 'A0A0F6B537' 
df_23_nl['Luciferase'] = 'NL'
df_23_nl['N-Terminus'] = 'GogA' 
df_23_nl['Normalisierung'] = 'PC7'

# red firefly luciferase plate 1 positiv Kontrollen + pja023

df_p26 = pd.read_csv(file26, skiprows=1)
df_p26 = df_p26.iloc[4:12,:12]
df_p26.columns = ['L1','L2','L3']*4
df_p26 = df_p26.reset_index(drop=True)

s76 = df_p26.iloc[:8,:3]
s77 = df_p26.iloc[:8,4:7]
s78 = df_p26.iloc[:8,8:11]

df_pc7_rfl = pd.concat([s76])
df_pc7_rfl['number'] = np.arange(0,8,1)
df_pc7_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.7'
df_pc7_rfl['Konstrukt'] = 'Positivkontrolle.7'
df_pc7_rfl['Luciferase'] = 'RFL'
df_pc7_rfl['N-Terminus'] = 'PC7 SptP'
df_pc7_rfl['Normalisierung'] = 'PC7'

df_pc8_rfl = pd.concat([s77])
df_pc8_rfl['number'] = np.arange(0,8,1)
df_pc8_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.8'
df_pc8_rfl['Konstrukt'] = 'Positivkontrolle.8'
df_pc8_rfl['Luciferase'] = 'RFL'
df_pc8_rfl['N-Terminus'] = 'PC8 SptP'
df_pc8_rfl['Normalisierung'] = 'PC8'

df_23_rfl = pd.concat([s78])
df_23_rfl['number'] = np.arange(0,8,1)
df_23_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA023'
df_23_rfl['Konstrukt'] = 'A0A0F6B537' 
df_23_rfl['Luciferase'] = 'RFL'
df_23_rfl['N-Terminus'] = 'GogA' 
df_23_rfl['Normalisierung'] = 'PC7'

# nanoluc luciferase plate 2 negativ Kontrollen + pja024

df_p27 = pd.read_csv(file27, skiprows=1)
df_p27 = df_p27.iloc[4:12,:12]
df_p27.columns = ['L1','L2','L3']*4
df_p27 = df_p27.reset_index(drop=True)

s79 = df_p27.iloc[:8,:3]
s80 = df_p27.iloc[:8,4:7]
s81 = df_p27.iloc[:8,8:11]

df_pc_leer4_nl = pd.concat([s79])
df_pc_leer4_nl['number'] = np.arange(0,8,1)
df_pc_leer4_nl['Linie'] = 'SB905 dSipA dSptP pHilA.4'
df_pc_leer4_nl['Konstrukt'] = 'Leerkontrolle4'
df_pc_leer4_nl['Luciferase'] = 'NL'
df_pc_leer4_nl['N-Terminus'] = 'leer4'
df_pc_leer4_nl['Normalisierung'] = 'PC7'

df_nc4_nl = pd.concat([s80])
df_nc4_nl['number'] = np.arange(0,8,1)
df_nc4_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.4'
df_nc4_nl['Konstrukt'] = 'Negativkontrolle.4'
df_nc4_nl['Luciferase'] = 'NL'
df_nc4_nl['N-Terminus'] = 'NC4 SptP'
df_nc4_nl['Normalisierung'] = 'PC7'

df_24_nl = pd.concat([s81])
df_24_nl['number'] = np.arange(0,8,1)
df_24_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA024'
df_24_nl['Konstrukt'] = 'Q8ZN18' ### name Linie=protein ID
df_24_nl['Luciferase'] = 'NL'
df_24_nl['N-Terminus'] = 'GogB'
df_24_nl['Normalisierung']='PC7'

# red firefly luciferase plate 2 negativ Kontrollen + pja024

df_p28 = pd.read_csv(file28, skiprows=1)
df_p28 = df_p28.iloc[4:12,:12]
df_p28.columns = ['L1','L2','L3']*4
df_p28 = df_p28.reset_index(drop=True)

s82 = df_p28.iloc[:8,:3]
s83 = df_p28.iloc[:8,4:7]
s84 = df_p28.iloc[:8,8:11]
  
df_pc_leer4_rfl = pd.concat([s82])
df_pc_leer4_rfl['number'] = np.arange(0,8,1)
df_pc_leer4_rfl['Linie'] = 'SB905 dSipA dSptP pHilA.4'
df_pc_leer4_rfl['Konstrukt'] = 'Leerkontrolle4'
df_pc_leer4_rfl['Luciferase'] = 'RFL'
df_pc_leer4_rfl['N-Terminus'] = 'leer4'
df_pc_leer4_rfl['Normalisierung'] = 'PC7'

df_nc4_rfl = pd.concat([s83])
df_nc4_rfl['number'] = np.arange(0,8,1)
df_nc4_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.4'
df_nc4_rfl['Konstrukt'] = 'Negativkontrolle.4'
df_nc4_rfl['Luciferase'] = 'RFL'
df_nc4_rfl['N-Terminus'] = 'NC4 SptP'
df_nc4_rfl['Normalisierung'] = 'PC7'

df_24_rfl = pd.concat([s84])
df_24_rfl['number'] = np.arange(0,8,1)
df_24_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA024'
df_24_rfl['Konstrukt'] = 'Q8ZN18' ### name Linie=protein ID
df_24_rfl['Luciferase'] = 'NL'
df_24_rfl['N-Terminus'] = 'GogB'
df_24_rfl['Normalisierung']='PC7'

# nanoluc luciferase plate 3 pJA0025, 027, 030

df_p29 = pd.read_csv(file29, skiprows=1)
df_p29 = df_p29.iloc[4:12,:12]
df_p29.columns = ['L1','L2','L3']*4
df_p29 = df_p29.reset_index(drop=True)

s85 = df_p29.iloc[:8,:3]
s101 = df_p29.iloc[:8,4:7]
s86 = df_p29.iloc[:8,8:11]

df_25_nl = pd.concat([s85])
df_25_nl['number'] = np.arange(0,8,1)
df_25_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA025'
df_25_nl['Konstrukt'] = 'A0A0F6AZI6' ### name category=protein ID
df_25_nl['Luciferase'] = 'NL'
df_25_nl['N-Terminus'] = 'GtgA'
df_25_nl['Normalisierung'] = 'PC7'

df_27_nl = pd.concat([s101])
df_27_nl['number'] = np.arange(0,8,1)
df_27_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA027'
df_27_nl['Konstrukt'] = 'A0A0H3NF83' # protein ID
df_27_nl['Luciferase'] = 'NL'
df_27_nl['N-Terminus'] = 'OrgC'
df_27_nl['Normalisierung']='PC8'

df_30_nl = pd.concat([s86])
df_30_nl['number'] = np.arange(0,8,1)
df_30_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA030'
df_30_nl['Konstrukt'] = 'Q8ZMM8' ### name category=protein ID
df_30_nl['Luciferase'] = 'NL'
df_30_nl['N-Terminus'] = 'PipB2'
df_30_nl['Normalisierung'] = 'PC8'

# red firefly luciferase plate 3 pJA023, 024, 025, 030

df_p30 = pd.read_csv(file30, skiprows=1)
df_p30 = df_p30.iloc[4:12,:12]
df_p30.columns = ['L1','L2','L3']*4
df_p30 = df_p30.reset_index(drop=True)

s87 = df_p30.iloc[:8,:3]
s104 = df_p30.iloc[:8,4:7]
s88 = df_p30.iloc[:8,8:11]

df_25_rfl = pd.concat([s87])
df_25_rfl['number'] = np.arange(0,8,1)
df_25_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA025'
df_25_rfl['Konstrukt'] = 'A0A0F6AZI6' ### name category=protein ID
df_25_rfl['Luciferase'] = 'RFL'
df_25_rfl['N-Terminus'] = 'GtgA'
df_25_rfl['Normalisierung'] = 'PC7'

df_27_rfl = pd.concat([s104])
df_27_rfl['number'] = np.arange(0,8,1)
df_27_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA027'
df_27_rfl['Konstrukt'] = 'A0A0H3NF83' # protein ID
df_27_rfl['Luciferase'] = 'RFL'
df_27_rfl['N-Terminus'] = 'OrgC'
df_27_rfl['Normalisierung']='PC8'

df_30_rfl = pd.concat([s88])
df_30_rfl['number'] = np.arange(0,8,1)
df_30_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA030'
df_30_rfl['Konstrukt'] = 'Q8ZMM8' ### name category=protein ID
df_30_rfl['Luciferase'] = 'RFL'
df_30_rfl['N-Terminus'] = 'PipB2'
df_30_rfl['Normalisierung'] = 'PC8'


######
#Assay 20.09.22
######

file31 = '220920_nl_plate_01_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pJA0222
file32 = '220920_rfl_plate_01_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen + pJA022
file33 = '220920_nl_plate_02_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pJA026
file34 = '220920_rfl_plate_02_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pJA026
file35 = '220920_nl_plate_03_001.csv' # nanoluc luciferase plate 3 pJA027, 028, 041
file36 = '220920_rfl_plate_03_001.csv' # red firefly luciferase plate 3 pJA027, 028, 041
file37 = '220920_nl_plate_04_001.csv' # nanoluc luciferase plate 4 pJA042, 045, 046
file38 = '220920_rfl_plate_04_001.csv' # red firefly luciferase plate 4 pJA042, 045, 046

# nanoluc luciferase plate 1 positiv Kontrollen + pJA022

df_p31 = pd.read_csv(file31, skiprows=1)
df_p31 = df_p31.iloc[4:12,:12]
df_p31.columns = ['L1','L2','L3']*4
df_p31 = df_p31.reset_index(drop=True)

s89 = df_p31.iloc[:8,:3]
s90 = df_p31.iloc[:8,4:7]
s91 = df_p31.iloc[:8,8:11]

df_pc9_nl = pd.concat([s89])
df_pc9_nl['number'] = np.arange(0,8,1)
df_pc9_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.9'
df_pc9_nl['Konstrukt'] = 'Positivkontrolle.9'
df_pc9_nl['Luciferase'] = 'NL'
df_pc9_nl['N-Terminus'] = 'PC9 SptP'
df_pc9_nl['Normalisierung']='PC9'

df_pc10_nl = pd.concat([s90])
df_pc10_nl['number'] = np.arange(0,8,1)
df_pc10_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.10'
df_pc10_nl['Konstrukt'] = 'Positivkontrolle.10'
df_pc10_nl['Luciferase'] = 'NL'
df_pc10_nl['N-Terminus'] = 'PC10 SptP'
df_pc10_nl['Normalisierung']='PC10'

df_22_nl = pd.concat([s91])
df_22_nl['number'] = np.arange(0,8,1)
df_22_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA022'
df_22_nl['Konstrukt'] = 'Q8ZMI3' # protein ID
df_22_nl['Luciferase'] = 'NL'
df_22_nl['N-Terminus'] = 'AvrA'
df_22_nl['Normalisierung']='PC9'

# red firefly luciferase plate 1 positiv Kontrollen + pMP053

df_p32 = pd.read_csv(file32, skiprows=1)
df_p32 = df_p32.iloc[4:12,:12]
df_p32.columns = ['L1','L2','L3']*4
df_p32 = df_p32.reset_index(drop=True)

s92 = df_p32.iloc[:8,:3]
s93= df_p32.iloc[:8,4:7]
s94 = df_p32.iloc[:8,8:11]
  
df_pc9_rfl = pd.concat([s92])
df_pc9_rfl['number'] = np.arange(0,8,1)
df_pc9_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.9'
df_pc9_rfl['Konstrukt'] = 'Positivkontrolle.9'
df_pc9_rfl['Luciferase'] = 'RFL'
df_pc9_rfl['N-Terminus'] = 'PC9 SptP'
df_pc9_rfl['Normalisierung']='PC9'

df_pc10_rfl = pd.concat([s93])
df_pc10_rfl['number'] = np.arange(0,8,1)
df_pc10_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.10'
df_pc10_rfl['Konstrukt'] = 'Positivkontrolle.10'
df_pc10_rfl['Luciferase'] = 'RFL'
df_pc10_rfl['N-Terminus'] = 'PC10 SptP'
df_pc10_rfl['Normalisierung']='PC10'

df_22_rfl = pd.concat([s94])
df_22_rfl['number'] = np.arange(0,8,1)
df_22_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA022'
df_22_rfl['Konstrukt'] = 'Q8ZMI3' # protein ID
df_22_rfl['Luciferase'] = 'RFL'
df_22_rfl['N-Terminus'] = 'AvrA'
df_22_rfl['Normalisierung']='PC9'

# nanoluc luciferase plate 2 negativ Kontrollen + pJA026

df_p33 = pd.read_csv(file33, skiprows=1)
df_p33 = df_p33.iloc[4:12,:12]
df_p33.columns = ['L1','L2','L3']*4
df_p33 = df_p33.reset_index(drop=True)

s95 = df_p33.iloc[:8,:3]
s96 = df_p33.iloc[:8,4:7]
s97 = df_p33.iloc[:8,8:11]

df_pc_leer5_nl = pd.concat([s95])
df_pc_leer5_nl['number'] = np.arange(0,8,1)
df_pc_leer5_nl['Linie'] = 'SB905 dSipA dSptP pHilA.5'
df_pc_leer5_nl['Konstrukt'] = 'Leerkontrolle5'
df_pc_leer5_nl['Luciferase'] = 'NL'
df_pc_leer5_nl['N-Terminus'] = 'leer5'
df_pc_leer5_nl['Normalisierung']='PC9'

df_nc5_nl = pd.concat([s96])
df_nc5_nl['number'] = np.arange(0,8,1)
df_nc5_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.5'
df_nc5_nl['Konstrukt'] = 'Negativkontrolle.5'
df_nc5_nl['Luciferase'] = 'NL'
df_nc5_nl['N-Terminus'] = 'NC5 SptP'
df_nc5_nl['Normalisierung']='PC9'

df_26_nl = pd.concat([s97])
df_26_nl['number'] = np.arange(0,8,1)
df_26_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA026'
df_26_nl['Konstrukt'] = 'A0A0H3N9Y3' # protein ID
df_26_nl['Luciferase'] = 'NL'
df_26_nl['N-Terminus'] = 'GtgE'
df_26_nl['Normalisierung']='PC9'

# red firefly luciferase plate 2 negativ Kontrollen + pJA026

df_p34 = pd.read_csv(file34, skiprows=1)
df_p34 = df_p34.iloc[4:12,:12]
df_p34.columns = ['L1','L2','L3']*4
df_p34 = df_p34.reset_index(drop=True)

s98 = df_p34.iloc[:8,:3]
s99 = df_p34.iloc[:8,4:7]
s100 = df_p34.iloc[:8,8:11]

df_pc_leer5_rfl = pd.concat([s98])
df_pc_leer5_rfl['number'] = np.arange(0,8,1)
df_pc_leer5_rfl['Linie'] = 'SB905 dSipA dSptP pHilA.5'
df_pc_leer5_rfl['Konstrukt'] = 'Leerkontrolle5'
df_pc_leer5_rfl['Luciferase'] = 'RFL'
df_pc_leer5_rfl['N-Terminus'] = 'leer5'
df_pc_leer5_rfl['Normalisierung']='PC9'

df_nc5_rfl = pd.concat([s99])
df_nc5_rfl['number'] = np.arange(0,8,1)
df_nc5_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.5'
df_nc5_rfl['Konstrukt'] = 'Negativkontrolle.5'
df_nc5_rfl['Luciferase'] = 'RFL'
df_nc5_rfl['N-Terminus'] = 'NC5 SptP'
df_nc5_rfl['Normalisierung']='PC9'

df_26_rfl = pd.concat([s100])
df_26_rfl['number'] = np.arange(0,8,1)
df_26_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA026'
df_26_rfl['Konstrukt'] = 'A0A0H3N9Y3' # protein ID
df_26_rfl['Luciferase'] = 'RFL'
df_26_rfl['N-Terminus'] = 'GtgE'
df_26_rfl['Normalisierung']='PC9'

# nanoluc luciferase plate 3 pJA027, 028, 041

df_p35 = pd.read_csv(file35, skiprows=1)
df_p35 = df_p35.iloc[4:12,:12]
df_p35.columns = ['L1','L2','L3']*4
df_p35 = df_p35.reset_index(drop=True)

s102 = df_p35.iloc[:8,4:7]
s103 = df_p35.iloc[:8,8:11]

df_28_nl = pd.concat([s102])
df_28_nl['number'] = np.arange(0,8,1)
df_28_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA028'
df_28_nl['Konstrukt'] = 'A0A0F6AZQ0' # protein ID
df_28_nl['Luciferase'] = 'NL'
df_28_nl['N-Terminus'] = 'PipA'
df_28_nl['Normalisierung']='PC10'

df_41_nl = pd.concat([s103])
df_41_nl['number'] = np.arange(0,8,1)
df_41_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA041'
df_41_nl['Konstrukt'] = 'P0CL47' # protein ID
df_41_nl['Luciferase'] = 'NL'
df_41_nl['N-Terminus'] = 'SipC'
df_41_nl['Normalisierung']='PC10'

# red firefly luciferase plate 3 pJA027, 028, 041

df_p36 = pd.read_csv(file36, skiprows=1)
df_p36 = df_p36.iloc[4:12,:12]
df_p36.columns = ['L1','L2','L3']*4
df_p36 = df_p36.reset_index(drop=True)

s105 = df_p36.iloc[:8,4:7]
s106 = df_p36.iloc[:8,8:11]

df_28_rfl = pd.concat([s105])
df_28_rfl['number'] = np.arange(0,8,1)
df_28_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA028'
df_28_rfl['Konstrukt'] = 'A0A0F6AZQ0' # protein ID
df_28_rfl['Luciferase'] = 'RFL'
df_28_rfl['N-Terminus'] = 'PipA'
df_28_rfl['Normalisierung']='PC10'

df_41_rfl = pd.concat([s106])
df_41_rfl['number'] = np.arange(0,8,1)
df_41_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA041'
df_41_rfl['Konstrukt'] = 'P0CL47' # protein ID
df_41_rfl['Luciferase'] = 'RFL'
df_41_rfl['N-Terminus'] = 'SipC'
df_41_rfl['Normalisierung']='PC10'

# nanoluc luciferase plate 4 pJA042, 045, 046

df_p37 = pd.read_csv(file37, skiprows=1)
df_p37 = df_p37.iloc[4:12,:12]
df_p37.columns = ['L1','L2','L3']*4
df_p37 = df_p37.reset_index(drop=True)

s107 = df_p37.iloc[:8,:3]
s108 = df_p37.iloc[:8,4:7]
s109 = df_p37.iloc[:8,8:11]

df_42_nl = pd.concat([s107])
df_42_nl['number'] = np.arange(0,8,1)
df_42_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA042'
df_42_nl['Konstrukt'] = 'Q56026' # protein ID
df_42_nl['Luciferase'] = 'NL'
df_42_nl['N-Terminus'] = 'SipD'
df_42_nl['Normalisierung']='PC10'

df_45_nl = pd.concat([s108])
df_45_nl['number'] = np.arange(0,8,1)
df_45_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA045'
df_45_nl['Konstrukt'] = 'O30916' # protein ID
df_45_nl['Luciferase'] = 'NL'
df_45_nl['N-Terminus'] = 'SopB'
df_45_nl['Normalisierung']='PC10'

df_46_nl = pd.concat([s109])
df_46_nl['number'] = np.arange(0,8,1)
df_46_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA046'
df_46_nl['Konstrukt'] = 'P40722' # protein ID
df_46_nl['Luciferase'] = 'NL'
df_46_nl['N-Terminus'] = 'SopD'
df_46_nl['Normalisierung']='PC10'

# red firefly luciferase plate 4 pJA042, 045, 046

df_p38 = pd.read_csv(file38, skiprows=1)
df_p38 = df_p38.iloc[4:12,:12]
df_p38.columns = ['L1','L2','L3']*4
df_p38 = df_p38.reset_index(drop=True)

s110 = df_p38.iloc[:8,:3]
s111 = df_p38.iloc[:8,4:7]
s112 = df_p38.iloc[:8,8:11]

df_42_rfl = pd.concat([s110])
df_42_rfl['number'] = np.arange(0,8,1)
df_42_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA042'
df_42_rfl['Konstrukt'] = 'Q56026' # protein ID
df_42_rfl['Luciferase'] = 'RFL'
df_42_rfl['N-Terminus'] = 'SipD'
df_42_rfl['Normalisierung']='PC10'

df_45_rfl = pd.concat([s111])
df_45_rfl['number'] = np.arange(0,8,1)
df_45_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA045'
df_45_rfl['Konstrukt'] = 'O30916' # protein ID
df_45_rfl['Luciferase'] = 'RFL'
df_45_rfl['N-Terminus'] = 'SopB'
df_45_rfl['Normalisierung']='PC10'

df_46_rfl = pd.concat([s112])
df_46_rfl['number'] = np.arange(0,8,1)
df_46_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA046'
df_46_rfl['Konstrukt'] = 'P40722' # protein ID
df_46_rfl['Luciferase'] = 'RFL'
df_46_rfl['N-Terminus'] = 'SopD'
df_46_rfl['Normalisierung']='PC10'


########
#Assay 21.9.22
########

file39 = '220921_nl_plate01_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen 
file40 = '220921_rfl_plate01_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen
file41 = '220921_nl_plate02_001.csv' # nanoluc luciferase plate 2 negativ kontrollen
file42 = '220921_rfl_plate02_001.csv' # red firefly luciferase plate 2 negativ kontrollen
file43 = '220921_nl_plate03_001.csv' # nanoluc luciferase plate 3 pJA049, 050, 053
file44 = '220921_rfl_plate03_001.csv' # red firefly luciferase plate 3 pJA049, 050, 053
file45 = '220921_nl_plate04_001.csv' # nanoluc luciferase plate 4 pJA056, 057, 061
file46 = '220921_rfl_plate04_001.csv' # red firefly luciferase plate 4 pJA056, 057, 061
file47 = '220921_nl_plate05_001.csv' # nanoluc luciferase plate 4 pJA063, 064
file48 = '220921_rfl_plate05_001.csv' # red firefly luciferase plate 4 pJA063, 064

# nanoluc luciferase plate 1 positiv Kontrollen

df_p39 = pd.read_csv(file39, skiprows=1)
df_p39 = df_p39.iloc[4:12,:12]
df_p39.columns = ['L1','L2','L3']*4
df_p39 = df_p39.reset_index(drop=True)

s113 = df_p39.iloc[:8,:3]
s114 = df_p39.iloc[:8,4:7]
s115 = df_p39.iloc[:8,8:11]

df_pc11_nl = pd.concat([s113])
df_pc11_nl['number'] = np.arange(0,8,1)
df_pc11_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.11'
df_pc11_nl['Konstrukt'] = 'PC.11'
df_pc11_nl['Luciferase'] = 'NL'
df_pc11_nl['N-Terminus'] = 'PC11 SptP'
df_pc11_nl['Normalisierung']='PC11'

df_pc12_nl = pd.concat([s114])
df_pc12_nl['number'] = np.arange(0,8,1)
df_pc12_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.12'
df_pc12_nl['Konstrukt'] = 'PC.12'
df_pc12_nl['Luciferase'] = 'NL'
df_pc12_nl['N-Terminus'] = 'PC12 SptP'
df_pc12_nl['Normalisierung']='PC12'

df_pc13_nl = pd.concat([s115])
df_pc13_nl['number'] = np.arange(0,8,1)
df_pc13_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.13'
df_pc13_nl['Konstrukt'] = 'PC.13' 
df_pc13_nl['Luciferase'] = 'NL'
df_pc13_nl['N-Terminus'] = 'PC13 SptP'
df_pc13_nl['Normalisierung']='PC13'

# red firefly luciferase plate 1 positiv Kontrollen

df_p40 = pd.read_csv(file40, skiprows=1)
df_p40 = df_p40.iloc[4:12,:12]
df_p40.columns = ['L1','L2','L3']*4
df_p40 = df_p40.reset_index(drop=True)

s116 = df_p40.iloc[:8,:3]
s117 = df_p40.iloc[:8,4:7]
s118 = df_p40.iloc[:8,8:11]
  
df_pc11_rfl = pd.concat([s116])
df_pc11_rfl['number'] = np.arange(0,8,1)
df_pc11_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.11'
df_pc11_rfl['Konstrukt'] = 'PC.11'
df_pc11_rfl['Luciferase'] = 'RFL'
df_pc11_rfl['N-Terminus'] = 'PC11 SptP'
df_pc11_rfl['Normalisierung']='PC11'

df_pc12_rfl = pd.concat([s117])
df_pc12_rfl['number'] = np.arange(0,8,1)
df_pc12_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.12'
df_pc12_rfl['Konstrukt'] = 'PC.12'
df_pc12_rfl['Luciferase'] = 'RFL'
df_pc12_rfl['N-Terminus'] = 'PC12 SptP'
df_pc12_rfl['Normalisierung']='PC12'

df_pc13_rfl = pd.concat([s118])
df_pc13_rfl['number'] = np.arange(0,8,1)
df_pc13_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.13'
df_pc13_rfl['Konstrukt'] = 'PC.13' 
df_pc13_rfl['Luciferase'] = 'RFL'
df_pc13_rfl['N-Terminus'] = 'PC13 SptP'
df_pc13_rfl['Normalisierung']='PC13'

# nanoluc luciferase plate 2 negativ Kontrollen + pJA047

df_p41 = pd.read_csv(file41, skiprows=1)
df_p41 = df_p41.iloc[4:12,:12]
df_p41.columns = ['L1','L2','L3']*4
df_p41 = df_p41.reset_index(drop=True)

s119 = df_p41.iloc[:8,:3]
s120 = df_p41.iloc[:8,4:7]

df_pc_leer6_nl = pd.concat([s119])
df_pc_leer6_nl['number'] = np.arange(0,8,1)
df_pc_leer6_nl['Linie'] = 'SB905 dSipA dSptP pHilA.6'
df_pc_leer6_nl['Konstrukt'] = 'Leerkontrolle6'
df_pc_leer6_nl['Luciferase'] = 'NL'
df_pc_leer6_nl['N-Terminus'] = 'leer6'
df_pc_leer6_nl['Normalisierung']='PC11'

df_nc6_nl = pd.concat([s120])
df_nc6_nl['number'] = np.arange(0,8,1)
df_nc6_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.6'
df_nc6_nl['Konstrukt'] = 'NC.6'
df_nc6_nl['Luciferase'] = 'NL'
df_nc6_nl['N-Terminus'] = 'NC6 SptP'
df_nc6_nl['Normalisierung']='PC11'

# red firefly luciferase plate 2 negativ Kontrollen + pJA047

df_p42 = pd.read_csv(file42, skiprows=1)
df_p42 = df_p42.iloc[4:12,:12]
df_p42.columns = ['L1','L2','L3']*4
df_p42 = df_p42.reset_index(drop=True)

s122 = df_p42.iloc[:8,:3]
s123 = df_p42.iloc[:8,4:7]

df_pc_leer6_rfl = pd.concat([s122])
df_pc_leer6_rfl['number'] = np.arange(0,8,1)
df_pc_leer6_rfl['Linie'] = 'SB905 dSipA dSptP pHilA.6'
df_pc_leer6_rfl['Konstrukt'] = 'Leerkontrolle6'
df_pc_leer6_rfl['Luciferase'] = 'RFL'
df_pc_leer6_rfl['N-Terminus'] = 'leer6'
df_pc_leer6_rfl['Normalisierung']='PC11'

df_nc6_rfl = pd.concat([s123])
df_nc6_rfl['number'] = np.arange(0,8,1)
df_nc6_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049.6'
df_nc6_rfl['Konstrukt'] = 'NC.6'
df_nc6_rfl['Luciferase'] = 'RFL'
df_nc6_rfl['N-Terminus'] = 'NC6 SptP'
df_nc6_rfl['Normalisierung']='PC11'

# nanoluc luciferase plate 3 pJA049, 050, 053

df_p43 = pd.read_csv(file43, skiprows=1)
df_p43 = df_p43.iloc[4:12,:12]
df_p43.columns = ['L1','L2','L3']*4
df_p43 = df_p43.reset_index(drop=True)

s125 = df_p43.iloc[:8,:3]
s126 = df_p43.iloc[:8,4:7]
s127 = df_p43.iloc[:8,8:11]

df_49_nl = pd.concat([s125])
df_49_nl['number'] = np.arange(0,8,1)
df_49_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA049'
df_49_nl['Konstrukt'] = 'Q7CQD4' ### name Linie=protein ID
df_49_nl['Luciferase'] = 'NL'
df_49_nl['N-Terminus'] = 'SopE2'
df_49_nl['Normalisierung']='PC11'

df_50_nl = pd.concat([s126])
df_50_nl['number'] = np.arange(0,8,1)
df_50_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA050'
df_50_nl['Konstrukt'] = 'Q8ZPY9' ### name Linie=protein ID
df_50_nl['Luciferase'] = 'NL'
df_50_nl['N-Terminus'] = 'SopF'
df_50_nl['Normalisierung']='PC11'

df_53_nl = pd.concat([s127])
df_53_nl['number'] = np.arange(0,8,1)
df_53_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA053'
df_53_nl['Konstrukt'] = 'P0A2N2' ### name Linie=protein ID
df_53_nl['Luciferase'] = 'NL'
df_53_nl['N-Terminus'] = 'SpvD'
df_53_nl['Normalisierung']='PC12'

# red firefly luciferase plate 3 pJA049, 050, 053

df_p44 = pd.read_csv(file44, skiprows=1)
df_p44 = df_p44.iloc[4:12,:12]
df_p44.columns = ['L1','L2','L3']*4
df_p44 = df_p44.reset_index(drop=True)

s128 = df_p44.iloc[:8,:3]
s129 = df_p44.iloc[:8,4:7]
s130 = df_p44.iloc[:8,8:11]

df_49_rfl = pd.concat([s128])
df_49_rfl['number'] = np.arange(0,8,1)
df_49_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA049'
df_49_rfl['Konstrukt'] = 'Q7CQD4' ### name Linie=protein ID
df_49_rfl['Luciferase'] = 'RFL'
df_49_rfl['N-Terminus'] = 'SopE2'
df_49_rfl['Normalisierung']='PC11'

df_50_rfl = pd.concat([s129])
df_50_rfl['number'] = np.arange(0,8,1)
df_50_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA050'
df_50_rfl['Konstrukt'] = 'Q8ZPY9' ### name Linie=protein ID
df_50_rfl['Luciferase'] = 'RFL'
df_50_rfl['N-Terminus'] = 'SopF'
df_50_rfl['Normalisierung']='PC11'

df_53_rfl = pd.concat([s130])
df_53_rfl['number'] = np.arange(0,8,1)
df_53_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA053'
df_53_rfl['Konstrukt'] = 'P0A2N2' ### name Linie=protein ID
df_53_rfl['Luciferase'] = 'RFL'
df_53_rfl['N-Terminus'] = 'SpvD'
df_53_rfl['Normalisierung']='PC12'

# nanoluc luciferase plate 4 pJA056, 057, 061

df_p45 = pd.read_csv(file45, skiprows=1)
df_p45 = df_p45.iloc[4:12,:12]
df_p45.columns = ['L1','L2','L3']*4
df_p45= df_p45.reset_index(drop=True)

s132 = df_p45.iloc[:8,4:7]

df_57_nl = pd.concat([s132])
df_57_nl['number'] = np.arange(0,8,1)
df_57_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA057'
df_57_nl['Konstrukt'] = 'H9L496' ### name Linie=protein ID
df_57_nl['Luciferase'] = 'NL'
df_57_nl['N-Terminus'] = 'SsaL'
df_57_nl['Normalisierung']='PC12'

# red firefly luciferase plate 4 pJA056, 057, 061

df_p46 = pd.read_csv(file46, skiprows=1)
df_p46 = df_p46.iloc[4:12,:12]
df_p46.columns = ['L1','L2','L3']*4
df_p46 = df_p46.reset_index(drop=True)

s135 = df_p46.iloc[:8,4:7]

df_57_rfl = pd.concat([s135])
df_57_rfl['number'] = np.arange(0,8,1)
df_57_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA057'
df_57_rfl['Konstrukt'] = 'H9L496' ### name Linie=protein ID
df_57_rfl['Luciferase'] = 'RFL'
df_57_rfl['N-Terminus'] = 'SsaL'
df_57_rfl['Normalisierung']='PC12'

# nanoluc luciferase plate 5 pJA063, 064

df_p47 = pd.read_csv(file47, skiprows=1)
df_p47 = df_p47.iloc[4:12,:12]
df_p47.columns = ['L1','L2','L3']*4
df_p47 = df_p47.reset_index(drop=True)

s137 = df_p47.iloc[:8,:3]
s138 = df_p47.iloc[:8,4:7]

df_63_nl = pd.concat([s137])
df_63_nl['number'] = np.arange(0,8,1)
df_63_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA063'
df_63_nl['Konstrukt'] = 'A0A0F6AZL3' ### name Linie=protein ID
df_63_nl['Luciferase'] = 'NL'
df_63_nl['N-Terminus'] = 'SseI'
df_63_nl['Normalisierung']='PC12'

df_64_nl = pd.concat([s138])
df_64_nl['number'] = np.arange(0,8,1)
df_64_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA064'
df_64_nl['Konstrukt'] = 'Q9FD10' ### name Linie=protein ID
df_64_nl['Luciferase'] = 'NL'
df_64_nl['N-Terminus'] = 'SseJ'
df_64_nl['Normalisierung']='PC13'

# red firefly luciferase plate 5 pJA063, 064

df_p48 = pd.read_csv(file48, skiprows=1)
df_p48 = df_p48.iloc[4:12,:12]
df_p48.columns = ['L1','L2','L3']*4
df_p48 = df_p48.reset_index(drop=True)

s139 = df_p48.iloc[:8,:3]
s140 = df_p48.iloc[:8,4:7]

df_63_rfl = pd.concat([s139])
df_63_rfl['number'] = np.arange(0,8,1)
df_63_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA063'
df_63_rfl['Konstrukt'] = 'A0A0F6AZL3' ### name Linie=protein ID
df_63_rfl['Luciferase'] = 'RFL'
df_63_rfl['N-Terminus'] = 'SseI'
df_63_rfl['Normalisierung']='PC12'

df_64_rfl = pd.concat([s140])
df_64_rfl['number'] = np.arange(0,8,1)
df_64_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA064'
df_64_rfl['Konstrukt'] = 'Q9FD10' ### name Linie=protein ID
df_64_rfl['Luciferase'] = 'RFL'
df_64_rfl['N-Terminus'] = 'SseJ'
df_64_rfl['Normalisierung']='PC13'

#####
#Assay 5.10.22
####
file49 = '221005_plate04_nl_001.csv' # nanoluc luciferase plate 4 pJA047, 056, 061
file50 = '221005_plate04_rfl_001.csv' # red firefly luciferase plate 4 pJA47, 056, 061

# Nanluc luciferase plate 4 pJA047, 056, 061

df_p49 = pd.read_csv(file49, skiprows=1)
df_p49 = df_p49.iloc[4:12,:12]
df_p49.columns = ['L1','L2','L3']*4
df_p49 = df_p49.reset_index(drop=True)

s121 = df_p49.iloc[:8,:3]
s131 = df_p49.iloc[:8,4:7]
s133 = df_p49.iloc[:8,8:11]

df_47_nl = pd.concat([s121])
df_47_nl['number'] = np.arange(0,8,1)
df_47_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA047'
df_47_nl['Konstrukt'] = 'Q8ZQC8' ### name Linie=protein ID
df_47_nl['Luciferase'] = 'NL'
df_47_nl['N-Terminus'] = 'SopD2'
df_47_nl['Normalisierung']='PC8'

df_56_nl = pd.concat([s131])
df_56_nl['number'] = np.arange(0,8,1)
df_56_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA056'
df_56_nl['Konstrukt'] = 'A0A0H3NKW1' ### name Linie=protein ID
df_56_nl['Luciferase'] = 'NL'
df_56_nl['N-Terminus'] = 'SsaG'
df_56_nl['Normalisierung']='PC8'

df_61_nl = pd.concat([s133])
df_61_nl['number'] = np.arange(0,8,1)
df_61_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA061'
df_61_nl['Konstrukt'] = 'H9L407' ### name Linie=protein ID
df_61_nl['Luciferase'] = 'NL'
df_61_nl['N-Terminus'] = 'SseF'
df_61_nl['Normalisierung']='PC8'

# RFL plate 4 pJA047, 056, 061

df_p50 = pd.read_csv(file50, skiprows=1)
df_p50 = df_p50.iloc[4:12,:12]
df_p50.columns = ['L1','L2','L3']*4
df_p50 = df_p50.reset_index(drop=True)

s124 = df_p44.iloc[:8,:3]
s134 = df_p44.iloc[:8,4:7]
s136 = df_p44.iloc[:8,8:11]

df_47_rfl = pd.concat([s124])
df_47_rfl['number'] = np.arange(0,8,1)
df_47_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA047'
df_47_rfl['Konstrukt'] = 'Q8ZQC8' ### name Linie=protein ID
df_47_rfl['Luciferase'] = 'RFL'
df_47_rfl['N-Terminus'] = 'SopD2'
df_47_rfl['Normalisierung']='PC8'

df_56_rfl = pd.concat([s134])
df_56_rfl['number'] = np.arange(0,8,1)
df_56_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA056'
df_56_rfl['Konstrukt'] = 'A0A0H3NKW1' ### name Linie=protein ID
df_56_rfl['Luciferase'] = 'RFL'
df_56_rfl['N-Terminus'] = 'SsaG'
df_56_rfl['Normalisierung']='PC8'

df_61_rfl = pd.concat([s136])
df_61_rfl['number'] = np.arange(0,8,1)
df_61_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA061'
df_61_rfl['Konstrukt'] = 'H9L407' ### name Linie=protein ID
df_61_rfl['Luciferase'] = 'RFL'
df_61_rfl['N-Terminus'] = 'SseF'
df_61_rfl['Normalisierung']='PC8'

df_data_rfl = pd.concat([df_pc1_rfl, df_pc2_rfl, df_pc3_rfl, df_pc4_rfl, df_pc5_rfl, df_pc6_rfl, df_pc7_rfl, df_pc8_rfl, df_pc9_rfl, df_pc10_rfl, df_pc11_rfl, df_pc12_rfl, df_pc13_rfl, df_pc_leer_rfl, df_pc_leer2_rfl, df_pc_leer3_rfl, df_pc_leer4_rfl, df_pc_leer5_rfl, df_pc_leer6_rfl, df_nc_rfl, df_nc2_rfl, df_nc3_rfl, df_nc4_rfl, df_nc5_rfl, df_nc6_rfl, df_67_rfl, df_74_rfl, df_76_rfl, df_35_rfl, df_36_rfl, df_37_rfl, df_38_rfl, df_39_rfl, df_31_rfl, df_32_rfl, df_33_rfl, df_34_rfl, df_68_rfl, df_70_rfl, df_72_rfl, df_73_rfl, df_29_rfl, df_40_rfl, df_43_rfl, df_44_rfl, df_48_rfl, df_58_rfl, df_65_rfl, df_66_rfl, df_23_rfl, df_24_rfl, df_25_rfl, df_30_rfl, df_22_rfl, df_26_rfl, df_27_rfl, df_28_rfl, df_41_rfl, df_42_rfl, df_45_rfl, df_46_rfl, df_47_rfl, df_49_rfl, df_50_rfl, df_53_rfl, df_56_rfl, df_57_rfl, df_61_rfl, df_63_rfl, df_64_rfl])
df_data_nl = pd.concat ([df_pc1_nl, df_pc2_nl, df_pc3_nl, df_pc4_nl, df_pc5_nl, df_pc6_nl, df_pc7_nl, df_pc8_nl, df_pc9_nl, df_pc10_nl, df_pc11_nl, df_pc12_nl, df_pc13_nl, df_pc_leer_nl, df_pc_leer2_nl, df_pc_leer3_nl, df_pc_leer4_nl, df_pc_leer5_nl, df_pc_leer6_nl, df_nc_nl, df_nc2_nl, df_nc3_nl, df_nc4_nl, df_nc5_nl, df_nc6_nl, df_67_nl, df_74_nl, df_76_nl, df_35_nl, df_36_nl, df_37_nl, df_38_nl, df_39_nl, df_31_nl, df_32_nl, df_33_nl, df_34_nl, df_68_nl, df_70_nl, df_72_nl, df_73_nl, df_29_nl, df_40_nl, df_43_nl, df_44_nl, df_48_nl, df_58_nl, df_65_nl, df_66_nl, df_23_nl, df_24_nl, df_25_nl, df_30_nl, df_22_nl, df_26_nl, df_27_nl, df_28_nl, df_41_nl, df_42_nl, df_45_nl, df_46_nl, df_47_nl, df_49_nl, df_50_nl, df_53_nl, df_56_nl, df_57_nl, df_61_nl, df_63_nl, df_64_nl])

df_data_rfl.loc[:,['L1','L2','L3']] = df_data_rfl.loc[:,['L1','L2','L3']].astype('float32')
df_data_nl.loc[:,['L1','L2','L3']] = df_data_nl.loc[:,['L1','L2','L3']].astype('float32')

df_data_rfl['L_mean'] = df_data_rfl[['L1','L2','L3']].mean(axis=1)
df_data_nl['L_mean'] = df_data_nl[['L1', 'L2', 'L3']].mean(axis=1)

df_data_rfl['L_std'] = df_data_rfl[['L1','L2','L3']].std(axis=1)
df_data_rfl['L_sem'] = sp.sem(df_data_rfl[['L1','L2','L3']], axis=1, ddof=1)

df_data_nl['L_std'] = df_data_nl[['L1','L2','L3']].std(axis=1)
df_data_nl['L_sem'] = sp.sem(df_data_nl[['L1','L2','L3']], axis=1, ddof=1)

#Verrechnung der Lumineszenzen mit der OD pro Kultur

od_file = 'OD Tabelle alle platten mean, standardabweichung, standarderror of mean.xlsx'

df_od = pd.read_excel(od_file)
df_od.to_csv('CSV_OD_all.csv')

df_od['Linie'] = df_od['Linie'].str.rstrip()
df_od['Linie'] = df_od['Linie'].str.lstrip()

df_od.columns = ['number','Linie','OD1','OD2','OD_mean','OD_std','OD_sem']

df_full_data_rfl= pd.merge(df_data_rfl, df_od, on=['number','Linie'])
df_full_data_nl= pd.merge(df_data_nl, df_od, on=['number','Linie'])

df_long_L_rfl = pd.melt(df_full_data_rfl.loc[:,['L1','L2','L3','Linie','Konstrukt','Luciferase','N-Terminus','Normalisierung','number']], id_vars=['Linie','Konstrukt','Luciferase','N-Terminus','Normalisierung','number'], value_vars=['L1','L2','L3'], value_name='Lumineszenz')

df_long_L_nl = pd.melt(df_full_data_nl.loc[:,['L1','L2','L3','Linie','Konstrukt','Luciferase','N-Terminus','Normalisierung','number']], id_vars=['Linie','Konstrukt','Luciferase','N-Terminus','Normalisierung','number'], value_vars=['L1','L2','L3'], value_name='Lumineszenz')


plt.style.use('tableau-colorblind10')

df_od['OD_mean'] = df_od[['OD1','OD2']].mean(axis=1)
df_full_data_rfl['Lumineszenz_zu_OD'] = df_full_data_rfl['L_mean']/df_full_data_rfl['OD_mean']
df_full_data_nl['Lumineszenz_zu_OD'] = df_full_data_nl['L_mean']/df_full_data_nl['OD_mean']

#Normalisierung der NanoLuc Lumineszenzwerte mithilfe der Positivkontrolle

df_Normalisierungsdaten_PC1_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC1']
df_Normalisierung_PC1daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.1']
df_Normalisierungsdaten_PC1_nl['PC1'] = ['7.351089e+06', '8.581908e+06', '9.061059e+06', '8.851249e+06', '8.411158e+06', '8.937055e+06', '7.473660e+06', '8.735347e+06']*6
df_Normalisierungsdaten_PC1_nl.loc[:,['PC1']] = df_Normalisierungsdaten_PC1_nl.loc[:,['PC1']].astype('float32')
df_Normalisierungsdaten_PC1_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC1_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC1_nl['PC1']

df_Normalisierungsdaten_PC2_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC2']
df_Normalisierung_PC2daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.2']
df_Normalisierungsdaten_PC2_nl['PC2'] = ['7.576740e+06', '6.979040e+06', '7.802230e+06', '6.892805e+06', '6.636355e+06', '6.762475e+06', '7.594412e+06', '7.812635e+06']*6
df_Normalisierungsdaten_PC2_nl.loc[:,['PC2']] = df_Normalisierungsdaten_PC2_nl.loc[:,['PC2']].astype('float32')
df_Normalisierungsdaten_PC2_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC2_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC2_nl['PC2']

df_Normalisierungsdaten_PC3_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC3']
df_Normalisierung_PC3daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.3']
df_Normalisierungsdaten_PC3_nl['PC3'] = ['8.796157e+06', '8.855410e+06', '8.465863e+06', '8.401365e+06', '8.156243e+06', '8.429330e+06', '8.146677e+06', '8.187267e+06']*6
df_Normalisierungsdaten_PC3_nl.loc[:,['PC3']] = df_Normalisierungsdaten_PC3_nl.loc[:,['PC3']].astype('float32')
df_Normalisierungsdaten_PC3_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC3_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC3_nl['PC3']

df_Normalisierungsdaten_PC4_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC4']
df_Normalisierung_PC4daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.4']
df_Normalisierungsdaten_PC4_nl['PC4'] = ['6.850469e+06', '6.333446e+06', '6.179800e+06', '6.958541e+06', '6.744160e+06', '6.530330e+06', '6.544405e+06', '7.017183e+06']*6
df_Normalisierungsdaten_PC4_nl.loc[:,['PC4']] = df_Normalisierungsdaten_PC4_nl.loc[:,['PC4']].astype('float32')
df_Normalisierungsdaten_PC4_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC4_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC4_nl['PC4']

df_Normalisierungsdaten_PC5_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC5']
df_Normalisierung_PC5daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.5']
df_Normalisierungsdaten_PC5_nl['PC5'] = ['5.694271e+06', '6.873909e+06', '6.553678e+06', '6.841687e+06', '6.447520e+06', '7.667787e+06', '8.439253e+06', '7.444012e+06']*6
df_Normalisierungsdaten_PC5_nl.loc[:,['PC5']] = df_Normalisierungsdaten_PC5_nl.loc[:,['PC5']].astype('float32')
df_Normalisierungsdaten_PC5_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC5_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC5_nl['PC5']

df_Normalisierungsdaten_PC6_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC6']
df_Normalisierung_PC6daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.6']
df_Normalisierungsdaten_PC6_nl['PC6'] = ['8.739107e+06', '8.287730e+06', '8.876445e+06', '8.000914e+06', '8.751654e+06', '8.589806e+06', '9.857817e+06', '9.648221e+06']*6
df_Normalisierungsdaten_PC6_nl.loc[:,['PC6']] = df_Normalisierungsdaten_PC6_nl.loc[:,['PC6']].astype('float32')
df_Normalisierungsdaten_PC6_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC6_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC6_nl['PC6']

df_Normalisierungsdaten_PC7_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC7']
df_Normalisierung_PC7daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.7']
df_Normalisierungsdaten_PC7_nl['PC7'] = ['4.417342e+06', '6.074383e+06', '5.608530e+06', '4.481768e+06', '5.306905e+06', '4.767854e+06', '4.678091e+06', '5.664196e+06']*6
df_Normalisierungsdaten_PC7_nl.loc[:,['PC7']] = df_Normalisierungsdaten_PC7_nl.loc[:,['PC7']].astype('float32')
df_Normalisierungsdaten_PC7_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC7_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC7_nl['PC7']

df_Normalisierungsdaten_PC8_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC8']
df_Normalisierung_PC8daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.8']
df_Normalisierungsdaten_PC8_nl['PC8'] = ['3.974966e+06', '4.593724e+06', '4.079018e+06', '4.371810e+06', '4.361777e+06', '4.587055e+06', '3.562189e+06', '4.370518e+06']*6
df_Normalisierungsdaten_PC8_nl.loc[:,['PC8']] = df_Normalisierungsdaten_PC8_nl.loc[:,['PC8']].astype('float32')
df_Normalisierungsdaten_PC8_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC8_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC8_nl['PC8']

df_Normalisierungsdaten_PC9_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC9']
df_Normalisierung_PC9daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.9']
df_Normalisierungsdaten_PC9_nl['PC9'] = ['1.155006e+07', '1.201901e+07', '1.290061e+07', '1.263058e+07', '6.733377e+06', '8.121070e+06', '7.610060e+06', '9.137835e+06']*5
df_Normalisierungsdaten_PC9_nl.loc[:,['PC9']] = df_Normalisierungsdaten_PC9_nl.loc[:,['PC9']].astype('float32')
df_Normalisierungsdaten_PC9_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC9_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC9_nl['PC9']

df_Normalisierungsdaten_PC10_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC10']
df_Normalisierung_PC10daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'Positivkontrolle.10']
df_Normalisierungsdaten_PC10_nl['PC10'] = ['1.008432e+07', '5.537886e+06', '7.395897e+06', '7.136168e+06', '5.504203e+06', '4.412011e+06', '4.531671e+06', '5.058699e+06']*6
df_Normalisierungsdaten_PC10_nl.loc[:,['PC10']] = df_Normalisierungsdaten_PC10_nl.loc[:,['PC10']].astype('float32')
df_Normalisierungsdaten_PC10_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC10_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC10_nl['PC10']

df_Normalisierungsdaten_PC11_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC11']
df_Normalisierung_PC11daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'PC.11']
df_Normalisierungsdaten_PC11_nl['PC11'] = ['1.228867e+07', '1.044546e+07', '1.157442e+07', '5.375442e+06', '6.790207e+06', '6.843869e+06', '8.666467e+06', '7.551865e+06']*5
df_Normalisierungsdaten_PC11_nl.loc[:,['PC11']] = df_Normalisierungsdaten_PC11_nl.loc[:,['PC11']].astype('float32')
df_Normalisierungsdaten_PC11_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC11_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC11_nl['PC11']

df_Normalisierungsdaten_PC12_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC12']
df_Normalisierung_PC12daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'PC.12']
df_Normalisierungsdaten_PC12_nl['PC12'] = ['2.780614e+07', '1.660570e+07', '2.062624e+07', '1.131910e+07', '1.560799e+07', '1.346875e+07', '1.148249e+07', '1.254861e+07']*4
df_Normalisierungsdaten_PC12_nl.loc[:,['PC12']] = df_Normalisierungsdaten_PC12_nl.loc[:,['PC12']].astype('float32')
df_Normalisierungsdaten_PC12_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC12_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC12_nl['PC12']

df_Normalisierungsdaten_PC13_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='PC13']
df_Normalisierung_PC13daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'PC.13']
df_Normalisierungsdaten_PC13_nl['PC13'] = ['2.545186e+07', '3.125881e+07', '3.005353e+07', '4.105982e+07', '2.164075e+07', '1.285664e+07', '1.549629e+07', '2.003955e+07']*2
df_Normalisierungsdaten_PC13_nl.loc[:,['PC13']] = df_Normalisierungsdaten_PC13_nl.loc[:,['PC13']].astype('float32')
df_Normalisierungsdaten_PC13_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC13_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC13_nl['PC13']

df_normalized_data_nl = pd.concat([df_Normalisierungsdaten_PC1_nl, df_Normalisierungsdaten_PC2_nl, df_Normalisierungsdaten_PC3_nl, df_Normalisierungsdaten_PC4_nl, df_Normalisierungsdaten_PC5_nl, df_Normalisierungsdaten_PC6_nl, df_Normalisierungsdaten_PC7_nl, df_Normalisierungsdaten_PC8_nl, df_Normalisierungsdaten_PC9_nl, df_Normalisierungsdaten_PC10_nl, df_Normalisierungsdaten_PC11_nl, df_Normalisierungsdaten_PC12_nl, df_Normalisierungsdaten_PC13_nl])
df_normalized_data_nl = df_normalized_data_nl.sort_values(['Linie']).reset_index(drop=True)

#alle salmonella konstrukte
df_normalized_data_nl_ohne_PC2 =  df_normalized_data_nl.loc[df_normalized_data_nl['Konstrukt'] != 'Positivkontrolle.2']
df_normalized_data_nl_ohne_PC2_PC3 =  df_normalized_data_nl_ohne_PC2.loc[df_normalized_data_nl_ohne_PC2['Konstrukt'] != 'Positivkontrolle.3']
df_normalized_data_nl_ohne_PC2_PC3_PC4 =  df_normalized_data_nl_ohne_PC2_PC3.loc[df_normalized_data_nl_ohne_PC2_PC3['Konstrukt'] != 'Positivkontrolle.4']
df_normalized_data_nl_ohne_PC2345 =  df_normalized_data_nl_ohne_PC2_PC3_PC4.loc[df_normalized_data_nl_ohne_PC2_PC3_PC4['Konstrukt'] != 'Positivkontrolle.5']
df_normalized_data_nl_ohne_PC23456 =  df_normalized_data_nl_ohne_PC2345.loc[df_normalized_data_nl_ohne_PC2345['Konstrukt'] != 'Positivkontrolle.6']
df_normalized_data_nl_ohne_PC234567 =  df_normalized_data_nl_ohne_PC23456.loc[df_normalized_data_nl_ohne_PC23456['Konstrukt'] != 'Positivkontrolle.7']
df_normalized_data_nl_ohne_PC2345678 =  df_normalized_data_nl_ohne_PC234567.loc[df_normalized_data_nl_ohne_PC234567['Konstrukt'] != 'Positivkontrolle.8']
df_normalized_data_nl_ohne_PC23456789  = df_normalized_data_nl_ohne_PC2345678.loc[df_normalized_data_nl_ohne_PC2345678['Konstrukt'] != 'Positivkontrolle.9']
df_normalized_data_nl_ohne_PC2345678910  = df_normalized_data_nl_ohne_PC23456789.loc[df_normalized_data_nl_ohne_PC23456789['Konstrukt'] != 'Positivkontrolle.10']
df_normalized_data_nl_ohne_PC234567891011 = df_normalized_data_nl_ohne_PC2345678910.loc[df_normalized_data_nl_ohne_PC2345678910['Konstrukt'] != 'PC.11']
df_normalized_data_nl_ohne_PC23456789101112 = df_normalized_data_nl_ohne_PC234567891011.loc[df_normalized_data_nl_ohne_PC234567891011['Konstrukt'] != 'PC.12']
df_normalized_data_nl_ohne_PC2345678910111213 = df_normalized_data_nl_ohne_PC23456789101112.loc[df_normalized_data_nl_ohne_PC23456789101112['Konstrukt'] != 'PC.13']
df_normalized_data_nl_ohne_PC_leer1 =  df_normalized_data_nl_ohne_PC2345678910111213.loc[df_normalized_data_nl_ohne_PC2345678910111213['Konstrukt'] != 'Leerkontrolle1']
df_normalized_data_nl_ohne_PC_leer13 =  df_normalized_data_nl_ohne_PC_leer1.loc[df_normalized_data_nl_ohne_PC_leer1['Konstrukt'] != 'Leerkontrolle3']
df_normalized_data_nl_ohne_PC_leer134 =  df_normalized_data_nl_ohne_PC_leer13.loc[df_normalized_data_nl_ohne_PC_leer13['Konstrukt'] != 'Leerkontrolle4']
df_normalized_data_nl_ohne_PC_leer1345 =  df_normalized_data_nl_ohne_PC_leer134.loc[df_normalized_data_nl_ohne_PC_leer134['Konstrukt'] != 'Leerkontrolle5']
df_normalized_data_nl_ohne_PC_leer13456 = df_normalized_data_nl_ohne_PC_leer1345.loc[df_normalized_data_nl_ohne_PC_leer1345['Konstrukt'] != 'Leerkontrolle6']  
df_normalized_data_nl_ohne_PC_leer_nc1 =  df_normalized_data_nl_ohne_PC_leer13456.loc[df_normalized_data_nl_ohne_PC_leer13456['Konstrukt'] != 'Negativkontrolle.1']
df_normalized_data_nl_ohne_PC_leer_nc13 =  df_normalized_data_nl_ohne_PC_leer_nc1.loc[df_normalized_data_nl_ohne_PC_leer_nc1['Konstrukt'] != 'Negativkontrolle.3']
df_normalized_data_nl_ohne_PC_leer_nc134 =  df_normalized_data_nl_ohne_PC_leer_nc13.loc[df_normalized_data_nl_ohne_PC_leer_nc13['Konstrukt'] != 'Negativkontrolle.4']
df_normalized_data_nl_ohne_PC_leer_nc1345 =  df_normalized_data_nl_ohne_PC_leer_nc134.loc[df_normalized_data_nl_ohne_PC_leer_nc134['Konstrukt'] != 'Negativkontrolle.5']
df_normalized_data_nl_ohne_PC_leer_nc13456 = df_normalized_data_nl_ohne_PC_leer_nc1345.loc[df_normalized_data_nl_ohne_PC_leer_nc1345['Konstrukt'] != 'NC.6']

df_normalized_data_nl_ohne_PC_leer_nc13456 = df_normalized_data_nl_ohne_PC_leer_nc13456.sort_values(['Linie']).reset_index(drop=True)

###nur stark sekretierte mit über 60% vom WT in ein Plot 
df_normalized_data_nl_nur_stark_sekretierte_GtgA = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'GtgA']
df_normalized_data_nl_nur_stark_sekretierte_SipC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SipC']
df_normalized_data_nl_nur_stark_sekretierte_SlrP = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SlrP']
df_normalized_data_nl_nur_stark_sekretierte_SopE2 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SopE2']
df_normalized_data_nl_nur_stark_sekretierte_SseL = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseL']
df_normalized_data_nl_nur_stark_sekretierte_SboH = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SboH']
df_normalized_data_nl_nur_stark_sekretierte_SopF = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SopF']
df_normalized_data_nl_nur_stark_sekretierte_SteC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SteC']
df_normalized_data_nl_nur_stark_sekretierte_PC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'PK']
df_normalized_data_nl_nur_stark_sekretierte_NC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'NK']
df_normalized_data_nl_nur_stark_sekretierte_leer = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'LK']

df_normalized_data_nl_stark_sekretierte = pd.concat([df_normalized_data_nl_nur_stark_sekretierte_GtgA, df_normalized_data_nl_nur_stark_sekretierte_SipC, df_normalized_data_nl_nur_stark_sekretierte_SlrP, df_normalized_data_nl_nur_stark_sekretierte_SopE2, df_normalized_data_nl_nur_stark_sekretierte_SseL, df_normalized_data_nl_nur_stark_sekretierte_SboH, df_normalized_data_nl_nur_stark_sekretierte_SopF, df_normalized_data_nl_nur_stark_sekretierte_SteC, df_normalized_data_nl_nur_stark_sekretierte_PC, df_normalized_data_nl_nur_stark_sekretierte_NC, df_normalized_data_nl_nur_stark_sekretierte_leer])
df_normalized_data_nl_stark_sekretierte = df_normalized_data_nl_stark_sekretierte.sort_values(['Linie']).reset_index(drop=True)

### nur gering sekretierte mit unter 20% vom WT im Plot
df_normalized_data_nl_nur_gering_sekretierte_PC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'PK']
df_normalized_data_nl_nur_gering_sekretierte_NC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'NK']
df_normalized_data_nl_nur_gering_sekretierte_leer = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'LK']
df_normalized_data_nl_nur_gering_sekretierte_OrgC = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'OrgC']
df_normalized_data_nl_nur_gering_sekretierte_PipB = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'PipB']
df_normalized_data_nl_nur_gering_sekretierte_PipB2 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'PipB2']
df_normalized_data_nl_nur_gering_sekretierte_PrgI = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'PrgI']
df_normalized_data_nl_nur_gering_sekretierte_SifA = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SifA']
df_normalized_data_nl_nur_gering_sekretierte_SipA = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SipA']
df_normalized_data_nl_nur_gering_sekretierte_SopD2 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SopD2']
df_normalized_data_nl_nur_gering_sekretierte_SpvD = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SpvD']
df_normalized_data_nl_nur_gering_sekretierte_SsaL = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SsaL']
df_normalized_data_nl_nur_gering_sekretierte_SseB = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseB']
df_normalized_data_nl_nur_gering_sekretierte_SseI = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseI']
df_normalized_data_nl_nur_gering_sekretierte_SseJ = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseJ']
df_normalized_data_nl_nur_gering_sekretierte_SseK1 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseK1']
df_normalized_data_nl_nur_gering_sekretierte_SseK2 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseK2']
df_normalized_data_nl_nur_gering_sekretierte_SseK3 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SseK3']
df_normalized_data_nl_nur_gering_sekretierte_SspH2 = df_normalized_data_nl.loc[df_normalized_data_nl['N-Terminus'] == 'SspH2']

df_normalized_data_nl_gering_sekretierte = pd.concat([df_normalized_data_nl_nur_gering_sekretierte_PC, df_normalized_data_nl_nur_gering_sekretierte_NC, df_normalized_data_nl_nur_gering_sekretierte_leer, df_normalized_data_nl_nur_gering_sekretierte_OrgC, df_normalized_data_nl_nur_gering_sekretierte_PipB, df_normalized_data_nl_nur_gering_sekretierte_PipB2, df_normalized_data_nl_nur_gering_sekretierte_PrgI, df_normalized_data_nl_nur_gering_sekretierte_SifA, df_normalized_data_nl_nur_gering_sekretierte_SipA, df_normalized_data_nl_nur_gering_sekretierte_SopD2, df_normalized_data_nl_nur_gering_sekretierte_SpvD, df_normalized_data_nl_nur_gering_sekretierte_SsaL, df_normalized_data_nl_nur_gering_sekretierte_SseB, df_normalized_data_nl_nur_gering_sekretierte_SseI, df_normalized_data_nl_nur_gering_sekretierte_SseJ, df_normalized_data_nl_nur_gering_sekretierte_SseK1, df_normalized_data_nl_nur_gering_sekretierte_SseK2, df_normalized_data_nl_nur_gering_sekretierte_SseK3, df_normalized_data_nl_nur_gering_sekretierte_SspH2])
df_normalized_data_nl_gering_sekretierte = df_normalized_data_nl_gering_sekretierte.sort_values(['Linie']).reset_index(drop=True)

### plot mit mittel sekretierten zwischen 20% und 60% vom WT
df_normalized_data_nl_mittel_sekretierte1 = df_normalized_data_nl_ohne_PC_leer_nc13456.loc[df_normalized_data_nl_ohne_PC_leer_nc13456['N-Terminus'] != 'GtgA']
df_normalized_data_nl_mittel_sekretierte2 = df_normalized_data_nl_mittel_sekretierte1.loc[df_normalized_data_nl_mittel_sekretierte1['N-Terminus'] != 'SipC']
df_normalized_data_nl_mittel_sekretierte3 = df_normalized_data_nl_mittel_sekretierte2.loc[df_normalized_data_nl_mittel_sekretierte2['N-Terminus'] != 'SlrP']
df_normalized_data_nl_mittel_sekretierte4 = df_normalized_data_nl_mittel_sekretierte3.loc[df_normalized_data_nl_mittel_sekretierte3['N-Terminus'] != 'SopE2']
df_normalized_data_nl_mittel_sekretierte5 = df_normalized_data_nl_mittel_sekretierte4.loc[df_normalized_data_nl_mittel_sekretierte4['N-Terminus'] != 'SseL']
df_normalized_data_nl_mittel_sekretierte6 = df_normalized_data_nl_mittel_sekretierte5.loc[df_normalized_data_nl_mittel_sekretierte5['N-Terminus'] != 'SboH']
df_normalized_data_nl_mittel_sekretierte7 = df_normalized_data_nl_mittel_sekretierte6.loc[df_normalized_data_nl_mittel_sekretierte6['N-Terminus'] != 'SopF']
df_normalized_data_nl_mittel_sekretierte8 = df_normalized_data_nl_mittel_sekretierte7.loc[df_normalized_data_nl_mittel_sekretierte7['N-Terminus'] != 'SteC']
df_normalized_data_nl_mittel_sekretierte9 = df_normalized_data_nl_mittel_sekretierte8.loc[df_normalized_data_nl_mittel_sekretierte8['N-Terminus'] != 'OrgC']
df_normalized_data_nl_mittel_sekretierte10 = df_normalized_data_nl_mittel_sekretierte9.loc[df_normalized_data_nl_mittel_sekretierte9['N-Terminus'] != 'PipB']
df_normalized_data_nl_mittel_sekretierte11 = df_normalized_data_nl_mittel_sekretierte10.loc[df_normalized_data_nl_mittel_sekretierte10['N-Terminus'] != 'PipB2']
df_normalized_data_nl_mittel_sekretierte12 = df_normalized_data_nl_mittel_sekretierte11.loc[df_normalized_data_nl_mittel_sekretierte11['N-Terminus'] != 'PrgI']
df_normalized_data_nl_mittel_sekretierte13 = df_normalized_data_nl_mittel_sekretierte12.loc[df_normalized_data_nl_mittel_sekretierte12['N-Terminus'] != 'SifA']
df_normalized_data_nl_mittel_sekretierte14 = df_normalized_data_nl_mittel_sekretierte13.loc[df_normalized_data_nl_mittel_sekretierte13['N-Terminus'] != 'SipA']
df_normalized_data_nl_mittel_sekretierte15 = df_normalized_data_nl_mittel_sekretierte14.loc[df_normalized_data_nl_mittel_sekretierte14['N-Terminus'] != 'SopD2']
df_normalized_data_nl_mittel_sekretierte16 = df_normalized_data_nl_mittel_sekretierte15.loc[df_normalized_data_nl_mittel_sekretierte15['N-Terminus'] != 'SpvD']
df_normalized_data_nl_mittel_sekretierte17 = df_normalized_data_nl_mittel_sekretierte16.loc[df_normalized_data_nl_mittel_sekretierte16['N-Terminus'] != 'SsaL']
df_normalized_data_nl_mittel_sekretierte18 = df_normalized_data_nl_mittel_sekretierte17.loc[df_normalized_data_nl_mittel_sekretierte17['N-Terminus'] != 'SseB']
df_normalized_data_nl_mittel_sekretierte19 = df_normalized_data_nl_mittel_sekretierte18.loc[df_normalized_data_nl_mittel_sekretierte18['N-Terminus'] != 'SseI']
df_normalized_data_nl_mittel_sekretierte20 = df_normalized_data_nl_mittel_sekretierte19.loc[df_normalized_data_nl_mittel_sekretierte19['N-Terminus'] != 'SseJ']
df_normalized_data_nl_mittel_sekretierte21 = df_normalized_data_nl_mittel_sekretierte20.loc[df_normalized_data_nl_mittel_sekretierte20['N-Terminus'] != 'SseK1']
df_normalized_data_nl_mittel_sekretierte22 = df_normalized_data_nl_mittel_sekretierte21.loc[df_normalized_data_nl_mittel_sekretierte21['N-Terminus'] != 'SseK2']
df_normalized_data_nl_mittel_sekretierte23 = df_normalized_data_nl_mittel_sekretierte22.loc[df_normalized_data_nl_mittel_sekretierte22['N-Terminus'] != 'SseK3']
df_normalized_data_nl_mittel_sekretierte24 = df_normalized_data_nl_mittel_sekretierte23.loc[df_normalized_data_nl_mittel_sekretierte23['N-Terminus'] != 'SspH2']

df_normalized_data_nl_mittel_sekretierte24 = df_normalized_data_nl_mittel_sekretierte24.sort_values(['Linie']).reset_index(drop=True)

colors = ["b"]
sns.set(font_scale=1.5)
plt.figure(figsize=(15, 10))
graph=sns.barplot(data=df_normalized_data_nl_ohne_PC_leer_nc13456, x='N-Terminus', y='Lumineszenz/OD_normalisiert', palette=colors, ci = 'sd', dodge=False)
graph.axhline(1)
plt.title('Mittlere Lumineszenz/OD normalisiert NanoLuc Luciferase (+Standardabweichung)', y=1.02)
plt.xticks(rotation=45)
plt.show()
plt.savefig('Nanoluc_alle_salmonella_konstrukte.svg', dpi=600)

#Normalisierung der Red Firefly Lumineszenzwerte mithilfe der Positivkontrolle

df_Normalisierungsdaten_PC1_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC1']
df_Normalisierung_PC1daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.1']
df_Normalisierungsdaten_PC1_rfl['PC1'] = ['1195.498361', '1502.239122', '2063.645310', '1342.507904', '1358.549856', '1500.846369', '1145.264921', '1020.328267']*6
df_Normalisierungsdaten_PC1_rfl.loc[:,['PC1']] = df_Normalisierungsdaten_PC1_rfl.loc[:,['PC1']].astype('float32')
df_Normalisierungsdaten_PC1_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC1_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC1_rfl['PC1']

df_Normalisierungsdaten_PC2_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC2']
df_Normalisierung_PC2daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.2']
df_Normalisierungsdaten_PC2_rfl['PC2'] = ['1972.260851', '1700.842840', '2000.148472', '1930.648329', '1849.184874', '1658.594132', '2025.588030', '1580.623483']*6
df_Normalisierungsdaten_PC2_rfl.loc[:,['PC2']] = df_Normalisierungsdaten_PC2_rfl.loc[:,['PC2']].astype('float32')
df_Normalisierungsdaten_PC2_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC2_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC2_rfl['PC2']

df_Normalisierungsdaten_PC3_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC3']
df_Normalisierung_PC3daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.3']
df_Normalisierungsdaten_PC3_rfl['PC3'] = ['1464.888835', '1745.777832', '1546.666667', '1873.777832', '1933.333282', '1628.444499', '2043.333282', '2183.111165']*6
df_Normalisierungsdaten_PC3_rfl.loc[:,['PC3']] = df_Normalisierungsdaten_PC3_rfl.loc[:,['PC3']].astype('float32')
df_Normalisierungsdaten_PC3_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC3_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC3_rfl['PC3']

df_Normalisierungsdaten_PC4_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC4']
df_Normalisierung_PC4daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.4']
df_Normalisierungsdaten_PC4_rfl['PC4'] = ['1025.882353', '1032.982435', '792.000000', '956.862721', '1160.784338', '1250.000000', '1223.529412', '1326.666718']*6
df_Normalisierungsdaten_PC4_rfl.loc[:,['PC4']] = df_Normalisierungsdaten_PC4_rfl.loc[:,['PC4']].astype('float32')
df_Normalisierungsdaten_PC4_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC4_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC4_rfl['PC4']

df_Normalisierungsdaten_PC5_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC5']
df_Normalisierung_PC5daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.5']
df_Normalisierungsdaten_PC5_rfl['PC5'] = ['1358.431325', '1511.111165', '1295.686322', '1406.666718', '1449.411765', '1671.111165', '1752.888835', '1179.607867']*6
df_Normalisierungsdaten_PC5_rfl.loc[:,['PC5']] = df_Normalisierungsdaten_PC5_rfl.loc[:,['PC5']].astype('float32')
df_Normalisierungsdaten_PC5_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC5_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC5_rfl['PC5']

df_Normalisierungsdaten_PC6_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC6']
df_Normalisierung_PC6daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.6']
df_Normalisierungsdaten_PC6_rfl['PC6'] = ['3011.555664', '2599.111165', '2869.333333', '2876.444336', '2925.714286', '2571.428571', '3127.618931', '3015.111003']*6
df_Normalisierungsdaten_PC6_rfl.loc[:,['PC6']] = df_Normalisierungsdaten_PC6_rfl.loc[:,['PC6']].astype('float32')
df_Normalisierungsdaten_PC6_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC6_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC6_rfl['PC6']

df_Normalisierungsdaten_PC7_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC7']
df_Normalisierung_PC7daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.7']
df_Normalisierungsdaten_PC7_rfl['PC7'] = ['2768.000000', '4033.333435', '3814.902057', '2890.158808', '4094.117647', '2986.666581', '3073.684211', '2800.000000']*6
df_Normalisierungsdaten_PC7_rfl.loc[:,['PC7']] = df_Normalisierungsdaten_PC7_rfl.loc[:,['PC7']].astype('float32')
df_Normalisierungsdaten_PC7_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC7_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC7_rfl['PC7']

df_Normalisierungsdaten_PC8_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC8']
df_Normalisierung_PC8daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.8']
df_Normalisierungsdaten_PC8_rfl['PC8'] = ['1493.333333', '1768.888889', '1953.684211', '2311.111111', '2298.947368', '2082.962918', '1871.515225', '1786.666626']*6
df_Normalisierungsdaten_PC8_rfl.loc[:,['PC8']] = df_Normalisierungsdaten_PC8_rfl.loc[:,['PC8']].astype('float32')
df_Normalisierungsdaten_PC8_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC8_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC8_rfl['PC8']

df_Normalisierungsdaten_PC9_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC9']
df_Normalisierung_PC9daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.9']
df_Normalisierungsdaten_PC9_rfl['PC9'] = ['3192.452120', '6930.773276', '6742.875895', '4557.178985', '1730.307701', '3013.872923', '929.832984', '966.375493']*5
df_Normalisierungsdaten_PC9_rfl.loc[:,['PC9']] = df_Normalisierungsdaten_PC9_rfl.loc[:,['PC9']].astype('float32')
df_Normalisierungsdaten_PC9_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC9_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC9_rfl['PC9']

df_Normalisierungsdaten_PC10_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC10']
df_Normalisierung_PC10daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'Positivkontrolle.10']
df_Normalisierungsdaten_PC10_rfl['PC10'] = ['2623.626083', '1439.730478', '1115.770045', '1659.160779', '1028.704530', '1113.268463', '696.784121', '814.825517']*6
df_Normalisierungsdaten_PC10_rfl.loc[:,['PC10']] = df_Normalisierungsdaten_PC10_rfl.loc[:,['PC10']].astype('float32')
df_Normalisierungsdaten_PC10_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC10_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC10_rfl['PC10']

df_Normalisierungsdaten_PC11_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC11']
df_Normalisierung_PC11daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'PC.11']
df_Normalisierungsdaten_PC11_rfl['PC11'] = ['1910.617201', '1763.148712', '2753.665672', '1601.081483', '1502.649074', '1524.962822', '2401.269815', '1549.309748']*5
df_Normalisierungsdaten_PC11_rfl.loc[:,['PC11']] = df_Normalisierungsdaten_PC11_rfl.loc[:,['PC11']].astype('float32')
df_Normalisierungsdaten_PC11_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC11_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC11_rfl['PC11']

df_Normalisierungsdaten_PC12_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC12']
df_Normalisierung_PC12daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'PC.12']
df_Normalisierungsdaten_PC12_rfl['PC12'] = ['2841.176747', '1815.063210', '2419.383545', '1352.705148', '2223.078658', '1823.439239', '1723.865243', '1583.932899']*4
df_Normalisierungsdaten_PC12_rfl.loc[:,['PC12']] = df_Normalisierungsdaten_PC12_rfl.loc[:,['PC12']].astype('float32')
df_Normalisierungsdaten_PC12_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC12_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC12_rfl['PC12']

df_Normalisierungsdaten_PC13_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='PC13']
df_Normalisierung_PC13daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'PC.13']
df_Normalisierungsdaten_PC13_rfl['PC13'] = ['2834.557338', '3532.044597', '3758.830472', '4715.701289', '3305.793045', '1706.495383', '1818.860078', '2618.204689']*2
df_Normalisierungsdaten_PC13_rfl.loc[:,['PC13']] = df_Normalisierungsdaten_PC13_rfl.loc[:,['PC13']].astype('float32')
df_Normalisierungsdaten_PC13_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_PC13_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_PC13_rfl['PC13']

df_normalized_data_rfl = pd.concat([df_Normalisierungsdaten_PC1_rfl, df_Normalisierungsdaten_PC2_rfl, df_Normalisierungsdaten_PC3_rfl, df_Normalisierungsdaten_PC4_rfl, df_Normalisierungsdaten_PC5_rfl, df_Normalisierungsdaten_PC6_rfl, df_Normalisierungsdaten_PC7_rfl, df_Normalisierungsdaten_PC8_rfl, df_Normalisierungsdaten_PC9_rfl, df_Normalisierungsdaten_PC10_rfl, df_Normalisierungsdaten_PC11_rfl, df_Normalisierungsdaten_PC12_rfl, df_Normalisierungsdaten_PC13_rfl])
df_normalized_data_rfl = df_normalized_data_rfl.sort_values(['Linie']).reset_index(drop=True)

#alle salmonella konstrukte
df_normalized_data_rfl_ohne_PC2 =  df_normalized_data_rfl.loc[df_normalized_data_rfl['Konstrukt'] != 'Positivkontrolle.2']
df_normalized_data_rfl_ohne_PC2_PC3 =  df_normalized_data_rfl_ohne_PC2.loc[df_normalized_data_rfl_ohne_PC2['Konstrukt'] != 'Positivkontrolle.3']
df_normalized_data_rfl_ohne_PC2_PC3_PC4 =  df_normalized_data_rfl_ohne_PC2_PC3.loc[df_normalized_data_rfl_ohne_PC2_PC3['Konstrukt'] != 'Positivkontrolle.4']
df_normalized_data_rfl_ohne_PC2345 =  df_normalized_data_rfl_ohne_PC2_PC3_PC4.loc[df_normalized_data_rfl_ohne_PC2_PC3_PC4['Konstrukt'] != 'Positivkontrolle.5']
df_normalized_data_rfl_ohne_PC23456 =  df_normalized_data_rfl_ohne_PC2345.loc[df_normalized_data_rfl_ohne_PC2345['Konstrukt'] != 'Positivkontrolle.6']
df_normalized_data_rfl_ohne_PC234567 =  df_normalized_data_rfl_ohne_PC23456.loc[df_normalized_data_rfl_ohne_PC23456['Konstrukt'] != 'Positivkontrolle.7']
df_normalized_data_rfl_ohne_PC2345678 =  df_normalized_data_rfl_ohne_PC234567.loc[df_normalized_data_rfl_ohne_PC234567['Konstrukt'] != 'Positivkontrolle.8']
df_normalized_data_rfl_ohne_PC23456789 =  df_normalized_data_rfl_ohne_PC2345678.loc[df_normalized_data_rfl_ohne_PC2345678['Konstrukt'] != 'Positivkontrolle.9']
df_normalized_data_rfl_ohne_PC2345678910 =  df_normalized_data_rfl_ohne_PC23456789.loc[df_normalized_data_rfl_ohne_PC23456789['Konstrukt'] != 'Positivkontrolle.10']
df_normalized_data_rfl_ohne_PC234567891011 = df_normalized_data_rfl_ohne_PC2345678910.loc[df_normalized_data_rfl_ohne_PC2345678910['Konstrukt'] != 'PC.11']
df_normalized_data_rfl_ohne_PC23456789101112 = df_normalized_data_rfl_ohne_PC234567891011.loc[df_normalized_data_rfl_ohne_PC234567891011['Konstrukt'] != 'PC.12']
df_normalized_data_rfl_ohne_PC2345678910111213 = df_normalized_data_rfl_ohne_PC23456789101112.loc[df_normalized_data_rfl_ohne_PC23456789101112['Konstrukt'] != 'PC.13']
df_normalized_data_rfl_ohne_PC_leer1 =  df_normalized_data_rfl_ohne_PC2345678910111213.loc[df_normalized_data_rfl_ohne_PC2345678910111213['Konstrukt'] != 'Leerkontrolle1']
df_normalized_data_rfl_ohne_PC_leer13 =  df_normalized_data_rfl_ohne_PC_leer1.loc[df_normalized_data_rfl_ohne_PC_leer1['Konstrukt'] != 'Leerkontrolle3']
df_normalized_data_rfl_ohne_PC_leer134 =  df_normalized_data_rfl_ohne_PC_leer13.loc[df_normalized_data_rfl_ohne_PC_leer13['Konstrukt'] != 'Leerkontrolle4']
df_normalized_data_rfl_ohne_PC_leer1345 =  df_normalized_data_rfl_ohne_PC_leer134.loc[df_normalized_data_rfl_ohne_PC_leer134['Konstrukt'] != 'Leerkontrolle5']
df_normalized_data_rfl_ohne_PC_leer13456 = df_normalized_data_rfl_ohne_PC_leer1345.loc[df_normalized_data_rfl_ohne_PC_leer1345['Konstrukt'] != 'Leerkontrolle6']
df_normalized_data_rfl_ohne_PC_leer_nc1 =  df_normalized_data_rfl_ohne_PC_leer13456.loc[df_normalized_data_rfl_ohne_PC_leer13456['Konstrukt'] != 'Negativkontrolle.1']
df_normalized_data_rfl_ohne_PC_leer_nc13 =  df_normalized_data_rfl_ohne_PC_leer_nc1.loc[df_normalized_data_rfl_ohne_PC_leer_nc1['Konstrukt'] != 'Negativkontrolle.3']
df_normalized_data_rfl_ohne_PC_leer_nc134 =  df_normalized_data_rfl_ohne_PC_leer_nc13.loc[df_normalized_data_rfl_ohne_PC_leer_nc13['Konstrukt'] != 'Negativkontrolle.4']
df_normalized_data_rfl_ohne_PC_leer_nc1345 =  df_normalized_data_rfl_ohne_PC_leer_nc134.loc[df_normalized_data_rfl_ohne_PC_leer_nc134['Konstrukt'] != 'Negativkontrolle.5']
df_normalized_data_rfl_ohne_PC_leer_nc13456 = df_normalized_data_rfl_ohne_PC_leer_nc1345.loc[df_normalized_data_rfl_ohne_PC_leer_nc1345['Konstrukt'] != 'NC.6']

df_normalized_data_rfl_ohne_PC_leer_nc13456 = df_normalized_data_rfl_ohne_PC_leer_nc13456.sort_values(['Linie']).reset_index(drop=True)

###nur stark sekretierte mit über 60% vom WT in ein Plot 
df_normalized_data_rfl_nur_stark_sekretierte_GtgA = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'GtgA']
df_normalized_data_rfl_nur_stark_sekretierte_SipC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SipC']
df_normalized_data_rfl_nur_stark_sekretierte_SlrP = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SlrP']
df_normalized_data_rfl_nur_stark_sekretierte_SopE2 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SopE2']
df_normalized_data_rfl_nur_stark_sekretierte_SseL = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseL']
df_normalized_data_rfl_nur_stark_sekretierte_SboH = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SboH']
df_normalized_data_rfl_nur_stark_sekretierte_SopF = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SopF']
df_normalized_data_rfl_nur_stark_sekretierte_SteC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SteC']
df_normalized_data_rfl_nur_stark_sekretierte_PC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'PK']
df_normalized_data_rfl_nur_stark_sekretierte_NC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'NK']
df_normalized_data_rfl_nur_stark_sekretierte_leer = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'LK']

df_normalized_data_rfl_stark_sekretierte = pd.concat([df_normalized_data_rfl_nur_stark_sekretierte_GtgA, df_normalized_data_rfl_nur_stark_sekretierte_SipC, df_normalized_data_rfl_nur_stark_sekretierte_SlrP, df_normalized_data_rfl_nur_stark_sekretierte_SopE2, df_normalized_data_rfl_nur_stark_sekretierte_SseL, df_normalized_data_rfl_nur_stark_sekretierte_SboH, df_normalized_data_rfl_nur_stark_sekretierte_SopF, df_normalized_data_rfl_nur_stark_sekretierte_SteC, df_normalized_data_rfl_nur_stark_sekretierte_PC, df_normalized_data_rfl_nur_stark_sekretierte_NC, df_normalized_data_rfl_nur_stark_sekretierte_leer])
df_normalized_data_rfl_stark_sekretierte = df_normalized_data_rfl_stark_sekretierte.sort_values(['Linie']).reset_index(drop=True)

### nur gering sekretierte mit unter 20% vom WT im Plot
df_normalized_data_rfl_nur_gering_sekretierte_PC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'PK']
df_normalized_data_rfl_nur_gering_sekretierte_NC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'NK']
df_normalized_data_rfl_nur_gering_sekretierte_leer = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'LK']
df_normalized_data_rfl_nur_gering_sekretierte_OrgC = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'OrgC']
df_normalized_data_rfl_nur_gering_sekretierte_PipB = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'PipB']
df_normalized_data_rfl_nur_gering_sekretierte_PipB2 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'PipB2']
df_normalized_data_rfl_nur_gering_sekretierte_PrgI = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'PrgI']
df_normalized_data_rfl_nur_gering_sekretierte_SifA = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SifA']
df_normalized_data_rfl_nur_gering_sekretierte_SipA = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SipA']
df_normalized_data_rfl_nur_gering_sekretierte_SopD2 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SopD2']
df_normalized_data_rfl_nur_gering_sekretierte_SpvD = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SpvD']
df_normalized_data_rfl_nur_gering_sekretierte_SsaL = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SsaL']
df_normalized_data_rfl_nur_gering_sekretierte_SseB = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseB']
df_normalized_data_rfl_nur_gering_sekretierte_SseI = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseI']
df_normalized_data_rfl_nur_gering_sekretierte_SseJ = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseJ']
df_normalized_data_rfl_nur_gering_sekretierte_SseK1 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseK1']
df_normalized_data_rfl_nur_gering_sekretierte_SseK2 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseK2']
df_normalized_data_rfl_nur_gering_sekretierte_SseK3 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SseK3']
df_normalized_data_rfl_nur_gering_sekretierte_SspH2 = df_normalized_data_rfl.loc[df_normalized_data_rfl['N-Terminus'] == 'SspH2']

df_normalized_data_rfl_gering_sekretierte = pd.concat([df_normalized_data_rfl_nur_gering_sekretierte_PC, df_normalized_data_rfl_nur_gering_sekretierte_NC, df_normalized_data_rfl_nur_gering_sekretierte_leer, df_normalized_data_rfl_nur_gering_sekretierte_OrgC, df_normalized_data_rfl_nur_gering_sekretierte_PipB, df_normalized_data_rfl_nur_gering_sekretierte_PipB2, df_normalized_data_rfl_nur_gering_sekretierte_PrgI, df_normalized_data_rfl_nur_gering_sekretierte_SifA, df_normalized_data_rfl_nur_gering_sekretierte_SipA, df_normalized_data_rfl_nur_gering_sekretierte_SopD2, df_normalized_data_rfl_nur_gering_sekretierte_SpvD, df_normalized_data_rfl_nur_gering_sekretierte_SsaL, df_normalized_data_rfl_nur_gering_sekretierte_SseB, df_normalized_data_rfl_nur_gering_sekretierte_SseI, df_normalized_data_rfl_nur_gering_sekretierte_SseJ, df_normalized_data_rfl_nur_gering_sekretierte_SseK1, df_normalized_data_rfl_nur_gering_sekretierte_SseK2, df_normalized_data_rfl_nur_gering_sekretierte_SseK3, df_normalized_data_rfl_nur_gering_sekretierte_SspH2])
df_normalized_data_rfl_gering_sekretierte = df_normalized_data_rfl_gering_sekretierte.sort_values(['Linie']).reset_index(drop=True)

### plot mit mittel sekretierten zwischen 20% und 60% vom WT
df_normalized_data_rfl_mittel_sekretierte1 = df_normalized_data_rfl_ohne_PC_leer_nc13456.loc[df_normalized_data_rfl_ohne_PC_leer_nc13456['N-Terminus'] != 'GtgA']
df_normalized_data_rfl_mittel_sekretierte2 = df_normalized_data_rfl_mittel_sekretierte1.loc[df_normalized_data_rfl_mittel_sekretierte1['N-Terminus'] != 'SipC']
df_normalized_data_rfl_mittel_sekretierte3 = df_normalized_data_rfl_mittel_sekretierte2.loc[df_normalized_data_rfl_mittel_sekretierte2['N-Terminus'] != 'SlrP']
df_normalized_data_rfl_mittel_sekretierte4 = df_normalized_data_rfl_mittel_sekretierte3.loc[df_normalized_data_rfl_mittel_sekretierte3['N-Terminus'] != 'SopE2']
df_normalized_data_rfl_mittel_sekretierte5 = df_normalized_data_rfl_mittel_sekretierte4.loc[df_normalized_data_rfl_mittel_sekretierte4['N-Terminus'] != 'SseL']
df_normalized_data_rfl_mittel_sekretierte6 = df_normalized_data_rfl_mittel_sekretierte5.loc[df_normalized_data_rfl_mittel_sekretierte5['N-Terminus'] != 'SboH']
df_normalized_data_rfl_mittel_sekretierte7 = df_normalized_data_rfl_mittel_sekretierte6.loc[df_normalized_data_rfl_mittel_sekretierte6['N-Terminus'] != 'SopF']
df_normalized_data_rfl_mittel_sekretierte8 = df_normalized_data_rfl_mittel_sekretierte7.loc[df_normalized_data_rfl_mittel_sekretierte7['N-Terminus'] != 'SteC']
df_normalized_data_rfl_mittel_sekretierte9 = df_normalized_data_rfl_mittel_sekretierte8.loc[df_normalized_data_rfl_mittel_sekretierte8['N-Terminus'] != 'OrgC']
df_normalized_data_rfl_mittel_sekretierte10 = df_normalized_data_rfl_mittel_sekretierte9.loc[df_normalized_data_rfl_mittel_sekretierte9['N-Terminus'] != 'PipB']
df_normalized_data_rfl_mittel_sekretierte11 = df_normalized_data_rfl_mittel_sekretierte10.loc[df_normalized_data_rfl_mittel_sekretierte10['N-Terminus'] != 'PipB2']
df_normalized_data_rfl_mittel_sekretierte12 = df_normalized_data_rfl_mittel_sekretierte11.loc[df_normalized_data_rfl_mittel_sekretierte11['N-Terminus'] != 'PrgI']
df_normalized_data_rfl_mittel_sekretierte13 = df_normalized_data_rfl_mittel_sekretierte12.loc[df_normalized_data_rfl_mittel_sekretierte12['N-Terminus'] != 'SifA']
df_normalized_data_rfl_mittel_sekretierte14 = df_normalized_data_rfl_mittel_sekretierte13.loc[df_normalized_data_rfl_mittel_sekretierte13['N-Terminus'] != 'SipA']
df_normalized_data_rfl_mittel_sekretierte15 = df_normalized_data_rfl_mittel_sekretierte14.loc[df_normalized_data_rfl_mittel_sekretierte14['N-Terminus'] != 'SopD2']
df_normalized_data_rfl_mittel_sekretierte16 = df_normalized_data_rfl_mittel_sekretierte15.loc[df_normalized_data_rfl_mittel_sekretierte15['N-Terminus'] != 'SpvD']
df_normalized_data_rfl_mittel_sekretierte17 = df_normalized_data_rfl_mittel_sekretierte16.loc[df_normalized_data_rfl_mittel_sekretierte16['N-Terminus'] != 'SsaL']
df_normalized_data_rfl_mittel_sekretierte18 = df_normalized_data_rfl_mittel_sekretierte17.loc[df_normalized_data_rfl_mittel_sekretierte17['N-Terminus'] != 'SseB']
df_normalized_data_rfl_mittel_sekretierte19 = df_normalized_data_rfl_mittel_sekretierte18.loc[df_normalized_data_rfl_mittel_sekretierte18['N-Terminus'] != 'SseI']
df_normalized_data_rfl_mittel_sekretierte20 = df_normalized_data_rfl_mittel_sekretierte19.loc[df_normalized_data_rfl_mittel_sekretierte19['N-Terminus'] != 'SseJ']
df_normalized_data_rfl_mittel_sekretierte21 = df_normalized_data_rfl_mittel_sekretierte20.loc[df_normalized_data_rfl_mittel_sekretierte20['N-Terminus'] != 'SseK1']
df_normalized_data_rfl_mittel_sekretierte22 = df_normalized_data_rfl_mittel_sekretierte21.loc[df_normalized_data_rfl_mittel_sekretierte21['N-Terminus'] != 'SseK2']
df_normalized_data_rfl_mittel_sekretierte23 = df_normalized_data_rfl_mittel_sekretierte22.loc[df_normalized_data_rfl_mittel_sekretierte22['N-Terminus'] != 'SseK3']
df_normalized_data_rfl_mittel_sekretierte24 = df_normalized_data_rfl_mittel_sekretierte23.loc[df_normalized_data_rfl_mittel_sekretierte23['N-Terminus'] != 'SspH2']

df_normalized_data_rfl_mittel_sekretierte24 = df_normalized_data_rfl_mittel_sekretierte24.sort_values(['Linie']).reset_index(drop=True)

colors = ["b"]
sns.set(font_scale=1.5)
plt.figure(figsize=(15, 10))
graph2=sns.barplot(data=df_normalized_data_rfl_ohne_PC_leer_nc13456, x='N-Terminus', y='Lumineszenz/OD_normalisiert', palette=colors, ci = 'sd', dodge=False)
graph2.axhline(1)
plt.title('Mittlere Lumineszenz/OD normalisiert Red Firefly Luciferase (+Standardabweichung)', y=1.02)
plt.xticks(rotation=45)
plt.show()
plt.savefig('RFL_alle_salmonella_konstrukte.svg', dpi=600)



