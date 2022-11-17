import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import stats
import scipy.stats as sp 

file1 = '220915_plate_01_nl_001.csv' # nanoluc luciferase plate 1 Positiv Kontrollen + pMP053
file2 = '220915_plate_01_rfl_001.csv' # red firefly luciferase plate 1 Positiv Kontrollen + pMP053
file3 = '220915_plate_02_nl_001.csv' # nanoluc luciferase plate 2 negativ kontrollen + pMP057
file4 = '220915_plate_02_rfl_001.csv' # red firefly luciferase plate 2 negativ kontrollen + pMP057
file5 = '220915_plate_03_nl_001.csv' # nanoluc luciferase plate 3 pJA003, 005, 008
file6 = '220915_plate_03_rfl_001.csv' # red firefly luciferase plate 3 pJA003, 005, 008
file7 = '220915_plate_04_nl_001.csv' # nanoluc luciferase plate 4 pJA009, 017, 019
file8 = '220915_plate_04_rfl_001.csv' # red firefly luciferase plate 4 pJA009, 017, 019

# nanoluc luciferase plate 1 positiv Kontrollen + pMP053

df_p1 = pd.read_csv(file1, skiprows=1)
df_p1 = df_p1.iloc[4:12,:12]
df_p1.columns = ['L1','L2','L3']*4
df_p1 = df_p1.reset_index(drop=True)

s1 = df_p1.iloc[:8,:3]
s2 = df_p1.iloc[:8,4:7]
s3 = df_p1.iloc[:8,8:11]

df_SptP1_nl = pd.concat([s1])
df_SptP1_nl['time'] = 4
df_SptP1_nl['number'] = np.arange(0,8,1)
df_SptP1_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.1'
df_SptP1_nl['Konstrukt'] = 'PK2'
df_SptP1_nl['Luciferase'] = 'NL'
df_SptP1_nl['Normalisierung']='SptP1'

df_SptP2_nl = pd.concat([s2])
df_SptP2_nl['time'] = 4
df_SptP2_nl['number'] = np.arange(0,8,1)
df_SptP2_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.2'
df_SptP2_nl['Konstrukt'] = 'PK'
df_SptP2_nl['Luciferase'] = 'NL'
df_SptP2_nl['Normalisierung']='SptP2'

df_53_nl = pd.concat([s3])
df_53_nl['time'] = 4
df_53_nl['number'] = np.arange(0,8,1)
df_53_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP053'
df_53_nl['Konstrukt'] = 'YopH' 
df_53_nl['Luciferase'] = 'NL'
df_53_nl['Normalisierung']='SptP1'

# red firefly luciferase plate 1 positiv Kontrollen + pMP053

df_p2 = pd.read_csv(file2, skiprows=1)
df_p2 = df_p2.iloc[4:12,:12]
df_p2.columns = ['L1','L2','L3']*4
df_p2 = df_p2.reset_index(drop=True)

s4 = df_p2.iloc[:8,:3]
s5 = df_p2.iloc[:8,4:7]
s6 = df_p2.iloc[:8,8:11]
  
df_SptP1_rfl = pd.concat([s4])
df_SptP1_rfl['time'] = 4
df_SptP1_rfl['number'] = np.arange(0,8,1)
df_SptP1_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.1'
df_SptP1_rfl['Konstrukt'] = 'PK2'
df_SptP1_rfl['Luciferase'] = 'RFL'
df_SptP1_rfl['Normalisierung']='SptP1'

df_SptP2_rfl = pd.concat([s5])
df_SptP2_rfl['time'] = 4
df_SptP2_rfl['number'] = np.arange(0,8,1)
df_SptP2_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP049.2'
df_SptP2_rfl['Konstrukt'] = 'PK'
df_SptP2_rfl['Luciferase'] = 'RFL'
df_SptP2_rfl['Normalisierung']='SptP2'

df_53_rfl = pd.concat([s6])
df_53_rfl['time'] = 4
df_53_rfl['number'] = np.arange(0,8,1)
df_53_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP053'
df_53_rfl['Konstrukt'] = 'YopH' 
df_53_rfl['Luciferase'] = 'RFL'
df_53_rfl['Normalisierung']='SptP1'

# nanoluc luciferase plate 2 negativ Kontrollen + pMP057

df_p3 = pd.read_csv(file3, skiprows=1)
df_p3 = df_p3.iloc[4:12,:12]
df_p3.columns = ['L1','L2','L3']*4
df_p3 = df_p3.reset_index(drop=True)

s7 = df_p3.iloc[:8,:3]
s8 = df_p3.iloc[:8,4:7]
s9 = df_p3.iloc[:8,8:11]

df_SptP_leer_nl = pd.concat([s7])
df_SptP_leer_nl['time'] = 4
df_SptP_leer_nl['number'] = np.arange(0,8,1)
df_SptP_leer_nl['Linie'] = 'SB905 dSipA dSptP pHilA'
df_SptP_leer_nl['Konstrukt'] = 'LK'
df_SptP_leer_nl['Luciferase'] = 'NL'
df_SptP_leer_nl['Normalisierung']='SptP1'

df_NK_nl = pd.concat([s8])
df_NK_nl['time'] = 4
df_NK_nl['number'] = np.arange(0,8,1)
df_NK_nl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049'
df_NK_nl['Konstrukt'] = 'NK'
df_NK_nl['Luciferase'] = 'NL'
df_NK_nl['Normalisierung']='SptP1'

df_57_nl = pd.concat([s9])
df_57_nl['time'] = 4
df_57_nl['number'] = np.arange(0,8,1)
df_57_nl['Linie'] = 'SB905 dSipA dSptP pHilA pMP057'
df_57_nl['Konstrukt'] = 'YopO' 
df_57_nl['Luciferase'] = 'NL'
df_57_nl['Normalisierung']='SptP1'

# red firefly luciferase plate 2 negativ Kontrollen + pMP057

df_p4 = pd.read_csv(file4, skiprows=1)
df_p4 = df_p4.iloc[4:12,:12]
df_p4.columns = ['L1','L2','L3']*4
df_p4 = df_p4.reset_index(drop=True)

s10 = df_p4.iloc[:8,:3]
s11 = df_p4.iloc[:8,4:7]
s12 = df_p4.iloc[:8,8:11]
  
df_SptP_leer_rfl = pd.concat([s10])
df_SptP_leer_rfl['time'] = 4
df_SptP_leer_rfl['number'] = np.arange(0,8,1)
df_SptP_leer_rfl['Linie'] = 'SB905 dSipA dSptP pHilA'
df_SptP_leer_rfl['Konstrukt'] = 'LK'
df_SptP_leer_rfl['Luciferase'] = 'RFL'
df_SptP_leer_rfl['Normalisierung']='SptP1'

df_NK_rfl = pd.concat([s11])
df_NK_rfl['time'] = 4
df_NK_rfl['number'] = np.arange(0,8,1)
df_NK_rfl['Linie'] = 'SB905 dSipA dSptP dInvA pHilA pMP049'
df_NK_rfl['Konstrukt'] = 'NK'
df_NK_rfl['Luciferase'] = 'RFL'
df_NK_rfl['Normalisierung']='SptP1'

df_57_rfl = pd.concat([s12])
df_57_rfl['time'] = 4
df_57_rfl['number'] = np.arange(0,8,1)
df_57_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pMP057'
df_57_rfl['Konstrukt'] = 'YopO' 
df_57_rfl['Luciferase'] = 'RFL'
df_57_rfl['Normalisierung']='SptP1'

# nanoluc luciferase plate 3 pJA003, 005, 008

df_p5 = pd.read_csv(file5, skiprows=1)
df_p5 = df_p5.iloc[4:12,:12]
df_p5.columns = ['L1','L2','L3']*4
df_p5 = df_p5.reset_index(drop=True)

s13 = df_p5.iloc[:8,:3]
s14 = df_p5.iloc[:8,4:7]
s15 = df_p5.iloc[:8,8:11]

df_03_nl = pd.concat([s13])
df_03_nl['time'] = 4
df_03_nl['number'] = np.arange(0,8,1)
df_03_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA003'
df_03_nl['Konstrukt'] = 'Map' 
df_03_nl['Luciferase'] = 'NL'
df_03_nl['Normalisierung']='SptP1'

df_05_nl = pd.concat([s14])
df_05_nl['time'] = 4
df_05_nl['number'] = np.arange(0,8,1)
df_05_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA005'
df_05_nl['Konstrukt'] = 'NleE' 
df_05_nl['Luciferase'] = 'NL'
df_05_nl['Normalisierung']='SptP2'

df_08_nl = pd.concat([s15])
df_08_nl['time'] = 4
df_08_nl['number'] = np.arange(0,8,1)
df_08_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA008'
df_08_nl['Konstrukt'] = 'IpaH9.8' 
df_08_nl['Luciferase'] = 'NL'
df_08_nl['Normalisierung']='SptP2'

# red firefly luciferase plate 3 pJA003, 005, 008

df_p6 = pd.read_csv(file6, skiprows=1)
df_p6 = df_p6.iloc[4:12,:12]
df_p6.columns = ['L1','L2','L3']*4
df_p6 = df_p6.reset_index(drop=True)

s16 = df_p6.iloc[:8,:3]
s17 = df_p6.iloc[:8,4:7]
s18 = df_p6.iloc[:8,8:11]

df_03_rfl = pd.concat([s16])
df_03_rfl['time'] = 4
df_03_rfl['number'] = np.arange(0,8,1)
df_03_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA003'
df_03_rfl['Konstrukt'] = 'Map' 
df_03_rfl['Luciferase'] = 'RFL'
df_03_rfl['Normalisierung']='SptP1'

df_05_rfl = pd.concat([s17])
df_05_rfl['time'] = 4
df_05_rfl['number'] = np.arange(0,8,1)
df_05_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA005'
df_05_rfl['Konstrukt'] = 'NleE' 
df_05_rfl['Luciferase'] = 'RFL'
df_05_rfl['Normalisierung']='SptP2'

df_08_rfl = pd.concat([s18])
df_08_rfl['time'] = 4
df_08_rfl['number'] = np.arange(0,8,1)
df_08_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA008'
df_08_rfl['Konstrukt'] = 'IpaH9.8' 
df_08_rfl['Luciferase'] = 'RFL'
df_08_rfl['Normalisierung']='SptP2'

# nanoluc luciferase plate 4 pJA009, 017, 019

df_p7 = pd.read_csv(file7, skiprows=1)
df_p7 = df_p7.iloc[4:12,:12]
df_p7.columns = ['L1','L2','L3']*4
df_p7 = df_p7.reset_index(drop=True)

s19 = df_p7.iloc[:8,:3]
s20 = df_p7.iloc[:8,4:7]
s21 = df_p7.iloc[:8,8:11]

df_09_nl = pd.concat([s19])
df_09_nl['time'] = 4
df_09_nl['number'] = np.arange(0,8,1)
df_09_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA009'
df_09_nl['Konstrukt'] = 'OspB' 
df_09_nl['Luciferase'] = 'NL'
df_09_nl['Normalisierung']='SptP2'

df_17_nl = pd.concat([s20])
df_17_nl['time'] = 4
df_17_nl['number'] = np.arange(0,8,1)
df_17_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA017'
df_17_nl['Konstrukt'] = 'CT_115' 
df_17_nl['Luciferase'] = 'NL'
df_17_nl['Normalisierung']='SptP2'

df_19_nl = pd.concat([s21])
df_19_nl['time'] = 4
df_19_nl['number'] = np.arange(0,8,1)
df_19_nl['Linie'] = 'SB905 dSipA dSptP pHilA pJA019'
df_19_nl['Konstrukt'] = 'CT_223' 
df_19_nl['Luciferase'] = 'NL'
df_19_nl['Normalisierung']='SptP2'

# red firefly luciferase plate 4 pJA009, 017, 019

df_p8 = pd.read_csv(file8, skiprows=1)
df_p8 = df_p8.iloc[4:12,:12]
df_p8.columns = ['L1','L2','L3']*4
df_p8 = df_p8.reset_index(drop=True)

s22 = df_p8.iloc[:8,:3]
s23 = df_p8.iloc[:8,4:7]
s24 = df_p8.iloc[:8,8:11]

df_09_rfl = pd.concat([s22])
df_09_rfl['time'] = 4
df_09_rfl['number'] = np.arange(0,8,1)
df_09_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA009'
df_09_rfl['Konstrukt'] = 'OspB' 
df_09_rfl['Luciferase'] = 'RFL'
df_09_rfl['Normalisierung']='SptP2'

df_17_rfl = pd.concat([s23])
df_17_rfl['time'] = 4
df_17_rfl['number'] = np.arange(0,8,1)
df_17_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA017'
df_17_rfl['Konstrukt'] = 'CT_115' 
df_17_rfl['Luciferase'] = 'RFL'
df_17_rfl['Normalisierung']='SptP2'

df_19_rfl = pd.concat([s24])
df_19_rfl['time'] = 4
df_19_rfl['number'] = np.arange(0,8,1)
df_19_rfl['Linie'] = 'SB905 dSipA dSptP pHilA pJA019'
df_19_rfl['Konstrukt'] = 'CT_223' 
df_19_rfl['Luciferase'] = 'RFL'
df_19_rfl['Normalisierung']='SptP2'

df_data_rfl = pd.concat([df_SptP1_rfl, df_SptP2_rfl, df_SptP_leer_rfl, df_NK_rfl, df_53_rfl, df_57_rfl, df_03_rfl, df_05_rfl, df_08_rfl, df_09_rfl, df_17_rfl, df_19_rfl])
df_data_nl = pd.concat ([df_SptP1_nl, df_SptP2_nl, df_SptP_leer_nl, df_NK_nl, df_53_nl, df_57_nl, df_03_nl, df_05_nl, df_08_nl, df_09_nl, df_17_nl, df_19_nl])

df_data_rfl.loc[:,['L1','L2','L3', 'time']] = df_data_rfl.loc[:,['L1','L2','L3','time']].astype('float32')
df_data_nl.loc[:,['L1','L2','L3', 'time']] = df_data_nl.loc[:,['L1','L2','L3','time']].astype('float32')

df_data_rfl['L_mean'] = df_data_rfl[['L1','L2','L3']].mean(axis=1)
df_data_nl['L_mean'] = df_data_nl[['L1', 'L2', 'L3']].mean(axis=1)

df_data_rfl['L_std'] = df_data_rfl[['L1','L2','L3']].std(axis=1)
df_data_rfl['L_sem'] = sp.sem(df_data_rfl[['L1','L2','L3']], axis=1, ddof=1)

df_data_nl['L_std'] = df_data_nl[['L1','L2','L3']].std(axis=1)
df_data_nl['L_sem'] = sp.sem(df_data_nl[['L1','L2','L3']], axis=1, ddof=1)

#Verrechnung der Lumineszenz mit OD pro Kultur

od1_file = 'OD Tabelle plate 1 mean, standardabweichung, standarderror of mean.xlsx'
od2_file = 'OD Tabelle plate 2 mean, standardabweichung, standarderror of mean.xlsx'

df_od1 = pd.read_excel(od1_file)
df_od1.to_csv('CSV_OD_plate_1.csv')
df_od2 = pd.read_excel(od2_file)
df_od2.to_csv('CSV_OD_plate_2.csv')

df_od1['Linie'] = df_od1['Linie'].str.rstrip()
df_od1['Linie'] = df_od1['Linie'].str.lstrip()

df_od2['Linie'] = df_od2['Linie'].str.rstrip()
df_od2['Linie'] = df_od2['Linie'].str.lstrip ()

df_od1.columns = ['number','Linie','OD1','OD2','OD_mean','OD_std','OD_sem']
df_od2.columns = ['number','Linie','OD1','OD2','OD_mean','OD_std','OD_sem']

df_od=pd.concat([df_od1,df_od2])

df_full_data_rfl= pd.merge(df_data_rfl, df_od, on=['number','Linie'])
df_full_data_nl= pd.merge(df_data_nl, df_od, on=['number','Linie'])

df_long_L_rfl = pd.melt(df_full_data_rfl.loc[:,['L1','L2','L3','time','Linie','Konstrukt','Luciferase','Normalisierung','number']], id_vars=['time','Linie','Konstrukt','Luciferase','Normalisierung','number'], value_vars=['L1','L2','L3'], value_name='Lumineszenz')
df_long_L_nl = pd.melt(df_full_data_nl.loc[:,['L1','L2','L3','time','Linie','Konstrukt','Luciferase','Normalisierung','number']], id_vars=['time','Linie','Konstrukt','Luciferase','Normalisierung','number'], value_vars=['L1','L2','L3'], value_name='Lumineszenz')

plt.style.use('tableau-colorblind10')
df_od['OD_mean'] = df_od[['OD1','OD2']].mean(axis=1)
df_full_data_rfl['Lumineszenz_zu_OD'] = df_full_data_rfl['L_mean']/df_full_data_rfl['OD_mean']
df_full_data_nl['Lumineszenz_zu_OD'] = df_full_data_nl['L_mean']/df_full_data_nl['OD_mean']

#Normalisierung der Lumineszenzwerte

df_Normalisierungsdaten_SptP1_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='SptP1']
df_Normalisierung_SptP1daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'PK2']
df_Normalisierungsdaten_SptP1_nl['SptP1'] = ['1.314060e+07', '1.200240e+07', '9.729743e+06', '1.249769e+07', '1.005833e+07', '9.949031e+06', '1.076356e+07', '1.042029e+07']*6
df_Normalisierungsdaten_SptP1_nl.loc[:,['SptP1']] = df_Normalisierungsdaten_SptP1_nl.loc[:,['SptP1']].astype('float32')
df_Normalisierungsdaten_SptP1_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_SptP1_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_SptP1_nl['SptP1']

df_Normalisierungsdaten_SptP2_nl=df_full_data_nl.loc[df_full_data_nl['Normalisierung']=='SptP2']
df_Normalisierung_SptP2daten_nl=df_full_data_nl.loc[df_full_data_nl['Konstrukt'] == 'PK']
df_Normalisierungsdaten_SptP2_nl['SptP2'] = ['1.229171e+07', '1.288936e+07', '1.318101e+07', '1.247667e+07', '1.482219e+07', '9.304068e+06', '8.842473e+06', '6.992757e+06']*6
df_Normalisierungsdaten_SptP2_nl.loc[:,['SptP2']] = df_Normalisierungsdaten_SptP2_nl.loc[:,['SptP2']].astype('float32')
df_Normalisierungsdaten_SptP2_nl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_SptP2_nl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_SptP2_nl['SptP2']

df_normalized_data_nl = pd.concat([df_Normalisierungsdaten_SptP1_nl, df_Normalisierungsdaten_SptP2_nl])

df_normalized_data_nl_ohne_SptP1 = df_normalized_data_nl.loc[df_normalized_data_nl['Konstrukt'] != 'PK2']

colors = ["b"]
sns.set(font_scale=1.5)
plt.figure(figsize=(15, 10))
graph=sns.barplot(data=df_normalized_data_nl_ohne_SptP1, x='Konstrukt', y='Lumineszenz/OD_normalisiert', palette=colors, ci = 'sd', dodge=False)
graph.axhline(1)
plt.title('Mittlere Lumineszenz/OD normalisiert NanoLuc Luciferase (+Standardabweichung)', y=1.02)
plt.xticks(rotation=45)
plt.show()

df_Normalisierungsdaten_SptP1_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='SptP1']
df_Normalisierung_SptP1daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'PK2']
df_Normalisierungsdaten_SptP1_rfl['SptP1'] = ['2698.937307', '3110.996067', '2468.721899', '2952.965506', '2864.291921', '2385.090889', '2817.954881', '2554.075514']*6
df_Normalisierungsdaten_SptP1_rfl.loc[:,['SptP1']] = df_Normalisierungsdaten_SptP1_rfl.loc[:,['SptP1']].astype('float32')
df_Normalisierungsdaten_SptP1_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_SptP1_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_SptP1_rfl['SptP1']

df_Normalisierungsdaten_SptP2_rfl=df_full_data_rfl.loc[df_full_data_rfl['Normalisierung']=='SptP2']
df_Normalisierung_SptP2daten_rfl=df_full_data_rfl.loc[df_full_data_rfl['Konstrukt'] == 'PK']
df_Normalisierungsdaten_SptP2_rfl['SptP2'] = ['2166.161363', '2499.304796', '2527.957631', '2658.445934', '3297.060283', '2435.118829', '2251.772763', '1304.941570']*6
df_Normalisierungsdaten_SptP2_rfl.loc[:,['SptP2']] = df_Normalisierungsdaten_SptP2_rfl.loc[:,['SptP2']].astype('float32')
df_Normalisierungsdaten_SptP2_rfl ['Lumineszenz/OD_normalisiert'] = df_Normalisierungsdaten_SptP2_rfl['Lumineszenz_zu_OD']/df_Normalisierungsdaten_SptP2_rfl['SptP2']

df_normalized_data_rfl = pd.concat([df_Normalisierungsdaten_SptP1_rfl, df_Normalisierungsdaten_SptP2_rfl])

df_normalized_data_rfl_ohne_SptP1 = df_normalized_data_rfl.loc[df_normalized_data_rfl['Konstrukt'] != 'PK2']

colors = ["b"]
sns.set(font_scale=1.5)
plt.figure(figsize=(15, 10))
graph2=sns.barplot(data=df_normalized_data_rfl_ohne_SptP1, x='Konstrukt', y='Lumineszenz/OD_normalisiert', palette=colors, ci = 'sd', dodge=False) 
graph2.axhline(1)
plt.title('Mittlere Lumineszenz/OD normalisiert Red Firefly Luciferase (+Standardabweichung)', y=1.02)
plt.xticks(rotation=45)
plt.show()

