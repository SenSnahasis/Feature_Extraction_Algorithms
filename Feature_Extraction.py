#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import itertools
from collections import Counter


# In[2]:


data = pd.read_csv('initial_input_data.csv')
data.head()


# In[3]:


data = data.drop(['Entry'], axis = 1)
data.head()


# In[4]:


data.info()


# In[5]:


df=data.dropna()
df.info()


# In[6]:


df


# In[7]:


amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
               'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


# In[8]:


count = 0;
for index, row in df.iterrows():
    valid = True
    
    for amino in row['Sequence']:
        if (amino not in amino_acids):
            valid = False
            break
            
    if not valid:       
        df.drop(index, inplace=True)
        count += 1
        
print(f'{count} data dropped')


# In[9]:


df.info()


# In[13]:


count = 0;
for index, row in df.iterrows():
    valid = True
    str = ''
    for i in row['Protein_families']:
        str = str+i;
        
    if('superfamily' in str):
        valid = False

    if not valid:       
        df.drop(index, inplace=True)
        count += 1
        
    print("yes")
        
print(f'{count} data dropped')


# In[16]:


df.info()


# In[17]:


df.to_excel("output.xlsx")


# In[18]:


df


# In[19]:


new_data = pd.read_csv('final_input.csv')


# In[20]:


new_data


# In[25]:


new_data.info()


# In[21]:


new_data = new_data.drop(['ID'], axis = 1)
new_data.head()


# In[22]:


new_data = new_data.drop(['Features'], axis = 1)
new_data.head()


# In[23]:


new_data = new_data.drop(['Protein_names'], axis = 1)
new_data.head()


# In[26]:


df2=[]
for index,row in new_data.iterrows():
    df2.append(row['Sequence'])


# In[27]:


df2


# # n-gram

# In[24]:


def generatePossiblePatterns(sequence, n):
    list_possible_patterns = []

    for pos in range(0, len(sequence) - n):

        for pattern_length in range(n):
            pattern = sequence[pos: pos+pattern_length+1]

            if pattern not in list_possible_patterns:
                list_possible_patterns.append(pattern)

    return list_possible_patterns

def findOccuranceProbabilities(sequence, n):
    list_occurance_prob = []

    list_possible_patterns = generatePossiblePatterns(sequence, n)

    for pattern in list_possible_patterns:
        match_count = sum(1 for i in range(len(sequence))
                          if sequence.startswith(pattern, i))

        occurance_prob = match_count / (len(sequence) - n - 1)
        list_occurance_prob.append(occurance_prob)

    return list_occurance_prob

def find_n_gram(n, sequence):

    list_occurance_prob = findOccuranceProbabilities(sequence, n)

    mean = np.mean(list_occurance_prob)
    sd = np.std(list_occurance_prob)
    cv = sd/mean

    return (mean, sd, cv)


# In[28]:


ng2=pd.DataFrame(columns=['serial_no','mean_2','sd_2','cv_2'])
for i in range(4850):
    a=find_n_gram(2, df2[i])
    ng2.loc[i,['serial_no']] = i
    ng2.loc[i,['mean_2']] = a[0]
    ng2.loc[i,['sd_2']] = a[1]
    ng2.loc[i,['cv_2']] = a[2]
    
    print(i)


# In[29]:


ng2


# In[30]:


ng3=pd.DataFrame(columns=['serial_no','mean_3','sd_3','cv_3'])
for i in range(4850):
    a=find_n_gram(3, df2[i])
    ng3.loc[i,['serial_no']] = i
    ng3.loc[i,['mean_3']] = a[0]
    ng3.loc[i,['sd_3']] = a[1]
    ng3.loc[i,['cv_3']] = a[2]

    print(i)


# In[31]:


feature=pd.merge(ng2,ng3,left_on='serial_no',right_on='serial_no')


# In[33]:


ng4=pd.DataFrame(columns=['serial_no','mean_4','sd_4','cv_4'])
for i in range(4850):
    a=find_n_gram(4, df2[i])
    ng4.loc[i,['serial_no']] = i
    ng4.loc[i,['mean_4']] = a[0]
    ng4.loc[i,['sd_4']] = a[1]
    ng4.loc[i,['cv_4']] = a[2]
    
    print(i)


# In[34]:


feature=pd.merge(feature,ng4,left_on='serial_no',right_on='serial_no')


# In[36]:


ng5=pd.DataFrame(columns=['serial_no','mean_5','sd_5','cv_5'])
for i in range(4850):
    a=find_n_gram(5, df2[i])
    ng5.loc[i,['serial_no']] = i
    ng5.loc[i,['mean_5']] = a[0]
    ng5.loc[i,['sd_5']] = a[1]
    ng5.loc[i,['cv_5']] = a[2]
    
    print(i)


# In[37]:


feature=pd.merge(feature,ng5,left_on='serial_no',right_on='serial_no')


# In[38]:


feature


# # 6-letter

# In[70]:


properties = {
    'Tag': amino_acids,

    'Protein_Name': ['Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid',
                     'Phenylalanine', 'Glycine', 'Histidine', 'Isoleucine',
                     'Lysine', 'Leucine', 'Methionine', 'Asparagine',
                     'Proline', 'Glutamine', 'Arginine', 'Serine',
                     'Threonine', 'Valine', 'Tryptophan', 'Tyrosine'],

    'Molecular_Weight': [89.1, 121.16, 133.11, 147.13, 165.19,
                         75.07, 155.16, 131.18, 146.19, 131.18,
                         149.21, 132.12, 115.13, 146.15, 174.2,
                         105.09, 119.12, 117.15, 204.23, 181.19],

    'IsoElectric_Point': [6, 5.07, 2.77, 3.22, 5.48,
                          5.97, 7.59, 6.02, 9.74, 5.98,
                          5.74, 5.41, 6.3, 5.65, 10.76,
                          5.68, 5.6, 5.96, 5.89, 5.66],

    'Hydropathy_Property': ['Hydrophobic', 'Hydrophobic', 'Hydrophilic', 'Neutral',
                            'Very Hydrophobic', 'Neutral', 'Hydrophilic', 'Very Hydrophobic',
                            'Hydrophilic', 'Very Hydrophobic', 'Very Hydrophobic', 'Neutral',
                            'Hydrophilic', 'Neutral', 'Hydrophilic', 'Neutral',
                            'Neutral', 'Very Hydrophobic', 'Very Hydrophobic', 'Hydrophobic'],

    'Hydropathy_Label': [-1, -1, 1, 0, -2,
                         0, 1, -2, 1, -2,
                         -2, 0, 1, 0, 1,
                         0, 0, -2, -2, -1],

    '6_Letter_Encoding': ['e4', 'e3', 'e2', 'e2', 'e6',
                          'e4', 'e1', 'e5', 'e1', 'e5',
                          'e5', 'e2', 'e4', 'e2', 'e1',
                          'e4', 'e4', 'e5', 'e6', 'e6']
}


properties = pd.DataFrame(properties)


# In[40]:


import math
count = 0

def perform6LetterEncoding(sequence):
    for amino in sequence:
        exchange_code = properties.loc[properties['Tag'] == amino]['6_Letter_Encoding'].values[0]
        
        sequence = sequence.replace(amino, exchange_code)
    return sequence


def perform6LetterExchange(sequence):

    sequence = perform6LetterEncoding(sequence)

    mean, sd, cv = find_n_gram(4, sequence)

    sigmoidmean = 1 / (1 + math.exp(-mean))
    sigmoidsd = 1 / (1 + math.exp(-sd))
    sigmoidcv = 1 / (1 + math.exp(-cv))

    return (sigmoidmean, sigmoidsd, sigmoidcv)


# In[41]:


sixletter=pd.DataFrame(columns=['serial_no','sigmoidmean','sigmoidsd','sigmoidcv'])
for i in range(4850):
    a=perform6LetterExchange(df2[i])
    sixletter.loc[i,['serial_no']] = i
    sixletter.loc[i,['sigmoidmean']] = a[0]
    sixletter.loc[i,['sigmoidsd']] = a[1]
    sixletter.loc[i,['sigmoidcv']] = a[2]
    
    print(i)


# In[42]:


feature=pd.merge(feature,sixletter,left_on='serial_no',right_on='serial_no')


# # Positional Molecular Weight

# In[43]:


def findPositionalMolecularWeight(amino, pos):
    return properties.loc[properties['Tag'] == amino]['Molecular_Weight'].values * pos


# In[47]:


mole=pd.DataFrame(columns=['serial_no','Pos_Mole_Weight'])
for i in range(4850):
    pos = 1
    ans=0
    for amino in df2[i]:
        ans=ans+findPositionalMolecularWeight(amino, pos)
        pos=pos+1
    ans=ans/new_data.Length[i]
    
    mole.loc[i,['serial_no']] = i
    mole.loc[i,['Pos_Mole_Weight']] = ans[0]
    print(i)


# In[48]:


feature=pd.merge(feature,mole,left_on='serial_no',right_on='serial_no')


# # Positional Isoelectric Weight

# In[49]:


def findPositionalIsoelectricWeight(amino, pos):
    return properties.loc[properties['Tag'] == amino]['IsoElectric_Point'].values * pos


# In[51]:


iso=pd.DataFrame(columns=['serial_no','Pos_Iso_Weight'])
for i in range(4850):
    pos = 1
    ans=0
    for amino in df2[i]:
        ans=ans+findPositionalIsoelectricWeight(amino, pos)
        pos=pos+1
    ans=ans/new_data.Length[i]
    
    iso.loc[i,['serial_no']] = i
    iso.loc[i,['Pos_Iso_Weight']] = ans[0]
    print(i)


# In[52]:


feature=pd.merge(feature,iso,left_on='serial_no',right_on='serial_no')


# In[53]:


feature


# # Distance Encoding

# In[54]:


def encodeDistances(sequence):
    encoded_distances = {}

    for amino in amino_acids:
        encoded_distances[f'{amino}'] = []

    for index, amino in enumerate(sequence):
        # dict[amino] =
        #               0th index of list = prev occurance of the amino
        #               from 1st index -> store index - dict[amino][0]
        if len(encoded_distances[f'{amino}']) == 0:
            encoded_distances[f'{amino}'].append(index+1)
        else:
            prev_distance = encoded_distances[f'{amino}'][0] - 1
            encoded_distances[f'{amino}'].append(index - prev_distance)

    for amino in amino_acids:
        if len(encoded_distances[f'{amino}']) > 1:
            encoded_distances[f'{amino}'].pop(0)

    return encoded_distances


def applyMeanNormalization(encoded_list, mean):
    if(len(encoded_list)==0):
        return encoded_list
    else:
        min_distance = min(encoded_list)
        max_distance = max(encoded_list)

        if (len(encoded_list) == 1):
            range_distance = 1
        else:
            range_distance = max_distance - min_distance

        for index, distance in enumerate(encoded_list):
            encoded_list[index] = (distance - min_distance)/range_distance

        return encoded_list


def performDistanceBasedEncoding(sequence):

    n_mean = {}
    n_SD = {}
    n_CV = {}

    distances = encodeDistances(sequence)

    for amino in amino_acids:
        # find mean
        mean = np.mean(distances[f'{amino}'])

        distances[f'{amino}'] = applyMeanNormalization(
            distances[f'{amino}'], mean)

        # find the normalized values
        n_mean[f'{amino}'] = np.mean(distances[f'{amino}'])
        n_SD[f'{amino}'] = np.std(distances[f'{amino}'])

        if n_mean[f'{amino}'] != 0:
            n_CV[f'{amino}'] = n_SD[f'{amino}'] / n_mean[f'{amino}']
        else:
            n_CV[f'{amino}'] = 0

    return (n_mean, n_SD, n_CV)


# In[56]:


distance=pd.DataFrame(columns=['A_nmean','C_nmean','D_nmean','E_nmean','F_nmean','G_nmean','H_nmean','I_nmean','K_nmean','L_nmean','M_nmean','N_nmean','P_nmean','Q_nmean','R_nmean','S_nmean','T_nmean','V_nmean','W_nmean','Y_nmean','A_nsd','C_nsd','D_nsd','E_nsd','F_nsd','G_nsd','H_nsd','I_nsd','K_nsd','L_nsd','M_nsd','N_nsd','P_nsd','Q_nsd','R_nsd','S_nsd','T_nsd','V_nsd','W_nsd','Y_nsd','A_ncv','C_ncv','D_ncv','E_ncv','F_ncv','G_ncv','H_ncv','I_ncv','K_ncv','L_ncv','M_ncv','N_ncv','P_ncv','Q_ncv','R_ncv','S_ncv','T_ncv','V_ncv','W_ncv','Y_ncv'])
for i in range(4850):
    a=performDistanceBasedEncoding(df2[i])
    
    distance.loc[i,['serial_no']] = i
    distance.loc[i,['A_nmean']] = a[0]['A']
    distance.loc[i,['C_nmean']] = a[0]['C']
    distance.loc[i,['D_nmean']] = a[0]['D']
    distance.loc[i,['E_nmean']] = a[0]['E']
    distance.loc[i,['F_nmean']] = a[0]['F']
    distance.loc[i,['G_nmean']] = a[0]['G']
    distance.loc[i,['H_nmean']] = a[0]['H']
    distance.loc[i,['I_nmean']] = a[0]['I']
    distance.loc[i,['K_nmean']] = a[0]['K']
    distance.loc[i,['L_nmean']] = a[0]['L']
    distance.loc[i,['M_nmean']] = a[0]['M']
    distance.loc[i,['N_nmean']] = a[0]['N']
    distance.loc[i,['P_nmean']] = a[0]['P']
    distance.loc[i,['Q_nmean']] = a[0]['Q']
    distance.loc[i,['R_nmean']] = a[0]['R']
    distance.loc[i,['S_nmean']] = a[0]['S']
    distance.loc[i,['T_nmean']] = a[0]['T']
    distance.loc[i,['V_nmean']] = a[0]['V']
    distance.loc[i,['W_nmean']] = a[0]['W']
    distance.loc[i,['Y_nmean']] = a[0]['Y']
    
    distance.loc[i,['A_nsd']] = a[1]['A']
    distance.loc[i,['C_nsd']] = a[1]['C']
    distance.loc[i,['D_nsd']] = a[1]['D']
    distance.loc[i,['E_nsd']] = a[1]['E']
    distance.loc[i,['F_nsd']] = a[1]['F']
    distance.loc[i,['G_nsd']] = a[1]['G']
    distance.loc[i,['H_nsd']] = a[1]['H']
    distance.loc[i,['I_nsd']] = a[1]['I']
    distance.loc[i,['K_nsd']] = a[1]['K']
    distance.loc[i,['L_nsd']] = a[1]['L']
    distance.loc[i,['M_nsd']] = a[1]['M']
    distance.loc[i,['N_nsd']] = a[1]['N']
    distance.loc[i,['P_nsd']] = a[1]['P']
    distance.loc[i,['Q_nsd']] = a[1]['Q']
    distance.loc[i,['R_nsd']] = a[1]['R']
    distance.loc[i,['S_nsd']] = a[1]['S']
    distance.loc[i,['T_nsd']] = a[1]['T']
    distance.loc[i,['V_nsd']] = a[1]['V']
    distance.loc[i,['W_nsd']] = a[1]['W']
    distance.loc[i,['Y_nsd']] = a[1]['Y']
    
    distance.loc[i,['A_ncv']] = a[2]['A']
    distance.loc[i,['C_ncv']] = a[2]['C']
    distance.loc[i,['D_ncv']] = a[2]['D']
    distance.loc[i,['E_ncv']] = a[2]['E']
    distance.loc[i,['F_ncv']] = a[2]['F']
    distance.loc[i,['G_ncv']] = a[2]['G']
    distance.loc[i,['H_ncv']] = a[2]['H']
    distance.loc[i,['I_ncv']] = a[2]['I']
    distance.loc[i,['K_ncv']] = a[2]['K']
    distance.loc[i,['L_ncv']] = a[2]['L']
    distance.loc[i,['M_ncv']] = a[2]['M']
    distance.loc[i,['N_ncv']] = a[2]['N']
    distance.loc[i,['P_ncv']] = a[2]['P']
    distance.loc[i,['Q_ncv']] = a[2]['Q']
    distance.loc[i,['R_ncv']] = a[2]['R']
    distance.loc[i,['S_ncv']] = a[2]['S']
    distance.loc[i,['T_ncv']] = a[2]['T']
    distance.loc[i,['V_ncv']] = a[2]['V']
    distance.loc[i,['W_ncv']] = a[2]['W']
    distance.loc[i,['Y_ncv']] = a[2]['Y']
    
    print(i)


# In[57]:


feature=pd.merge(feature,distance,left_on='serial_no',right_on='serial_no')


# In[58]:


feature


# # Hydropathy Property

# In[67]:


def findHydropathyComposition(data_row, index):
    freq_Hpho = 0
    freq_Hphi = 0
    freq_Neu = 0

    for amino in amino_acids:
        freq_amino_in_sequence = data_row[amino]
        amino_hydropathy_type = properties.loc[properties['Tag']
                                               == amino]['Hydropathy_Label'].values

        if amino_hydropathy_type < 0:
            freq_Hpho += freq_amino_in_sequence
        elif amino_hydropathy_type > 0:
            freq_Hphi += freq_amino_in_sequence
        else:
            freq_Neu += freq_amino_in_sequence

    length = len(data_row['Sequence'])

    freq_Hpho = (freq_Hpho * 100) / length
    freq_Hphi = (freq_Hphi * 100) / length
    freq_Neu = (freq_Neu * 100) / length

    return (freq_Hpho, freq_Hphi, freq_Neu)


def findHydropathyTransmission(sequence):

    freq_Relative_Hydropathy = {
        'Hpho_Hpho': 0,
        'Hpho_Hphi': 0,
        'Hpho_Neu': 0,
        'Hphi_Hpho': 0,
        'Hphi_Hphi': 0,
        'Hphi_Neu': 0,
        'Neu_Hpho': 0,
        'Neu_Hphi': 0,
        'Neu_Neu': 0
    }

    for pos, amino in enumerate(sequence):

        length = len(sequence)

        if pos == length - 1:
            break

        next_amino = sequence[pos+1]

        amino_hydropathy_type = properties.loc[properties['Tag']
                                               == amino]['Hydropathy_Label'].values
        next_amino_hydropathy_type = properties.loc[properties['Tag']
                                                    == next_amino]['Hydropathy_Label'].values

        if amino_hydropathy_type < 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Hpho_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Hpho_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Hpho_Neu'] += 1

        if amino_hydropathy_type > 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Hphi_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Hphi_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Hphi_Neu'] += 1

        if amino_hydropathy_type == 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Neu_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Neu_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Neu_Neu'] += 1

    freq_Relative_Hydropathy['Hpho_Hpho'] = (
        freq_Relative_Hydropathy['Hpho_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Hpho_Hphi'] = (
        freq_Relative_Hydropathy['Hpho_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Hpho_Neu'] = (
        freq_Relative_Hydropathy['Hpho_Neu'] * 100) / length

    freq_Relative_Hydropathy['Hphi_Hpho'] = (
        freq_Relative_Hydropathy['Hphi_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Hphi_Hphi'] = (
        freq_Relative_Hydropathy['Hphi_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Hphi_Neu'] = (
        freq_Relative_Hydropathy['Hphi_Neu'] * 100) / length

    freq_Relative_Hydropathy['Neu_Hpho'] = (
        freq_Relative_Hydropathy['Neu_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Neu_Hphi'] = (
        freq_Relative_Hydropathy['Neu_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Neu_Neu'] = (
        freq_Relative_Hydropathy['Neu_Neu'] * 100) / length

    return freq_Relative_Hydropathy


# In[71]:


hydro=pd.DataFrame(columns=['serial_no','Hpho_Hpho', 'Hpho_Hphi', 'Hpho_Neu', 'Hphi_Hpho', 'Hphi_Hphi', 'Hphi_Neu', 'Neu_Hpho', 'Neu_Hphi', 'Neu_Neu'])
for i in range(4850):
    m=findHydropathyTransmission(df2[i])
    
    hydro.loc[i,['serial_no']] = i
    hydro.loc[i,['Hpho_Hpho']] = m['Hpho_Hpho']
    hydro.loc[i,['Hpho_Hphi']] = m['Hpho_Hphi']
    hydro.loc[i,['Hpho_Neu']] = m['Hpho_Neu']
    hydro.loc[i,['Hphi_Hpho']] = m['Hphi_Hpho']
    hydro.loc[i,['Hphi_Hphi']] = m['Hphi_Hphi']
    hydro.loc[i,['Hphi_Neu']] = m['Hphi_Neu']
    hydro.loc[i,['Neu_Hpho']] = m['Neu_Hpho']
    hydro.loc[i,['Neu_Hphi']] = m['Neu_Hphi']
    hydro.loc[i,['Neu_Neu']] = m['Neu_Neu']
    
    print(i)


# In[72]:


feature=pd.merge(feature,hydro,left_on='serial_no',right_on='serial_no')


# In[73]:


feature


# In[74]:


Protein_families=[]
for index,row in new_data.iterrows():
    Protein_families.append(row['Protein_families'])


P_family=pd.DataFrame(columns=['serial_no','Protein_families'])
for i in range(4850):
    str = Protein_families[i]
    s = ""
    for j in range(len(str)):
        if(str[j]==','):
            break
        s = s+str[j]
    P_family.loc[i,['serial_no']] = i
    P_family.loc[i,['Protein_families']] = s
    
    print(i)


# In[75]:


feature=pd.merge(feature,P_family,left_on='serial_no',right_on='serial_no')


# In[76]:


feature.to_excel("final_output.xlsx")


# In[ ]:




