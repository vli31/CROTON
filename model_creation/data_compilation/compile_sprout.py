'''
Generate labels from Leenay data (del_freq, 1bp_ins_proba etc.)'''

import os
import pandas as pd
import numpy as np
from scipy.stats import entropy
from model_creation.data_compilation.read_sprout_data import final_df


key_df_path = 'key_df.csv'
final_df_path = 'final_df.csv'

counts_dir = 'data/Sprout/counts'
insertions_dir = 'data/Sprout/30insertions'
master_dir = 'data/Sprout/master'

#######################
## Create master_dfs ##
#######################

def get_master_df(cts_path, ins_path, maxlen=80): # maxlen = 80 [do not consider insertions > 20 bp]
    counts = pd.read_csv(cts_path)
    counts = counts.rename(columns={ counts.columns[0]: 'cigar' })

    counts = counts.dropna()
    no_variant_inds = []
    for ind in range(len(counts)):
        if (counts.cigar.iloc[ind].count('SNV') == 1): no_variant_inds.append(ind)
        elif (counts.cigar.iloc[ind] == 'no variant'): no_variant_inds.append(ind)
    master_df = counts[~counts.index.isin(no_variant_inds)] # only deletions

    no_include_inds = []
    for ind in range(len(counts)):
        if (counts.cigar.iloc[ind].count('SNV') == 1): no_include_inds.append(ind)
        elif (counts.cigar.iloc[ind] == 'no variant'): no_include_inds.append(ind)
        
        if (counts.cigar.iloc[ind].count('I') > 0) and (counts.cigar.iloc[ind].count('D') > 0): no_include_inds.append(ind)
        elif (counts.cigar.iloc[ind].count('I') > 1) or (counts.cigar.iloc[ind].count('D') > 1): no_include_inds.append(ind)
        elif(counts.cigar.iloc[ind] == 'Other'): no_include_inds.append(ind)
    
    master_df = counts[~counts.index.isin(no_include_inds)] # only deletions
    master_df = master_df[['cigar', 'total', 'seq', 'seqlen']]

    if os.path.exists(ins_path):
        ins = pd.read_csv(ins_path)
        ins = ins.dropna()
        if len(ins) != 0:
            grouped = ins.groupby(ins.Sample)
            groups = list(set(ins.Sample.tolist()))
            ins_df = grouped.get_group(groups[0])
            ins_df = ins_df[['Allele', 'Count', 'insseq', 'seqlen']]
            for i in range(len(groups) - 1):
                sample = grouped.get_group(groups[i+1])
                sample = sample[['Allele', 'Count', 'insseq', 'seqlen']]
                ins_df = pd.merge(ins_df, sample, on=['Allele', 'insseq', 'seqlen'], how='outer').fillna(0)
                ins_df['Count'] = ins_df['Count_x'] + ins_df['Count_y']
                ins_df = ins_df[['Allele', 'Count', 'insseq', 'seqlen']]
            ins_df = ins_df[ins_df['seqlen'] <= maxlen]
            ins_df = ins_df.rename(columns={'Allele': 'cigar', 'Count':'total', 'insseq':'seq'})
            master_df = master_df.append(ins_df)
    master_df['seq'] = master_df['seq'].str.upper()
    
    return master_df

def create_master_files():
    df = pd.read_csv(key_df_path) #get key df with no replicates
    for i in range(len(df)):
        genename = df['genename'][i]
        id_ending_lst = df['id_ending'][i]
        id_ending_lst = id_ending_lst.strip('][').split(', ')
        id_ending_lst = [end.replace("'", '') for end in id_ending_lst]
        gene_paths = [genename + end + '.txt' for end in id_ending_lst]
        
        for j in range(len(id_ending_lst)):
            cts_path = os.path.join(counts_dir, 'counts-' + gene_paths[j])
            ins_path = os.path.join(insertions_dir, 'insertions-' + gene_paths[j])
            master_df_ = get_master_df(cts_path, ins_path)
            if j == 0: master_df = master_df_ # first dataframe
            else:
                master_df = master_df.merge(master_df_, on=['cigar', 'seq', 'seqlen'], how='outer')
                master_df['total_x'] = master_df['total_x'].fillna(0.)
                master_df['total_y'] = master_df['total_y'].fillna(0.)
                master_df['total'] = master_df['total_x'] + master_df['total_y']
                master_df = master_df.drop(columns=['total_x', 'total_y'])
        master_df = master_df.rename(columns={'total': 'count'})
        master_df['seq'] = master_df['seq'].str.upper()
        master_name = 'master-' + genename + '.txt'
        master_df.to_csv(os.path.join(master_dir, master_name), index=False)

################################
## Add statcols to key_df.csv ##
################################

def add_proba_1bp_insertion():
    df = pd.read_csv(key_df_path)
    df['prob_1bpins']: df[col] = np.nan
    for i in range(len(df)):
        genename = df['genename'][i]
        total_indels = df['insertions'][i] + df['deletions'][i]

        master_name = 'master-' + genename + '.txt'
        master_df = pd.read_csv(os.path.join(master_dir, master_name))
        master_df = master_df[master_df['seqlen'] == 61] #master only has 1 bp insertions now
        total_1bpins = master_df['count'].to_numpy().sum()

        df.loc[i, 'prob_1bpins'] = total_1bpins / total_indels

    df.to_csv(key_df_path, index=False)
#proba_1bp_insertion()

def add_proba_1bp_deletion():
    df = pd.read_csv(key_df_path)
    df['prob_1bpdel'] = np.nan
    for i in range(len(df)):
        genename = df['genename'][i]
        total_indels = df['insertions'][i] + df['deletions'][i]

        master_name = 'master-' + genename + '.txt'
        master_df = pd.read_csv(os.path.join(master_dir, master_name))
        master_df = master_df[master_df['seqlen'] == 59] #master only has 1 bp insertions now
        total_1bpdel = master_df['count'].to_numpy().sum()

        df.loc[i, 'prob_1bpdel'] = total_1bpdel / total_indels

    df.to_csv(key_df_path, index=False)
#proba_1bp_deletion()

def add_frameshift_freq():
    df = pd.read_csv(key_df_path)
    onemod3_freq, twomod3_freq = [], []
    for col in ['onemod3_freq', 'twomod3_freq', 'frameshift_freq']: df[col] = np.nan
    for i in range(len(df)):
        genename = df['genename'][i]
        id_ending_lst = df['id_ending'][i]
        id_ending_lst = id_ending_lst.strip('][').split(', ')
        id_ending_lst = [end.replace("'", '') for end in id_ending_lst]
        gene_paths = [genename + end + '.txt' for end in id_ending_lst]
        total, total_onemod3, total_twomod3 = 0, 0, 0
        for j in range(len(id_ending_lst)):
            cts_path = os.path.join(counts_dir, 'counts-' + gene_paths[j])
            counts = pd.read_csv(cts_path)
            counts = counts.rename(columns={ counts.columns[0]: 'cigar' })
            exclude_inds = []
            for ind in range(len(counts)):
                if (counts.cigar.iloc[ind].count('SNV') == 1): exclude_inds.append(ind)
                elif (counts.cigar.iloc[ind] == 'no variant'): exclude_inds.append(ind)
                elif(counts.cigar.iloc[ind] == 'Other'): exclude_inds.append(ind)
                elif (counts.cigar.iloc[ind].count('I') > 0) and (counts.cigar.iloc[ind].count('D') > 0): exclude_inds.append(ind)
                elif (counts.cigar.iloc[ind].count('I') > 1) or (counts.cigar.iloc[ind].count('D') > 1): exclude_inds.append(ind)
            counts = counts[~counts.index.isin(exclude_inds)] #get rid of complex variants and no variants, only indels including those outside of 60bp range 
            
            seqlens = []
            for cigar in counts.cigar:
                if 'D' in cigar:
                    cigar_ = cigar.split(':')
                    dellen = cigar_[1].replace('D', '')
                    seqlens.append(60-int(dellen))
                if 'I' in cigar:
                    cigar_ = cigar.split(':')
                    inslen = cigar_[1].replace('I', '')
                    seqlens.append(60+int(inslen))
            counts['seqlen'] = seqlens
            counts['mod3'] = (counts['seqlen'] - 60) % 3
            counts_onemod3 = counts[counts['mod3'] == 1]
            counts_twomod3 = counts[counts['mod3'] == 2]

            total += counts['total'].to_numpy().sum()
            total_onemod3 += counts_onemod3['total'].to_numpy().sum()
            total_twomod3 += counts_twomod3['total'].to_numpy().sum()
        
        onemod3_freq.append(total_onemod3 / total) #computing total over all counts dfs for gene
        twomod3_freq.append(total_twomod3 / total)
    
    df['onemod3_freq'] = onemod3_freq
    df['twomod3_freq'] = twomod3_freq
    df['frameshift_freq'] = df['onemod3_freq'] + df['twomod3_freq']
    df.to_csv(key_df_path, index=False)
#frameshift_freq()
