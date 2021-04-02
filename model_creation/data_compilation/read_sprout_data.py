'''
Reading and cleaning the Leenay et. al. bam_df
'''
import os
import statistics 
import numpy as np
import pandas as pd

# Run load_sprout.R first to get 'summary_df.tsv' and 'Leenay_bam_df.csv'

######################################################
## Create final_df.csv, add cols with num of indels ##
######################################################

summary_df_path = 'summary_df.txt'
Leenay_bam_df_path = 'Leenay_bam_df.csv'
final_df_path = 'final_df.csv'

def create_clean_csv(): # clean up Leenay_bam_df.csv
    bam_df = pd.read_csv(Leenay_bam_df_path, low_memory=False)
    bam_df_ = bam_df[bam_df['reference'].notna()] #Get rid of rows where reference is NA
    bam_df_ = bam_df_.rename(columns={'Unnamed: 0': 'Number'})
    bam_df_reduce = bam_df_.drop(columns=['Number', 'guide', 'gdlrow', 'reference', 'genename']) 
    bam_files = []
    for i in range(len(bam_df_reduce)):
        row_bam_lst = []
        for j in range(len(bam_df_reduce.iloc[i].dropna())):
            row_bam_lst.append(bam_df_reduce.iloc[i].dropna()[j])
        row_bam_file = ','.join(row_bam_lst)
        bam_files.append(row_bam_file)
    bam_df_['BAM file'] = bam_files # col with name of cols without NaN
    bam_df_ = bam_df_[['Number', 'guide', 'gdlrow', 'reference', 'genename', 'BAM file']]
    return bam_df

def create_final_df():
    summary_df = pd.read_csv(summary_df_path, sep="\t")
    summary_df = summary_df[['genename', 'refseq', 'chrom', 'ranges', 'strand']]
    bam_df = create_clean_csv()
    bam_df = bam_df[['Number', 'genename', 'guide', 'BAM file']]
    final_df = summary_df.merge(bam_df, on='genename', how='inner')
    final_df = final_df.drop_duplicates()
    final_df = final_df.reset_index(drop=True)
    final_df = final_df[['Number', 'genename', 'refseq', 'chrom', 'ranges', 'strand', 'guide', 'BAM file']]
    final_df.rename(columns = {'BAM file':'bams'}, inplace = True)
    final_df.to_csv(final_df_path, index=False) #1987 rows, 8 cols

counts_dir = 'data/Sprout/counts'

def add_num_indels():
    df = pd.read_csv(final_df_path)
    for col in ['insertions', 'deletions', 'total_out']: df[col] = np.nan
    bam_lst = df['bams'].tolist()
    
    files = os.listdir(counts_dir)
    for file_num in [x for x in range(len(files)) if x != 822]: # files[822] is counts-YWHAG-01-1699.txt which only has NaNs
        filename = files[file_num]
        counts = pd.read_csv(os.path.join(counts_dir, filename), index_col=0) #make counts df
        counts = counts.drop(columns=['seq', 'seqlen', 'proba', 'total'], errors='ignore')

        complex_indices = []
        ins_indices = [] #list outcomes with 1 insertion, no deletions
        del_indices = [] #list outcomes with 1 deletion, no insertions
        for ind in range(len(counts.index)):
            if (counts.index[ind].count('I') == 1) and (counts.index[ind].count('D') == 0): ins_indices.append(ind)
            elif (counts.index[ind].count('D') == 1) and (counts.index[ind].count('I') == 0): del_indices.append(ind)

            if (counts.index[ind].count('I') > 0) and (counts.index[ind].count('D') > 0): complex_indices.append(ind)
            elif (counts.index[ind].count('I') > 1) or (counts.index[ind].count('D') > 1): complex_indices.append(ind)
            elif(counts.index[ind] == 'Other'): complex_indices.append(ind)

        complex_rows = counts.iloc[complex_indices].index
        sim_counts = counts[~counts.index.isin(complex_rows)]
        total_out = sim_counts.to_numpy().sum() # sum all values in sim_counts = simple outs [include SNVs, no variant, Other]

        ins_counts = counts.iloc[ins_indices]
        ins = ins_counts.to_numpy().sum() # sum all values in insertions df = # insertions
        
        del_counts = counts.iloc[del_indices]
        dels = del_counts.to_numpy().sum() # sum all values in deletions df = # deletions

        bam = list(counts.columns)[0] # get first bam file/col name 
        df_index = [i for i, s in enumerate(bam_lst) if bam in s]
            # find row of df with bams that contains first column name (bam file) from counts df 
        
        df.loc[df_index, 'total_out'] = total_out
        df.loc[df_index, 'insertions'] = ins
        df.loc[df_index, 'deletions'] = dels

    df.to_csv(final_df_path, index=False)

# Run:
# create_final_df()
# add_num_indels()

#######################
## Create key_df.csv ##
#######################

key_df_path = 'key_df.csv'

def get_refseq_dups(): # different genename, same refseq
    df = pd.read_csv(final_df_path)
    refseq_lst = df['refseq'].tolist()
    refseq_dups_lst = []
    for i in range(len(df)):
        refseq = refseq_lst[i]
        temp_df = df[df['refseq'] == refseq]
        genename_lst = temp_df['genename'].tolist()
        if len(list(set(genename_lst))) != 1: refseq_dups_lst.append(list(set(genename_lst)))
    norepeat_dup_lsts = []
    for elem in refseq_dups_lst:
        if elem not in norepeat_dup_lsts:
            norepeat_dup_lsts.append(elem)
    # ['CXCR4', 'CXCR4r-01', 'CXCR4r80-01', 'CXCR4r-02', 'CXCR4r80-02'] 5
    # ['LEDGF', 'LEDGFr-01', 'LEDGFr80-01', 'LEDGFr-02', 'LEDGFr80-02'] 5
    # ['CDK9',  'CDK9r-01',  'CDK9r80-01',  'CDK9r-02',  'CDK9r80-02'] 5
    return norepeat_dup_lsts

def create_key_df():
    df = pd.read_csv(final_df_path)
    key_df = df[['genename', 'refseq']]
    key_df['refseq'] = key_df['refseq'].str.upper()
    key_df = key_df.drop_duplicates().reset_index(drop=True)
    key_df['id_ending'] = np.nan
    norepeat_dup_lsts = get_refseq_dups()

    for i in range(len(key_df)):
        genename = key_df['genename'][i]
        temp_df = df[df['genename'] == genename]
        Numbers_lst = temp_df['Number']
        Numbers_lst = ['-' + str(x) for x in Numbers_lst]
        key_df['id_ending'][i] = Numbers_lst
    
    for dup_lst in norepeat_dup_lsts:
        genename = min(dup_lst, key=len) 
        extra_str_lst = [s.replace(genename, '') for s in dup_lst]
        id_ending_ = []
        for dup_num in range(len(dup_lst)):
            dup_name = dup_lst[dup_num]
            extra_str = extra_str_lst[dup_num]
            id_end = key_df.loc[key_df['genename'] == dup_name, 'id_ending'].tolist()[0]
            id_end = [extra_str + end for end in id_end]
            id_ending_ += id_end
            if extra_str != '': key_df = key_df[key_df['genename'] != dup_name]
        
        key_df.loc[key_df['genename'] == genename, 'id_ending'] = str(id_ending_)

    key_df.to_csv(key_df_path, index=False)

def get_indels_and_totalout(): #based on final_df numbers
    df = pd.read_csv(key_df_path)
    final_df = pd.read_csv(final_df_path)
    final_df = final_df[['genename', 'insertions', 'deletions', 'total_out']]

    for i in range(len(df)):
        genename = df['genename'][i]
        final_df_ = final_df[final_df['genename'] == genename]
        insertions = final_df_['insertions'].to_numpy().sum()
        deletions = final_df_['deletions'].to_numpy().sum()
        total_out = final_df_['total_out'].to_numpy().sum()
        df.loc[i, 'insertions'] = insertions
        df.loc[i, 'deletions'] = deletions
        df.loc[i, 'total_out'] = total_out
   
    df['delfreq'] = df['deletions'] / (df['deletions'] + df['insertions'])
    df.to_csv(key_df_path, index=False)
