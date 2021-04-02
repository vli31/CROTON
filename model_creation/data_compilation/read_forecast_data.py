"""data parser for reading FORECasT training data
"""

import os
import sys
from tqdm import tqdm
import pandas as pd
import numpy as np
import h5py
from collections import defaultdict
import pickle

# Download data from: 

def read_single_outcome(fp, label_dict, token):
    oligo_id = None
    with open(fp, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("@@@"): # a new oligo
                oligo_id = line.lstrip("@")
                continue
            ele = line.split("\t")
            count = int(ele[1])
            if token: 
                outcome_id = ele[0]
                label_dict[oligo_id][outcome_id] += count
            else: 
                seq = ele[2]
                label_dict[oligo_id][seq] += count
    return label_dict


def read_outcomes(par_dir, token, cell_line="K562", replicate=None, dpi="DPI7", coverage="800x"):
    subdirs = []
    for x in os.listdir(par_dir):
        try:
            _, _, _, this_cell, this_covg, this_inft, this_dpi = x.split("_")
        except ValueError:
            continue
        if replicate and not (this_inft.endswith(replicate)):
            continue
        if cell_line and this_cell != cell_line:
            continue
        if coverage and this_covg != coverage:
            continue
        if dpi and this_dpi != dpi:
            continue
        subdirs.extend([os.path.join(par_dir, x, y) for y in os.listdir(os.path.join(par_dir, x))] )
    fp_list = [os.path.join(x, p) for x in subdirs for p in os.listdir(x) if p.endswith("_processedindels.txt")]
    
    # label_dict: oligo_id -> cigar_string -> count
    label_dict = defaultdict(lambda: defaultdict(int))
    for fp in tqdm(fp_list):
        label_dict = read_single_outcome(fp, label_dict, token)
    return label_dict


def read_oligo_seq(grna_fp):
    """This function will do two things:
        - read in target sequence, and reverse complement it if needed
        - align the cut sites (=PAM upstream 3nt), pad Ns if necessary

    Note
    ------
    When considering cut and PAM sites, upstream 3nt needs to consider the strandness
    """
    oligo_df = pd.read_csv(grna_fp, sep="\t")
    n_rev = 0
    oligo_df['TargetSequence_forward_centered'] = ''
    for i in range(oligo_df.shape[0]):
        seq = get_recenter_pam(oligo_df['TargetSequence'][i],
                pam_index=oligo_df['PAM Index'][i],
                strand=oligo_df['Strand'][i],
                out_len=60)
        if oligo_df['Strand'][i] == 'REVERSE':
            seq = reverse_complement(seq)
            n_rev += 1
        oligo_df.at[i, 'TargetSequence_forward_centered'] = seq

    # print("total rev=%i"%n_rev)
    oligo_df.index = oligo_df.ID
    oligo_dict = oligo_df['TargetSequence_forward_centered'].to_dict()
    # print(oligo_dict)
    return oligo_dict


def get_recenter_pam(seq, pam_index, strand, out_len=60):
    flanking_len = out_len // 2
    cut_index = pam_index - 3 if strand=='FORWARD' else pam_index + 3
    left = seq[max(0, cut_index-flanking_len) : cut_index]
    right = seq[cut_index : min(len(seq), cut_index+flanking_len)]
    if len(left) < flanking_len:
        left = 'N'*(flanking_len - len(left)) + left
    if len(right) < flanking_len:
        right = right + 'N'*(flanking_len - len(right))
    return left + right


def reverse_complement(seq):
    letter_match = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_complement = "".join(letter_match[b] for b in reversed(seq))
    return reverse_complement

def main(token):
    oligo_dict = read_oligo_seq(grna_fp="/mnt/home/zzhang/workspace/src/cripsr-repair/resources/grna-target_list-pub.txt")
    label_dict = read_outcomes(par_dir="/mnt/home/zzhang/workspace/src/cripsr-repair/data/fa_2018_nbt", cell_line="K562", replicate=None,
            dpi="DPI7", coverage="800x", token=token)
    return oligo_dict, label_dict
