"""Convert cigar string to distribution
"""

from model_creation.data_compilation.read_forecast_data import main as read_data_main
import numpy as np
import re
from amber.utils.data_parser import seq_to_matrix
import pickle

def token_to_full_indel(indel):
    """From SelfTarget
    """
    indel_toks = indel.split('_')
    indel_type, indel_details = indel_toks[0], ''
    if len(indel_toks) > 1:
        indel_details =  indel_toks[1]
    cigar_toks = re.findall(r'([CLRDI]+)(-?\d+)', indel_details)
    details, muts = {'I':0,'D':0,'C':0}, []
    for (letter,val) in cigar_toks:
        details[letter] = eval(val)
    if len(indel_toks) > 2 or (indel_type == '-' and len(indel_toks) > 1):
        mut_toks = re.findall(r'([MNDSI]+)(-?\d+)(\[[ATGC]+\])?', indel_toks[-1])
        for (letter,val,nucl) in mut_toks:
            if nucl == '':
                nucl = '[]'
            muts.append((letter, eval(val), nucl[1:-1]))
    if indel_type[0] == '-':
        isize = 0
    else:
        isize = eval(indel_type[1:])
    return indel_type[0],isize,details, muts


def get_indel_freq(label_dict, label_oligos):
    #freq_dict = {}  # oligo -> [deletion_count, insertion_count]
    freq = []
    for oligo in label_oligos:
        this = [0,0]
        for cigar in label_dict[oligo]:
            if cigar.startswith("D"):
                this[0] += label_dict[oligo][cigar]
            elif cigar.startswith("I"):
                this[1] += label_dict[oligo][cigar]
            else:
                raise Exception("Error in cigar: %s"%cigar)
        #freq_dict[oligo] = this
        freq.append(this)
    freq = np.asarray(freq)
    freq = freq[:,0] / (freq[:,0] + freq[:,1]) # prop. of deletion / indels
    return freq


def get_1bp_insertion(label_dict, label_oligos):
    prob_1bp = []
    for oligo in label_oligos:
        this = [0,0]
        for cigar in label_dict[oligo]:
            full_indel = token_to_full_indel(cigar)
            count = label_dict[oligo][cigar]
            if full_indel[0] == "I" and full_indel[1] == 1:
                this[0] += count
            else:
                this[1] += count
        prob_1bp.append(this)
    prob_1bp = np.asarray(prob_1bp)
    prob_1bp = prob_1bp[:,0] / (prob_1bp[:,0] + prob_1bp[:,1])
    return prob_1bp


def split_train_val_test(label_oligos, test_prop=0.1, val_prop=0.1, seed=None):
    n_sample = len(label_oligos)
    test_num = int(n_sample*test_prop)
    val_num = int(n_sample*val_prop)
    np.random.seed(seed)
    idx = np.random.choice(np.arange(len(label_oligos)), len(label_oligos), replace=False)
    test_idx, val_idx, train_idx = idx[0:test_num], idx[test_num:(test_num+val_num)], idx[(test_num+val_num):]
    return test_idx, val_idx, train_idx 


def dump_pickle(seqs, labels, idx, fp=None):
    x= seqs[idx]
    if type(labels) is list:
        y = [y_[idx] for y_ in labels]
    elif type(labels) is np.ndarray:
        y = labels[idx]
    else:
        raise Exception("Error in checking labels type; must be list or np.array")

    if fp:
        pickle.dump((x,y), open(fp,"wb"), -1)

    return x, y


def main():
    oligo_dict, label_dict = read_data_main()
    label_oligos = [k for k in label_dict if k in oligo_dict]

    # compile many summary stats as labels
    indel_freq = get_indel_freq(label_dict, label_oligos)
    prob_ins1bp = get_1bp_insertion(label_dict, label_oligos)

    # compile the raw sequence one-hot encoded
    seqs = np.asarray([seq_to_matrix(oligo_dict[oligo]) for oligo in label_oligos])

    # split into test, val, train, and save pikle
    test_idx, val_idx, train_idx = split_train_val_test(label_oligos, seed=777)

    # dump to disk
    _ = dump_pickle(seqs, [indel_freq, prob_ins1bp], train_idx, fp="./data/01_train_data/forecast.train.pkl")
    _ = dump_pickle(seqs, [indel_freq, prob_ins1bp], val_idx, fp="./data/01_train_data/forecast.val.pkl")
    _ = dump_pickle(seqs, [indel_freq, prob_ins1bp], test_idx, fp="./data/01_train_data/forecast.test.pkl")
