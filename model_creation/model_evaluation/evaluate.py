
import os
import pickle
import numpy as np
import pandas as pd
from tensorflow.python.keras.models import load_model
from sklearn.metrics import roc_auc_score
import scipy.stats as ss

TASK_IDENTIFIER = {'delfreq':0, 'prob_1bpins': 1, 'prob_1bpdel': 2, 'onemod3_freq': 3, 'twomod3_freq': 4, 'frameshift_freq': 5}
statlst = ['delfreq','prob_1bpins','prob_1bpdel','onemod3_freq','twomod3_freq','frameshift_freq']

def one_hot_encode(seq, base_map): #takes upper case seq
    mapping = dict(zip(base_map, range(4)))
    split_seq = seq.split('N')
    map_seq = [mapping[i] for i in split_seq[0]]
    final_seq = np.eye(4)[map_seq]
    for n in range(len(split_seq) - 1):
        N_arr = np.array([0.25, 0.25, 0.25, 0.25])
        final_seq = np.vstack((final_seq, N_arr))
        map_seq = [mapping[i] for i in split_seq[n + 1]]
        final_seq_ = np.eye(4)[map_seq]
        final_seq = np.vstack((final_seq, final_seq_))
    return final_seq

#Function only works with sequences of the same length
def one_hot_encode_lst(seq_col, len_, base_map, save):
    number = 1
    stack = one_hot_encode(seq_col[0], base_map)
    for i in range(1, len(seq_col)):
        seq = one_hot_encode(seq_col[i], base_map)
        stack = np.concatenate([stack, seq])
        number += 1
    stack = np.reshape(stack, (number, len_, 4))
    if save:
        output_path = "data/npy/encoded_%ibp_%s" % (len_, base_map)
        np.save(output_path, stack, allow_pickle=True)
    return stack

def get_binary_cols(df, statlst=statlst):
    for stat in statlst:
        statlst = df[stat].tolist()
        stat_median = np.median(statlst)
        binary_statlst = np.where(statlst>=stat_median, 1.0, statlst)
        binary_statlst = np.where(statlst<stat_median, 0.0, binary_statlst)
        df.loc[:,stat+'_binary'] = list(binary_statlst)
    
    return df

def get_croton_obs(dataset):
    if dataset == 'forecast':
        x_test, y_test = pickle.load(open('./data/data/Forecast/test_.pkl', 'rb'))
        obs_df = pd.DataFrame({'delfreq':y_test[0], 'prob_1bpins':y_test[1], 'prob_1bpdel':y_test[2],
            'onemod3_freq':y_test[6], 'twomod3_freq':y_test[7], 'frameshift_freq':y_test[8]})
    
    if dataset == 'sprout':
        obs_df = pd.read_csv('key_df.csv')
        obs_df['refseq'] = obs_df['refseq'].str.upper()
        x_test = one_hot_encode_lst(obs_df['refseq'], 60, 'ACGT', False)

    obs_df = get_binary_cols(obs_df)
    
    return x_test, obs_df

def get_croton_pred(stat, x_test, croton_path='CROTON.h5'):
    pred_inx = TASK_IDENTIFIER[stat] 
    model = load_model(croton_path)
    pred = model.predict(x_test)
    pred = pred[:,pred_inx]
    pred_median = np.median(pred)
    pred_binary = np.where(pred>=pred_median, 1.0, pred)
    pred_binary = np.where(pred<pred_median, 0.0, pred_binary)

    return pred, pred_binary

def get_aucroc_stats(dataset): # *****ONLY CROTON SUPPORTS dataset = 'forecast'*****
    x_test, obs_df = get_croton_obs(dataset)
    
    modelpred_stats = {'df_label':[], 'stat':[], 'auc':[], 'pearson':[], 'kendall':[]}
    for i in range(len(statlst)):
        modelpred_stats['df_label'].append('croton_' + dataset)
        stat = statlst[i]
        obs, obs_binary = obs_df[stat], obs_df[stat+'_binary']
        
        # Get pred
        pred, pred_binary = get_croton_pred(stat, x_test)

        # Get aucroc statistics
        auc = roc_auc_score(obs_binary, pred)
        modelpred_stats['auc'].append(auc)

        # Get corr statistics
        modelpred_stats['pearson'].append(round(ss.pearsonr(obs, pred)[0], 6))
        modelpred_stats['kendall'].append(round(ss.kendalltau(obs, pred)[0], 6))
    
    # Put model prediction statistics into a dataframe
    modelpred_stats_df = pd.DataFrame.from_dict(modelpred_stats)
    try: 
        stats_df = pd.read_csv('modelpred_stats.csv') 
        modelpred_stats_df = pd.concat([stats_df, modelpred_stats_df], ignore_index=True) 
    except:
        pass
    modelpred_stats_df.to_csv('modelpred_stats.csv', index=False)

