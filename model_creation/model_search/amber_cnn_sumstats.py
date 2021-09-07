#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Search architecture for CRISPR outcomes prediction based on FORECasT
"""

# ZZ, Sept 18, 2020
# Revised Nov 28, 2020 : more control over arguments
# Revised Jan. 21, 2021 : more labels to predict

from amber import Amber
from amber.architect import ModelSpace, Operation
from amber.modeler import KerasResidualCnnBuilder
from amber.utils import run_from_ipython
import numpy as np
import pandas as pd
import pickle
import os
import sys
import scipy.stats as ss
import argparse
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.utils import plot_model
from amber.plots import plot_training_history

TASK_IDENTIFIER = {'del_freq': 0, '1ins_freq': 1, '1del_freq': 2, 'avgdel_len': 3, 'avgins_len': 4, 'entropy': 5,
                   'onemod3_freq': 6, 'twomod3_freq': 7, 'frameshift_freq': 8}


def load_pickle_data(pickle_fp, tasks):
    """Refer to data_compilation for how sequences and labels are compiled
    """
    data = pickle.load(open(pickle_fp, "rb"))
    x = data[0]
    ys = []
    for task in tasks:
        task_col = TASK_IDENTIFIER[task]
        ys.append(data[1][task_col])
    y = np.array(ys).transpose()
    return x, y


def random_sample_controller(skip_target, model_space):
    """Random sample architectures from a model space, with sampling residual connection at a
    specified skip_target in [0,1]

    Parameters
    ----------
    skip_target : float
        the probability of network-wise skip connections
    model_space : amber.architect.ModelSpace
        an instance of model search space, where computation operations are sampled

    Returns
    ----------
    list of int
        a sampled list of architecture tokens
    """
    arc_seq = []
    num_layers = len(model_space)
    num_skips = sum([layer_id for layer_id in range(num_layers)])
    skips = np.zeros(num_skips, dtype=int)
    skips_idx = np.random.choice(np.arange(num_skips), size=int(num_skips * skip_target), replace=False)
    skips[skips_idx] = 1
    skips = skips.tolist()
    print(np.mean(skips))
    for layer_id in range(num_layers):
        count = np.random.choice(np.arange(len(model_space[layer_id])), size=1, replace=False).tolist()
        if layer_id > 0:
            skip = [skips.pop(0) for _ in range(layer_id)]
            count.extend(skip)
        arc_seq.extend(count)
    return arc_seq


def read_controller_train_history(fn, last_only=None):
    """

    Parameters
    ----------
    fn : str
        filepath to "train_hist.tsv" generated by amber controller
    last_only : int, or None
        only read specified number of episodes from the last; if none,
        look for the best over controller train history

    Returns
    -------
    best_arc : list
        architecture tokens with the best reward
    best_auc : float
        the reward (i.e. auc) for the best architecture
    """
    d = pd.read_csv(fn, sep=",")
    if last_only is not None:
        d = d.loc[d.iloc[:, 0] >= (d.shape[0] - last_only)]
    best_auc = - np.inf
    best_arc = None
    for i in range(d.shape[0]):
        index, auc = d.iloc[i][0], d.iloc[i][2]
        arc = d.iloc[i][3:].to_list()
        if auc > best_auc:
            best_arc = arc
            best_auc = auc
    return best_arc, best_auc


def robust_spearmanr(y_true, y_score):
    """Reward function, Spearman correlation with resistence on all-constant predictions
    and avoids zero-divisions

    Parameters
    ----------
    y_true : np.array
        observed values
    y_score : np.array
        predicted values

    Returns
    -------
    rsp : float
        1 + spearman correlation coefficient
    """
    if len(set(y_score)) < 10:
        sp = 0
    else:
        sp = ss.spearmanr(y_true.flatten(), y_score.flatten()).correlation
    rsp = sp + 1  # avoid division by 0
    return rsp


def get_model_space(out_filters=64, num_layers=9, num_pool=4):
    """Model search space getter

    Parameters
    ----------
    out_filters : int
        number of filters for the base layer (the bottom layer); will multiple by a factor of 2 whenver a strided
        pooling is encounterd, as specified in num_pool
    num_layers : int
        number of layers in the model space
    num_pool : int
        number of strided pooling

    Returns
    -------
    model_space : amber.architect.ModelSpace
    """
    model_space = ModelSpace()
    expand_layers = [num_layers // num_pool * i - 1 for i in range(num_pool)]
    for i in range(num_layers):
        model_space.add_layer(i, [
            Operation('conv1d', filters=out_filters, kernel_size=8, activation='relu'),
            Operation('conv1d', filters=out_filters, kernel_size=4, activation='relu'),
            Operation('conv1d', filters=out_filters, kernel_size=8, activation='relu', dilation=4),
            Operation('conv1d', filters=out_filters, kernel_size=4, activation='relu', dilation=4),
            Operation('maxpool1d', filters=out_filters, pool_size=4, strides=1),
            Operation('avgpool1d', filters=out_filters, pool_size=4, strides=1),
            Operation('identity', filters=out_filters),
        ])
        if i in expand_layers:
            out_filters *= 2
    return model_space


def main(mode, wd, dataset, tasks, aux_reward_weight=1.0, enable_run=True):
    """Main wrapper for amber-croton creation

    Parameters
    ----------
    mode : str
        must be in ['search', 'train', 'random'].
    wd : str
        working directory
    dataset : str
        dataset identifier. only funcional for 'FORECasT'
    tasks : list of str
        list of task identifiers
    aux_reward_weight : float
        Auxilary reward for a second validation data. *Currently Non-functional*.
    enable_run : bool
        Only used in search mode. If true, will run amber search upon called; otherwise, return the contructed
        `amber.Amber` instance

    Returns
    -------
    None
    """
    if dataset == 'forecast':
        train_data = load_pickle_data(tasks=tasks, pickle_fp='./data/data/Forecast/train_.pkl')
    else:
        raise ValueError("Unknown dataset identifier: %s" % dataset)
    val_data = load_pickle_data(tasks=tasks, pickle_fp='./data/data/Forecast/val_.pkl')
    test_data = load_pickle_data(tasks=tasks, pickle_fp='./data/data/Forecast/test_.pkl')

    # First, define the components we need to use
    type_dict = {
        'controller_type': 'GeneralController',
        'modeler_type': 'EnasCnnModelBuilder',
        'knowledge_fn_type': 'zero',
        'reward_fn_type': 'LossAucReward',
        'manager_type': 'EnasManager',
        'env_type': 'EnasTrainEnv'
    }

    child_batchsize = 512
    # Next, define the specifics
    input_node = Operation('input', shape=(60, 4), name="input")
    output_node = Operation('dense', units=len(tasks), activation='sigmoid')
    model_compile_dict = {
        'loss': 'binary_crossentropy',
        'optimizer': 'adam'
    }
    if mode == 'search':
        model_space = get_model_space(out_filters=32, num_layers=8, num_pool=1)
        pickle.dump(model_space, open(os.path.join(wd, "model_space.pkl"), "wb"))
    else:
        model_space = pickle.load(open(os.path.join(wd, "model_space.pkl"), "rb"))

    fc_units = 32
    flatten_op = 'GAP'
    use_ppo_loss = False
    samps_per_controller_step = 100
    specs = {
        'model_space': model_space,

        'controller': {
            'share_embedding': {i: 0 for i in range(1, len(model_space))},
            'with_skip_connection': True,
            'skip_connection_unique_connection': False,
            'skip_weight': 0.3,
            'skip_target': 0.4,
            'lstm_size': 64,
            'lstm_num_layers': 1,
            'kl_threshold': 0.01,
            'train_pi_iter': 100 if use_ppo_loss else 10,
            'optim_algo': 'adam',
            'rescale_reward': False,
            'use_ppo_loss': use_ppo_loss,
            'temperature': 2.,
            'lr_init': 0.001,
            'tanh_constant': 1.5,
            'buffer_size': 1,
            'batch_size': 10
        },

        'model_builder': {
            'dag_func': 'EnasConv1dDAG',
            'batch_size': child_batchsize,
            'inputs_op': [input_node],
            'outputs_op': [output_node],
            'model_compile_dict': model_compile_dict,
            'dag_kwargs': {
                'stem_config': {
                    'flatten_op': flatten_op.lower(),
                    'fc_units': fc_units
                }
            }
        },

        'knowledge_fn': {
            'data': {},
            'params': {}
        },

        'reward_fn': {
            'method': robust_spearmanr,
        },

        'manager': {
            'data': {
                'train_data': train_data,
                'validation_data': val_data
            },
            'params': {
                'epochs': 2,
                'child_batchsize': child_batchsize,
                'store_fn': 'general',
                'working_dir': wd,
                'verbose': 2
            }
        },

        'train_env': {
            'max_episode': 120,
            'max_step_per_ep': samps_per_controller_step,
            'working_dir': wd,
            'time_budget': "24:00:00",
            'with_input_blocks': False,
            'with_skip_connection': True
        }
    }

    # finally, run program
    if mode == 'search':
        amb = Amber(types=type_dict, specs=specs)
        if enable_run is True:
            amb.run()
        return amb
    else:
        wsf = 6
        model_compile_dict['optimizer'] = 'adam'
        kmb = KerasResidualCnnBuilder(
            inputs_op=input_node,
            output_op=output_node,
            fc_units=fc_units * wsf,
            flatten_mode=flatten_op,
            model_compile_dict=model_compile_dict,
            model_space=model_space,
            dropout_rate=0.4,
            wsf=wsf,
            add_conv1_under_pool=True
        )
        if mode == 'train':
            best_arc, best_auc = read_controller_train_history(fn=os.path.join(wd, 'train_history.csv'),
                                                               last_only=samps_per_controller_step)
            print("best_arc=%s" % best_arc)
            print("best auc=%s" % best_auc)
            wd_ = os.path.join(wd, 'train')
        elif mode == 'random':
            # 25/28 is learned from amber
            best_arc = random_sample_controller(skip_target=25 / 28, model_space=model_space)
            wd_ = os.path.join(wd, 'random')
        verbose = 2
        os.makedirs(wd_, exist_ok=True)
        model = kmb(model_states=best_arc)
        plot_model(model, to_file=os.path.join(wd_, "model.png"))
        model_weight_fp = os.path.join(wd_, "bestmodel.h5")
        checkpointer = ModelCheckpoint(
            model_weight_fp,
            monitor='val_loss',
            save_best_only=True,
            save_weights_only=False,
            verbose=verbose
        )
        earlystopper = EarlyStopping(
            monitor='val_loss',
            patience=50,
            verbose=verbose
        )

        hist = model.fit(
            train_data[0], train_data[1],
            epochs=500,
            batch_size=child_batchsize,
            verbose=verbose,
            validation_data=val_data,
            callbacks=[checkpointer, earlystopper]
        )

        model.load_weights(model_weight_fp)
        val_pred = model.predict(val_data[0])
        val_df = {'obs_%s' % k: val_data[1][:, i] for i, k in enumerate(tasks)}
        val_df.update({
            'pred_%s' % k: val_pred[:, i] for i, k in enumerate(tasks)})
        val_df = pd.DataFrame(val_df)
        val_df.to_csv(os.path.join(wd_, "val.tsv"), sep="\t", index=False)
        fo = open(os.path.join(wd_, "metrics.txt"), "w")
        for task in tasks:
            print("%s pearson=%.5f" % (task, ss.pearsonr(val_df['obs_%s' % task], val_df['pred_%s' % task])[0]))
            fo.write(
                "VAL\t%s pearson\t%.5f\n" % (task, ss.pearsonr(val_df['obs_%s' % task], val_df['pred_%s' % task])[0]))
        test_pred = model.predict(test_data[0])
        test_df = {'obs_%s' % k: test_data[1][:, i] for i, k in enumerate(tasks)}
        test_df.update({
            'pred_%s' % k: test_pred[:, i] for i, k in enumerate(tasks)})
        test_df = pd.DataFrame(test_df)
        test_df.to_csv(os.path.join(wd_, "test.tsv"), sep="\t", index=False)
        print('TEST')
        for task in tasks:
            print("%s pearson=%.5f" % (task, ss.pearsonr(test_df['obs_%s' % task], test_df['pred_%s' % task])[0]))
            fo.write("TEST\t%s pearson\t%.5f\n" % (
                task, ss.pearsonr(test_df['obs_%s' % task], test_df['pred_%s' % task])[0]))
        plot_training_history(hist, wd_)
        fo.close()


if __name__ == "__main__":
    if not run_from_ipython():
        parser = argparse.ArgumentParser(description='AMBER search for CROTON')
        parser.add_argument("--wd", type=str, help="working dir")
        parser.add_argument("--mode", type=str, choices=['search', 'train', 'random'], required=True, help="run mode")
        parser.add_argument("--dataset", type=str, choices=['all', 'forecast', 'sprout'], help="training dataset")
        parser.add_argument("--tasks", type=str, nargs='+',
                            choices=['del_freq', '1ins_freq', '1del_freq', 'frameshift_freq',
                                     'onemod3_freq', 'twomod3_freq',
                                     'avgdel_len', 'avgins_len', 'entropy'],
                            help="tasks of interest")
        parser.add_argument("--aux-weight", type=float, default=1.0,
                            help="weight for auxilary reward; must be in [0, 10]")

        args = parser.parse_args()
        os.makedirs(args.wd, exist_ok=True)
        with open(os.path.join(args.wd, "args.txt"), "w") as f:
            f.write("\n".join(sys.argv))
        main(mode=args.mode, wd=args.wd, dataset=args.dataset, tasks=args.tasks, aux_reward_weight=args.aux_weight,
             enable_run=True)