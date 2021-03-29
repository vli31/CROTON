"""read in a collection of random-sample model architecture performances, and
plot the distribution agains amber-design
"""

# zzjfrank, jan 25 2021
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

target_dir = './outputs/forecast_freqs/'
tasks = ['1ins', '1del', 'del', 'onemod3', 'twomod3', 'frameshift']
tasks_to_label = {
    '1ins': '1 bp Insertion',
    '1del': '1 bp Deletion',
    'del': 'Deletion Freq.',
    'onemod3': '1 bp Frameshift Freq.',
    'twomod3': '2 bp Frameshift Freq.',
    'frameshift': 'Frameshift Freq.'
}


def read_bg():
    """Read a collection of random-sampled architecture performances as the background performance
    distribution in a model space

    The random-sampled architectures are created by setting `mode=random`. See the `main()` function
     in amber_cnn_sumstats.py for more information.

    Returns
    -------
    df : pandas.DataFrame
        a dataframe of testing performances
    """
    sub_dir = 'random_collections.denseResConn'
    cl = [os.path.join(target_dir, sub_dir, x, 'metrics.txt') for x in os.listdir(os.path.join(target_dir, sub_dir))]
    dfs = [pd.read_table(x, header=None, index_col=(0, 1), sep="\t") for x in cl if os.path.isfile(x)]
    df = pd.concat(dfs, axis=1)
    return df


df = read_bg()
df = df.iloc[6:]  # only testing
df.index = [x[1].split()[0].rstrip('_freq') for x in df.index]

amb = pd.read_table(os.path.join(target_dir, 'train', 'metrics.txt'), header=None, index_col=(0, 1), sep="\t")
amb = amb.iloc[6:]
amb.index = [x[1].split()[0].rstrip('_freq') for x in amb.index]

# performance over random
outperform_pct = {}
for idx in df.index:
    outperform_pct[idx] = np.mean(df.loc[idx] <= amb.loc[idx][2])
print(outperform_pct)

# plot
SIZE = 12
fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))
plot_df = pd.melt(df.transpose())
sns.boxplot(y='variable', x='value', order=tasks, data=plot_df,
            showfliers=False, width=.6, ax=ax)
for patch in ax.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .7))
ax.scatter(amb[2].loc[tasks], np.arange(0, amb.shape[0]), color='red', s=16, label='CROTON', zorder=10)  # *
sns.stripplot(x="value", y="variable", order=tasks, data=plot_df, color=".25",  # color=".3", marker='o',
              size=4, linewidth=0, label='Sampled CNNs')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[0:2], labels[0:2], fontsize=SIZE)
ax.xaxis.grid(True)
ax.set_xlim(0, 1)
ax.set_xlabel('Test Pearson Correlation', fontsize=SIZE)
ax.set_ylabel('')
ax.set_yticklabels([tasks_to_label[task] for task in tasks], fontsize=SIZE)
ax.set_xticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=SIZE - 2)

# for i, t in enumerate(tasks):
#    ax.text(x=1, y=i, s="%i %%" % (outperform_pct[t]*100), fontsize=12)
sns.despine(trim=True, left=True)
fig.tight_layout()
# fig.savefig('./images/ISMB_ECCB/amber/AMBER_vs_baselineCNN.DRN.pdf')
fig.savefig('./images/ISMB_ECCB/amber/AMBER_baselineCNN.png', dpi=350)
