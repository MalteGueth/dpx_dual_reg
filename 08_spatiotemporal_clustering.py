"""
=============================================================
Compute test statistics for effect of condition and moderator
=============================================================

Compute T- and F- for the effect of conditions and search for
significant spatio-temporal clusters of activity. Further,
estimate the moderating effect of behavioral performance measures
on a group-level.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import patsy

from sklearn.linear_model import LinearRegression

from mne.stats.cluster_level import _setup_adjacency, _find_clusters
from mne.decoding import get_coef
from mne.channels import find_ch_adjacency
from mne import read_epochs

# All parameters are defined in config.py
from config import subjects, fname, LoggingFormat
subjects = subjects[subjects != 5]

# experimental group for each subject
group = list()
for subject in subjects:
    # subject's experimental group
    if subject in [1, 4, 7, 8, 10, 13, 18, 19, 20, 21, 22, 23, 25]:
        group.append('Low')
    elif subject in [2, 3, 6, 9, 11, 12, 14, 15, 16, 17, 24, 26, 27]:
        group.append('High')
    elif subject in [5]:
        continue
group = pd.DataFrame({'group': group})

# load individual beta coefficients
betas_cue = np.load(fname.results + '/subj_betas_cue_by_block_m250_robust.npy')
generic_subj = subjects[0]

###############################################################################
# 1) import epochs to use as template
input_file = fname.output(subject=generic_subj,
                          processing_step='cue_epochs',
                          file_type='epo.fif')
cue_epo = read_epochs(input_file, preload=True)
cue_epo = cue_epo['Correct A', 'Correct B'].copy()
cue_epo = cue_epo.crop(tmin=-0.25, tmax=2.45, include_tmax=False)

# save the generic info structure of cue epochs (i.e., channel names, number of
# channels, etc.).
epochs_info = cue_epo.info
n_channels = len(epochs_info['ch_names'])
n_times = len(cue_epo.times)
times = cue_epo.times

###############################################################################
# 2) compute bootstrap confidence interval for cue betas and t-values

# create design matrix
group = patsy.dmatrix("group", group,
                      return_type='dataframe')

# set random state for replication
random_state = 42
random = np.random.RandomState(random_state)

# number of random samples
boot = 10000

# place holders for bootstrap samples
cluster_H0 = np.zeros(boot)
f_H0 = np.zeros(boot)

# setup adjacency matrix
n_tests = betas_cue.shape[2]
adjacency, ch_names = find_ch_adjacency(epochs_info, ch_type='eeg')
adjacency = _setup_adjacency(adjacency, n_tests, n_times)

# threshold parameters for clustering
threshold = dict(start=0.2, step=0.2)

# store a_bias (bootstrap) betas
group_pp_betas = np.zeros((boot, n_channels * n_times))

# center betas around zero
betas_null = betas_cue[0, ...] - betas_cue[0, ...].mean(axis=0)

# run bootstrap for regression coefficients
for i in range(boot):

    # log progress
    print(LoggingFormat.PURPLE +
          LoggingFormat.BOLD +
          'Running bootstrap sample %s of %s' % (i, boot) +
          LoggingFormat.END)

    # *** 2.1) create bootstrap sample ***
    # extract random subjects from overall sample
    resampled_subjects = random.choice(range(betas_cue[0, ...].shape[0]),
                                       betas_cue[0, ...].shape[0],
                                       replace=True)
    # resampled betas
    resampled_betas = betas_cue[0, ...][resampled_subjects, :]

    # *** 2.2) estimate effect of moderator (i.e., PBI) on group-level ***
    # set up and fit moderator (i.e., PBI) model using bootstrap sample
    model_boot = LinearRegression(fit_intercept=False)
    model_boot.fit(X=group.iloc[resampled_subjects], y=resampled_betas)

    # extract regression coefficients
    group_coefs = get_coef(model_boot, 'coef_')

    # save bootstrap betas
    for pred_i, predictor in enumerate(group.columns):
        if 'pbi_rt' in predictor:
            # store regression coefficient for moderator (i.e., PBI)
            group_pp_betas[i, :] = group_coefs[:, pred_i]

    # remove prev object
    del resampled_betas

    # *** 2.3) compute test statistic for bootstrap sample ***
    # compute standard error
    resampled_betas = betas_null[resampled_subjects, :]
    se = resampled_betas.std(axis=0) / np.sqrt(resampled_betas.shape[0])

    # compute t-values
    t_vals = resampled_betas.mean(axis=0) / se
    # transform to f-values
    f_vals = t_vals ** 2
    # save max f-value
    f_H0[i] = f_vals.max()

    # transpose for clustering
    t_vals = t_vals.reshape((n_channels, n_times))
    t_vals = np.transpose(t_vals, (1, 0))
    t_vals = t_vals.ravel()

    # compute clustering on squared t-values (i.e., f-values)
    clusters, cluster_stats = _find_clusters(t_vals,
                                             t_power=1,
                                             threshold=threshold,
                                             adjacency=adjacency,
                                             tail=0)

    # Save max cluster mass. Combined, the max cluster mass values
    # computed on the basis of the bootstrap samples provide an approximation
    # of the cluster mass distribution under H0
    if len(clusters):
        cluster_H0[i] = cluster_stats.max()
    else:
        cluster_H0[i] = 0.0

##############################################################################
# 3) Save results of bootstrap procedure

# save f-max distribution
np.save(fname.results + '/f_H0_10000b_2t_m250_null_robust_cue.npy', f_H0)
# save cluster mass distribution
np.save(fname.results + '/cluster_H0_10000b_2t_m250_null_robust_cue.npy',
        cluster_H0)
# save pbi_rt betas
np.save(fname.results + '/group_betas_m250_null_robust.npy', group_pp_betas)

# plot results
y_max = np.histogram(f_H0, bins=100)[0]
fig, ax = plt.subplots()
ax.hist(f_H0, ec='k', bins=100)
ax.axvline(np.quantile(f_H0, 0.95),
           color='red', linestyle='dashed', linewidth=2)
ax.text(x=np.quantile(f_H0, 0.95),
        y=np.max(y_max),
        s='p<0.05',
        fontdict=dict(fontsize=12),
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=1.0))
plt.savefig(fname.figures + '/F_max_distribution_cue_robust.pdf', dpi=300)
