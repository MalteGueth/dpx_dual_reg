"""
=====================================
Compute T-map for effect of condition
=====================================

Mass-univariate analysis of cue evoked activity.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mne.stats.cluster_level import _setup_adjacency, _find_clusters
from mne.channels import make_1020_channel_selections, find_ch_adjacency
from mne.evoked import EvokedArray
from mne import read_epochs, grand_average

from config import subjects, fname
# exclude subjects 51
subjects = subjects[subjects != 51]

###############################################################################
# 1) Load results of bootstrap procedure

# load f-max distribution
f_H0 = np.load(fname.results + '/f_H0_10000b_2t_m250_null_robust.npy')
# load cluster mass distribution
cluster_H0 = np.load(fname.results +
                     '/cluster_H0_10000b_2t_m250_null_robust.npy')

# also load individual beta coefficients
betas = np.load(fname.results + '/subj_betas_cue_m250_robust.npy')
r2 = np.load(fname.results + '/subj_r2_cue_m250_robust.npy')

###############################################################################
# 1) import epochs to use as template

# baseline to be applied
baseline = (-0.300, -0.050)

# import the output from previous processing step
input_file = fname.output(subject=subjects[0],
                          processing_step='cue_epochs',
                          file_type='epo.fif')
cue_epo = read_epochs(input_file, preload=True)
cue_epo = cue_epo['Correct A', 'Correct B'].copy()
cue_epo_nb = cue_epo.copy().crop(tmin=-0.250, tmax=2.450, include_tmax=False)
cue_epo = cue_epo.apply_baseline(baseline).crop(tmin=-0.300)

# save the generic info structure of cue epochs (i.e., channel names, number of
# channels, etc.).
epochs_info = cue_epo_nb.info
n_channels = len(epochs_info['ch_names'])
n_times = len(cue_epo_nb.times)
times = cue_epo_nb.times
tmin = cue_epo_nb.tmin

# split channels into ROIs for results section
selections = make_1020_channel_selections(epochs_info, midline='12z')

# placeholder for results
betas_evoked = dict()
r2_evoked = dict()

# ###############################################################################
# 2) loop through subjects and extract betas
for n_subj, subj in enumerate(subjects):
    subj_beta = betas[n_subj, :]
    subj_beta = subj_beta.reshape((n_channels, n_times))
    betas_evoked[str(subj)] = EvokedArray(subj_beta, epochs_info, tmin)

    subj_r2 = r2[n_subj, :]
    subj_r2 = subj_r2.reshape((n_channels, n_times))
    r2_evoked[str(subj)] = EvokedArray(subj_r2, epochs_info, tmin)

effect_of_cue = grand_average([betas_evoked[str(subj)] for subj in subjects])
cue_r2 = grand_average([r2_evoked[str(subj)] for subj in subjects])

###############################################################################
# 3) Plot beta weights for the effect of condition

# arguments fot the time-series maps
ts_args = dict(gfp=False,
               time_unit='s',
               ylim=dict(eeg=[-6.5, 6.5]),
               xlim=[-.25, 2.5])

# times to plot
ttp = [0.20, 0.35, 0.60, 0.70, 1.00, 1.25, 2.35]

# arguments fot the topographical maps
topomap_args = dict(sensors=False,
                    time_unit='ms',
                    vmin=7, vmax=-7,
                    average=0.05,
                    extrapolate='head')

# create plot
title = 'Regression coefficients (B - A, 64 EEG channels)'
fig = effect_of_cue.plot_joint(ttp,
                               ts_args=ts_args,
                               topomap_args=topomap_args,
                               title=title,
                               show=False)
fig.axes[-1].texts[0]._fontproperties._size = 12.0  # noqa
fig.axes[-1].texts[0]._fontproperties._weight = 'bold'  # noqa
fig.axes[0].set_xticks(list(np.arange(-0.25, 2.55, 0.25)), minor=False)
fig.axes[0].set_yticks(list(np.arange(-6.0, 6.5, 3.0)), minor=False)
fig.axes[0].set_xticklabels(list(np.arange(-250, 2550, 250)))
fig.axes[0].set_xlabel('Time (ms)')
fig.axes[0].axhline(y=0.0, xmin=-0.5, xmax=2.5,
                    color='black', linestyle='dashed', linewidth=0.8)
fig.axes[0].axvline(x=0.0, ymin=-6.0, ymax=6.0,
                    color='black', linestyle='dashed', linewidth=0.8)
fig.axes[0].spines['top'].set_visible(False)
fig.axes[0].spines['right'].set_visible(False)
fig.axes[0].spines['left'].set_bounds(-6.0, 6.0)
fig.axes[0].spines['bottom'].set_bounds(-0.25, 2.5)
fig.axes[0].xaxis.set_label_coords(0.5, -0.2)
w, h = fig.get_size_inches()
fig.set_size_inches(w * 1.15, h * 1.15)
fig_name = fname.figures + '/Evoked_average_betas.pdf'
fig.savefig(fig_name, dpi=300)

###############################################################################
# 4) Plot R-squared for the effect of condition

# arguments fot the time-series maps
ts_args = dict(gfp=False,
               time_unit='s',
               unit=False,
               ylim=dict(eeg=[-0.005, 0.06]),
               xlim=[-0.25, 2.5])

# times to plot
ttp = [0.20, 0.35, 0.60, 0.70, 1.00, 1.25, 2.35]

# arguments fot the topographical maps
topomap_args = dict(cmap='magma_r',
                    scalings=dict(eeg=1),
                    sensors=False,
                    time_unit='ms',
                    vmin=0.0, vmax=0.05,
                    average=0.05,
                    extrapolate='head')

# create R-squared plot
title = 'Proportion of variance explained by cue type'
fig = cue_r2.plot_joint(ttp,
                        ts_args=ts_args,
                        topomap_args=topomap_args,
                        title=title,
                        show=False)
fig.axes[-1].texts[0]._fontproperties._size = 12.0  # noqa
fig.axes[-1].texts[0]._fontproperties._weight = 'bold'  # noqa
fig.axes[0].set_xticks(list(np.arange(-0.25, 2.55, 0.25)), minor=False)
fig.axes[0].set_xticklabels(list(np.arange(-250, 2550, 250)))
fig.axes[0].set_xlabel('Time (ms)')
fig.axes[0].set_yticks(list(np.arange(0.0, 0.065, 0.03)), minor=False)
fig.axes[0].axvline(x=0.0, ymin=0.0, ymax=1,
                    color='black', linestyle='dashed', linewidth=.8)
fig.axes[0].spines['top'].set_visible(False)
fig.axes[0].spines['right'].set_visible(False)
fig.axes[0].spines['left'].set_bounds(0.0, 0.06)
fig.axes[0].spines['bottom'].set_bounds(-0.25, 2.5)
fig.axes[0].xaxis.set_label_coords(0.5, -0.2)
fig.axes[0].set_ylabel('Average R-squared')
w, h = fig.get_size_inches()
fig.set_size_inches(w * 1.15, h * 1.15)
fig_name = fname.figures + '/Evoked_average_R2.pdf'
fig.savefig(fig_name, dpi=300)

###############################################################################
# 5) Estimate t-test based on original condition betas
se = betas.std(axis=0) / np.sqrt(betas.shape[0])
t_vals = betas.mean(axis=0) / se
f_vals = t_vals ** 2

# transpose for later clustering
t_clust = t_vals.reshape((n_channels, n_times))
f_clust = np.transpose(t_clust, (1, 0))
t_clust = f_clust.ravel()

# get upper CI bound from cluster mass H0
# (equivalent to alpha < 0.01 sig. level)

# f values above alpha level (based on f-max statistics)
sig_mask = f_vals > np.quantile(f_H0, [.99], axis=0)
# clusters threshold
cluster_thresh = np.quantile(cluster_H0, [0.99], axis=0)

# ###############################################################################
# 6) Plot results

# back projection to channels x time points
t_vals = t_vals.reshape((n_channels, n_times))
f_vals = f_vals.reshape((n_channels, n_times))
sig_mask = sig_mask.reshape((n_channels, n_times))

# create evoked object containing the resulting t-values
group_t = dict()
group_t['effect of cue (B-A)'] = EvokedArray(t_vals, epochs_info, tmin)

# initialise plot
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(23, 5.5))

# plot channel ROIs
for s, selection in enumerate(selections):
    picks = selections[selection]

    group_t['effect of cue (B-A)'].plot_image(xlim=[-0.25, 2.5],
                                              picks=picks,
                                              clim=dict(eeg=[-12, 12]),
                                              colorbar=False,
                                              axes=ax[s],
                                              mask=sig_mask,
                                              mask_cmap='RdBu_r',
                                              mask_alpha=0.5,
                                              show=False,
                                              unit=False,
                                              # keep values scale
                                              scalings=dict(eeg=1)
                                              )

    # tweak plot appearance
    if selection in {'Left', 'Right'}:
        title = selection + ' hemisphere'
    else:
        title = 'Midline'
    ax[s].title._text = title # noqa
    ax[s].set_ylabel('Channels', labelpad=10.0,
                     fontsize=11.0, fontweight='bold')

    ax[s].set_xlabel('Time (s)',
                     labelpad=10.0, fontsize=11.0, fontweight='bold')

    ax[s].set_xticks(list(np.arange(-.25, 2.55, .25)), minor=False)
    ax[s].set_xticklabels(list(np.arange(-250, 2550, 250)), rotation=45)
    ax[s].set_xlabel('Time (ms)')
    ax[s].set_yticks(np.arange(len(picks)), minor=False)
    labels = [group_t['effect of cue (B-A)'].ch_names[i] for i in picks]
    ax[s].set_yticklabels(labels, minor=False)
    ax[s].spines['top'].set_visible(False)
    ax[s].spines['right'].set_visible(False)
    ax[s].spines['left'].set_bounds(-0.5, len(picks)-0.5)
    ax[s].spines['bottom'].set_bounds(-.25, 2.5)
    ax[s].texts = []

    # add intercept line (at 0 s) and customise figure boundaries
    ax[s].axvline(x=0, ymin=0, ymax=len(picks),
                  color='black', linestyle='dashed', linewidth=1.0)

    colormap = cm.get_cmap('RdBu_r')
    orientation = 'vertical'
    norm = Normalize(vmin=-12.0, vmax=12.0)
    divider = make_axes_locatable(ax[s])
    cax = divider.append_axes('right', size='2.5%', pad=0.2)
    cbar = ColorbarBase(cax, cmap=cm.get_cmap('RdBu_r'),
                        ticks=[-12.0, -6.0, 0.0, 6.0, 12.0], norm=norm,
                        label=r'Effect of cue (T-value B-A)',
                        orientation=orientation)
    cbar.outline.set_visible(False)
    cbar.ax.set_frame_on(True)
    label = r'Difference B-A (in $\mu V$)'
    for key in ('left', 'top',
                'bottom' if orientation == 'vertical' else 'right'):
        cbar.ax.spines[key].set_visible(False)

    fig.subplots_adjust(
        left=0.05, right=0.95, bottom=0.15, wspace=0.3, hspace=0.25)

# save figure
fig.savefig(fname.figures + '/T-map_image_effect_of_cue.pdf', dpi=300)

# # inspect topomaps
# group_t['effect of cue (B-A)'].plot_topomap(times=[0.20, 0.50, 1.3],
#                                             average=0.1,
#                                             mask=sig_mask,
#                                             units=None,
#                                             scalings=dict(eeg=1),
#                                             outlines='head',
#                                             sensors=True)

# ###############################################################################
# 7) Plot results

# set up channel adjacency matrix
n_tests = betas.shape[1]
adjacency, ch_names = find_ch_adjacency(epochs_info, ch_type='eeg')
adjacency = _setup_adjacency(adjacency, n_tests, n_times)

# threshold parameters for clustering
threshold = dict(start=0.2, step=0.2)

clusters, cluster_stats = _find_clusters(t_clust,
                                         t_power=1,
                                         threshold=threshold,
                                         adjacency=adjacency,
                                         tail=0)

# get significant clusters
cl_sig_mask = cluster_stats > cluster_thresh

cl_sig_mask = np.transpose(
    cl_sig_mask.reshape((n_times, n_channels)), (1, 0))
cluster_stats = np.transpose(
    cluster_stats.reshape((n_times, n_channels)), (1, 0))

# create evoked object containing the resulting t-values
cluster_map = dict()

cluster_map['ST-clustering effect of cue (B-A)'] = EvokedArray(cluster_stats,
                                                               epochs_info,
                                                               tmin)

# initialise plot
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 5))

# plot channel ROIs
for s, selection in enumerate(selections):
    picks = selections[selection]

    effect_of_cue.plot_image(xlim=[-0.25, 2.5],
                             picks=picks,
                             clim=dict(eeg=[-5.0, 5.0]),
                             colorbar=False,
                             axes=ax[s],
                             mask=cl_sig_mask,
                             mask_cmap='RdBu_r',
                             mask_alpha=0.5,
                             show=False)

    # tweak plot appearance
    if selection in {'Left', 'Right'}:
        title = selection + ' hemisphere'
    else:
        title = 'Midline'
    ax[s].title._text = title # noqa
    ax[s].set_ylabel('Channels', labelpad=10.0,
                     fontsize=11.0, fontweight='bold')

    ax[s].set_xlabel('Time (s)',
                     labelpad=10.0, fontsize=11.0, fontweight='bold')

    ax[s].set_xticks(list(np.arange(-.25, 2.55, .25)), minor=False)
    ax[s].set_yticks(np.arange(len(picks)), minor=False)
    labels = [cluster_map['ST-clustering effect of cue (B-A)'].ch_names[i]
              for i in picks]
    ax[s].set_yticklabels(labels, minor=False)
    ax[s].spines['top'].set_visible(False)
    ax[s].spines['right'].set_visible(False)
    ax[s].spines['left'].set_bounds(-0.5, len(picks)-0.5)
    ax[s].spines['bottom'].set_bounds(-.25, 2.5)
    ax[s].texts = []

    # add intercept line (at 0 s) and customise figure boundaries
    ax[s].axvline(x=0, ymin=0, ymax=len(picks),
                  color='black', linestyle='dashed', linewidth=1.0)

    colormap = cm.get_cmap('RdBu_r')
    orientation = 'vertical'
    norm = Normalize(vmin=-5.0, vmax=5.0)
    divider = make_axes_locatable(ax[s])
    cax = divider.append_axes('right', size='2.5%', pad=0.2)
    cbar = ColorbarBase(cax, cmap=cm.get_cmap('RdBu_r'),
                        ticks=[-5.0, 0, 5.0], norm=norm,
                        label=r'Effect of cue (ß-weight B-A)',
                        orientation=orientation)
    cbar.outline.set_visible(False)
    cbar.ax.set_frame_on(True)
    label = r'Difference B-A (in $\mu V$)'
    for key in ('left', 'top',
                'bottom' if orientation == 'vertical' else 'right'):
        cbar.ax.spines[key].set_visible(False)

    fig.subplots_adjust(
        left=0.05, right=0.95, bottom=0.15, wspace=0.3, hspace=0.25)

# save figure
fig.savefig(fname.figures + '/Cluster-map_image_effect_of_cue.pdf', dpi=300)
