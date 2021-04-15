# -*- coding: utf-8 -*-
"""Utility functions for plotting.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sklearn.preprocessing import normalize


def plot_z_scores(z_scores, channels, bads=None, cmap='inferno', show=False):

    cmap = cm.get_cmap(cmap)

    # plot results
    z_colors = normalize(
        np.abs(z_scores).reshape((1, z_scores.shape[0]))).ravel()

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    fig, ax = plt.subplots(figsize=(20, 6))
    if z_scores.max() < 5.0:
        y_lim = 5
    else:
        y_lim = int(z_scores.max() + 2)

    for i in range(z_scores.shape[0]):
        ch = channels[i]
        # show channel names in red if bad by correlation
        if ch in bads:
            col = 'crimson'
        else:
            col = 'k'
        ax.axhline(y=5.0, xmin=-1.0, xmax=65,
                   color='crimson', linestyle='dashed', linewidth=2.0)
        ax.text(-5.0, 5.0, 'crit. Z-score', fontsize=14,
                verticalalignment='center', horizontalalignment='center',
                color='crimson', bbox=props)
        ax.bar(i, np.abs(z_scores[i]), width=0.9, color=cmap(z_colors[i]))
        ax.text(i, np.abs(z_scores[i]) + 0.25, ch, color=col,
                fontweight='bold', fontsize=9,
                ha='center', va='center', rotation=45)
    ax.set_ylim(0, y_lim)
    ax.set_xlim(-1, 64)

    plt.title('EEG channel deviation', {'fontsize': 15, 'fontweight': 'bold'})
    plt.xlabel('Channels', {'fontsize': 13}, labelpad=10)
    plt.ylabel('Abs. Z-Score', {'fontsize': 13}, labelpad=10)

    plt.xticks([])
    plt.yticks(fontsize=12)

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_bounds(0, y_lim)

    plt.close(fig)

    return fig.show() if show else fig
