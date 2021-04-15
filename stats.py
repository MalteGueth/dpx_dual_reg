# -*- coding: utf-8 -*-
"""Utility functions for computing statistical measures.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import numpy as np

from scipy import stats

from mne import grand_average


def within_subject_cis(insts, ci=0.95):
    # see Morey (2008): Confidence Intervals from Normalized Data:
    # A correction to Cousineau (2005)

    # the type of the provided instances should be the same
    if not all(isinstance(inst, (dict, list)) for inst in insts):
        raise ValueError('instances must be of same type (either dict ot list)')

    # the number of subjects should be the same
    n_subj = np.unique([len(i) for i in insts])
    if len(n_subj) > 1:
        raise ValueError('inst must be of same length')

    if isinstance(insts[0], dict):
        subjs = insts[0].keys()
    else:
        subjs = np.arange(0, len(insts[0]))

    # correction factor for number of conditions
    n_cond = len(insts)
    corr_factor = np.sqrt(n_cond / (n_cond - 1))

    # compute overall subject ERPs
    subject_erp = {subj: grand_average([insts[cond][subj]
                                        for cond in range(n_cond)])
                   for subj in subjs}

    # compute condition grand averages
    grand_averages = [grand_average(list(cond.values()))
                      for cond in insts]

    # place holder for results
    n_channels = grand_averages[0].data.shape[0]
    n_times = grand_averages[0].data.shape[1]
    norm_erps = np.zeros((n_cond, int(n_subj), n_channels, n_times))

    # compute normed ERPs,
    # ((condition ERP - subject ERP) + grand average) * corr_factor
    for n_s, subj in enumerate(subjs):

        for ic, cond in enumerate(insts):
            erp_data = cond[subj].data.copy() - subject_erp[subj].data
            erp_data = (erp_data + grand_averages[ic].data) * corr_factor

            norm_erps[ic, n_s, :] = erp_data

    erp_sem = np.zeros((n_cond, n_channels, n_times))
    for n_c in range(n_cond):
        erp_sem[n_c, :] = stats.sem(norm_erps[n_c, :], axis=0)

    confidence = np.zeros((n_cond, n_channels, n_times))

    for n_c in range(n_cond):
        confidence[n_c, :] = erp_sem[n_c, :] * \
                             stats.t.ppf((1 + ci) / 2.0, int(n_subj) - 1)

    return confidence
