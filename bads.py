"""
=================
Find bad channels
=================

Methods for finding bad (e.g., noisy) channels in EEG data.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import warnings

import numpy as np
from scipy.stats import median_abs_deviation as mad

from mne.io.base import BaseRaw


# main function which implements different methods
def find_bad_channels(inst, picks='eeg',
                      method='correlation',
                      mad_threshold=1,
                      std_threshold=1,
                      r_threshold=0.4,
                      percent_threshold=0.1,
                      time_step=1.0,
                      sfreq=None,
                      return_z_scores=False,
                      channels=None):

    # arguments to be passed to pick_types
    kwargs = {pick: True for pick in [picks]}

    # check that tha input data can be handled by the function
    if isinstance(inst, BaseRaw):
        # only keep data from desired channels
        inst = inst.copy().pick_types(**kwargs)
        dat = inst.get_data() * 1e6  # to microvolt
        channels = inst.ch_names
        sfreq = inst.info['sfreq']
    elif isinstance(inst, np.ndarray):
        if not channels:
            raise ValueError('If "inst" is not an instance of BaseRaw a list '
                             'of channel names must be provided')
        dat = inst
    else:
        raise ValueError('inst must be an instance of BaseRaw or a numpy array')

    # save shape of data
    n_channels, n_samples = dat.shape
    if n_channels != len(channels):
        raise ValueError("Number and channels and data dimensions don't match")

    # make sure method arguments are in a list
    if not isinstance(method, list):
        method = [method]

    # place holder for results
    bad_channels = dict()

    # 1) find channels with zero or near zero activity
    if 'flat' in method:
        # compute estimates of channel activity
        mad_flats = mad(dat, scale=1, axis=1) < mad_threshold
        std_flats = np.std(dat, axis=1) < std_threshold

        # flat channels identified
        flats = np.argwhere(np.logical_or(mad_flats, std_flats))
        flats = np.asarray([channels[int(flat)] for flat in flats])

        # warn user if too many channels were identified as flat
        if len(flats) > (n_channels / 2):
            warnings.warn('Too many channels have been identified as "flat"! '
                          'Make sure the input values in "inst" are provided '
                          'on a volt scale. '
                          'Otherwise try choosing another (meaningful) '
                          'threshold for identification.')

        bad_channels.update(flat=flats)

    # 3) find bad channels by deviation (high variability in amplitude)
    if 'deviation' in method:

        # mean absolute deviation (MAD) scores for each channel
        mad_scores = \
            [mad(dat[i, :], scale=1) for i in range(n_channels)]

        # compute robust z-scores for each channel
        rz_scores = \
            0.6745 * (mad_scores - np.nanmedian(mad_scores)) / mad(mad_scores,
                                                                   scale=1)

        # channels identified by deviation criterion
        bad_deviation = \
            [channels[i] for i in np.where(np.abs(rz_scores) > 5.0)[0]]

        bad_channels.update(deviation=np.asarray(bad_deviation))

        if return_z_scores:
            bad_channels.update(deviation_z_scores=rz_scores)

    # 3) find channels with low correlation to other channels
    if 'correlation' in method:

        # check that sampling frequency argument was provided
        if not sfreq:
            raise ValueError('If "inst" is not an instance of BaseRaw a '
                             'sampling frequency must be provided. Usually '
                             'the sampling frequency of the EEG recording in'
                             'question.')

        # based on the length of the provided data,
        # determine size and amount of time windows for analyses
        corr_frames = time_step * sfreq
        corr_window = np.arange(corr_frames)

        # sample index (i.e., time offsets) for each window to time window
        # to use for correlation analyis
        corr_offsets = np.arange(1, (n_samples - corr_frames), corr_frames)
        n_corr_steps = corr_offsets.shape[0]
        # place holders for correlation coefficients
        max_r = np.ones((n_channels, n_corr_steps))
        channel_r = np.ones((n_corr_steps, n_channels))

        # create time windows for analysis
        dat_t = np.transpose(dat)
        dat_windowed = np.reshape(
            np.transpose(dat_t[0: corr_window.shape[0] * n_corr_steps, :]),
            (n_channels, corr_window.shape[0], n_corr_steps),
            order="F",)

        # compute (pearson) correlation coefficient across channels
        # (for each channel and analysis time window)
        # take the absolute of the 98th percentile of the correlations with
        # the other channels as a measure of how well that channel is correlated
        # to other channels
        for k in range(0, n_corr_steps):
            eeg_portion = np.transpose(np.squeeze(dat_windowed[:, :, k]))
            window_correlation = np.corrcoef(np.transpose(eeg_portion))
            abs_corr = \
                np.abs(
                    np.subtract(
                        window_correlation, np.diag(np.diag(window_correlation))
                    )
                )
            channel_r[k, :] = np.quantile(abs_corr, 0.98, axis=0)

        # fill in the actual correlations
        max_r[np.arange(0, n_channels), :] = np.transpose(channel_r)


        # check which channels correlate badly with the other channels (i.e.,
        # are below correlation threshold) in a certain fraction of windows
        # (bad_time_threshold)
        thresholded_correlations = max_r < r_threshold
        thresholded_correlations = thresholded_correlations.astype(int)
        frac_bad_corr_windows = np.mean(thresholded_correlations, axis=1)

        # find the corresponding channel names and return
        bad_idxs = np.argwhere(frac_bad_corr_windows > percent_threshold)
        uncorrelated_channels = [channels[int(bad)] for bad in bad_idxs]

        bad_channels.update(correlation=np.asarray(uncorrelated_channels))  # noqa: E501

    return bad_channels
