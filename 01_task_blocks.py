"""
=================================
Extract relevant segments of data
=================================

Segments that were recorded during that where recorded during task performance
are extracted and the self-paced breaks (in between experimental blocks) are
dropped.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import matplotlib.pyplot as plt

from mne.io import read_raw_fif
from mne import events_from_annotations, concatenate_raws, open_report

# All parameters are defined in config.py
from config import fname, n_jobs, parser, LoggingFormat

# Handle command line arguments
args = parser.parse_args()
subject = args.subject

print(LoggingFormat.PURPLE +
      LoggingFormat.BOLD +
      'Extracting task blocks for subject %s' % subject +
      LoggingFormat.END)

###############################################################################
# 1) Import the output from previous processing step
input_file = fname.output(subject=subject,
                          processing_step='raw_files',
                          file_type='raw.fif')
raw = read_raw_fif(input_file, preload=True)

# drop status channel
raw.drop_channels('Status')

###############################################################################
# 2) Find periods of time in the data with no presented stimuli (i.e., the
# self-paced breaks)

# relevant events
ids = {'127': 1,
       '245': 2,
       '13': 3,
       '12': 4,
       '113': 5,
       '112': 6,
       '70': 7,
       '71': 8,
       '72': 9,
       '73': 10,
       '74': 11,
       '75': 12,
       '76': 13,
       '77': 14,
       '78': 15,
       '79': 16,
       '80': 17,
       '81': 18
       }

# extract events
events = events_from_annotations(raw, event_id=ids)

# cue events
cue_evs = events[0]
cue_evs = cue_evs[(cue_evs[:, 2] >= 7) & (cue_evs[:, 2] <= 12)]

# latencies and difference between two consecutive cues
latencies = cue_evs[:, 0] / raw.info['sfreq']
diffs = [(y - x) for x, y in zip(latencies, latencies[1:])]

# get first event after a long break (i.e., when the time difference between
# stimuli is greater than 10 seconds). This should only be the case in between
# task blocks.
breaks = [diff for diff in range(len(diffs)) if diffs[diff] > 10]
print('\n Identified breaks at positions', breaks)

###############################################################################
# 3) Save start and end points of task blocks
# start of first block
b1s = latencies[breaks[0] + 1] - 2
# end of first block
b1e = latencies[breaks[1]] + 6

# start of second block
b2s = latencies[breaks[1] + 1] - 2
# end of second block
b2e = latencies[breaks[2]] + 6

# start of third block
b3s = latencies[breaks[2] + 1] - 2
# end of third block
b3e = latencies[breaks[3]] + 6

# start of fourth block
b4s = latencies[breaks[3] + 1] - 2
# end of fourth block
b4e = latencies[-1] + 6

###############################################################################
# 4) Extract data belonging to the task blocks an concatenate
# block 1
raw_bl1 = raw.copy().crop(tmin=b1s, tmax=b1e)
# block 2
raw_bl2 = raw.copy().crop(tmin=b2s, tmax=b2e)
# block 1
raw_bl3 = raw.copy().crop(tmin=b3s, tmax=b3e)
# block 2
raw_bl4 = raw.copy().crop(tmin=b4s, tmax=b4e)

# concatenate data
raw_bl = concatenate_raws([raw_bl1, raw_bl2, raw_bl3, raw_bl4])

###############################################################################
# 5) Remove slow drifts and line noise

# Setting up band-pass filter from 0.1 - 40 Hz
#
# FIR filter parameters
# ---------------------
# Designing a one-pass, zero-phase, non-causal bandpass filter:
# - Windowed time-domain design (firwin) method
# - Hamming window with 0.0194 passband ripple and 53 dB stopband attenuation
# - Lower passband edge: 0.10
# - Lower transition bandwidth: 0.10 Hz (-6 dB cutoff frequency: 0.05 Hz)
# - Upper passband edge: 40.00 Hz
# - Upper transition bandwidth: 10.00 Hz (-6 dB cutoff frequency: 45.00 Hz)
# - Filter length: 8449 samples (33.004 sec)
raw_bl_filt = raw_bl.copy().filter(l_freq=0.1, h_freq=40.,
                               picks=['eeg', 'eog'],
                               filter_length='auto',
                               l_trans_bandwidth='auto',
                               h_trans_bandwidth='auto',
                               method='fir',
                               phase='zero',
                               fir_window='hamming',
                               fir_design='firwin',
                               n_jobs=n_jobs)

# plot filtered data
filt_plot = raw_bl_filt.plot(scalings=dict(eeg=50e-6, eog=50e-6),
                             n_channels=len(raw_bl_filt.info['ch_names']),
                             show=False)

# plot power spectral density
fig, ax = plt.subplots(figsize=(10, 5))
raw_bl_filt.plot_psd(fmax=70, show=False, ax=ax)

###############################################################################
# 6) Export data to .fif for further processing
# output path
output_path = fname.output(processing_step='task_blocks',
                           subject=subject,
                           file_type='raw.fif')

# sample down and save file
raw_bl_filt.save(output_path, overwrite=True)

###############################################################################
# 7) Create HTML report
blocks_duration = '<p>Block 1 Duration from %s to %s seconds.<br>'\
                  'Block 1 length: %s seconds<p>'\
                  '<p>Block 2 Duration from %s to %s seconds.<br>'\
                  'Block 2 length: %s seconds<p>' \
                  '<p>Block 2 Duration from %s to %s seconds.<br>' \
                  'Block 3 length: %s seconds<p>' \
                  '<p>Block 2 Duration from %s to %s seconds.<br>' \
                  'Block 4 length: %s seconds<p>' \
                  % (round(b1s, 2), round(b1e, 2), round(b1e - b1s, 2),
                     round(b2s, 2), round(b2e, 2), round(b2e - b2s, 2),
                     round(b3s, 2), round(b3e, 2), round(b3e - b3s, 2),
                     round(b4s, 2), round(b4e, 2), round(b4e - b4s, 2))

with open_report(fname.report(subject=subject)[0]) as report:
    report.add_htmls_to_section(htmls=blocks_duration,
                                captions='Block durations',
                                section='Filtered data')
    report.add_figs_to_section(fig, 'Blocks PSD',
                               section='Filtered data',
                               replace=True)
    report.add_figs_to_section(filt_plot,
                               'Filtered data',
                               section='Filtered data',
                               replace=True)
    report.save(fname.report(subject=subject)[1], overwrite=True,
                open_browser=False)
