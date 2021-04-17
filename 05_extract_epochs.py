"""
==================================
Extract epochs from continuous EEG
==================================

Extract epochs for each experimental condition

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
from os import path as op

import pandas as pd
import numpy as np

from mne import events_from_annotations, Epochs, open_report
from mne.io import read_raw_fif

# All parameters are defined in config.py
from config import fname, parser, LoggingFormat

# Handle command line arguments
args = parser.parse_args()
subject = args.subject

print(LoggingFormat.PURPLE +
      LoggingFormat.BOLD +
      'Extracting epochs for subject %s' % subject +
      LoggingFormat.END)

###############################################################################
# 1) Import the output from previous processing step
input_file = fname.output(subject=subject,
                          processing_step='repaired_with_ica',
                          file_type='raw.fif')
raw = read_raw_fif(input_file, preload=True)

# only keep EEG channels
raw.pick_types(eeg=True)

###############################################################################
# 2) Get events from continuous EEG data

# create a dictionary with event IDs for standardised handling
ev_id = {'112': 1,
         '113': 2,
         '12': 3,
         '13': 4,
         '70': 5,
         '71': 6,
         '72': 7,
         '73': 8,
         '74': 9,
         '75': 10,
         '76': 11,
         '77': 12,
         '78': 13,
         '79': 14,
         '80': 15,
         '81': 16,
         'EDGE boundary': 17}
# extract events
events = events_from_annotations(raw, event_id=ev_id, regexp=None)

###############################################################################
# 3) Recode events into respective conditions and add information about valid
# and invalid responses

# copy of events
new_evs = events[0].copy()

# global variables
trial = 0
broken = []
sfreq = raw.info['sfreq']
block_end = events[0][events[0][:, 2] == 17, 0] / sfreq
# place holders for results
block = []
probe_ids = []
reaction = []
rt = []

# loop trough events and recode them
for event in range(len(new_evs[:, 2])):

    # - if event is a cue stimulus -
    if new_evs[event, 2] in {5, 6, 7, 8, 9, 10}:

        # save block based on onset (before or after break)
        if (new_evs[event, 0] / sfreq) < block_end[0]:
            block.append(0)
        elif ((new_evs[event, 0] / sfreq) > block_end[0]) & ((new_evs[event, 0] / sfreq) < block_end[1]):
            block.append(0)
        elif ((new_evs[event, 0] / sfreq) > block_end[1]) & ((new_evs[event, 0] / sfreq) < block_end[2]):
            block.append(1)
        elif (new_evs[event, 0] / sfreq) > block_end[2]:
            block.append(1)

        # -- if next event is a false reaction --
        if new_evs[event + 1, 2] in {1, 2}:
            # if event is an A-cue
            if new_evs[event, 2] == 5:
                # recode as too soon A-cue
                new_evs[event, 2] = 18
            # if event is a B-cue
            elif new_evs[event, 2] in {6, 7, 8, 9, 10}:
                # recode as too soon B-cue
                new_evs[event, 2] = 19

            # look for next probe
            i = 2
            while new_evs[event + i, 2] not in {11, 12, 13, 14, 15, 16}:
                # if no probe is found before next cue is presented
                if new_evs[event + i, 2] in {5, 6, 7, 8, 9, 10}:
                    broken.append(trial)
                    break
                i += 1

            # if probe is an X
            if new_evs[event + i, 2] == 11:
                # recode as too soon X-probe
                new_evs[event + i, 2] = 20
            # if probe is an Y
            elif new_evs[event + i, 2] in {12, 13, 14, 15, 16}:
                # recode as too soon Y-probe
                new_evs[event + i, 2] = 21

            # save trial information as NaN
            trial += 1
            rt.append(np.nan)
            reaction.append(np.nan)
            # go on to next trial
            continue

        # --- 2nd check: if next event is a probe stimulus ---
        elif new_evs[event + 1, 2] in {11, 12, 13, 14, 15, 16}:

            # if event after probe is a reaction
            if new_evs[event + 2, 2] in {1, 2, 3, 4}:

                # save reaction time
                rt.append(
                    (new_evs[event + 2, 0] - new_evs[event + 1, 0]) / sfreq)

                # if reaction is correct
                if new_evs[event + 2, 2] in {3, 4}:

                    # save response
                    reaction.append('Correct')

                    # if cue was an A
                    if new_evs[event, 2] == 5:
                        # recode as correct A-cue
                        new_evs[event, 2] = 22

                        # if probe was an X
                        if new_evs[event + 1, 2] == 11:
                            # recode as correct AX probe combination
                            new_evs[event + 1, 2] = 23

                        # if probe was a Y
                        else:
                            # recode as correct AY probe combination
                            new_evs[event + 1, 2] = 24

                        # go on to next trial
                        trial += 1
                        continue

                    # if cue was a B
                    else:
                        # recode as correct B-cue
                        new_evs[event, 2] = 25

                        # if probe was an X
                        if new_evs[event + 1, 2] == 11:
                            # recode as correct BX probe combination
                            new_evs[event + 1, 2] = 26
                        # if probe was a Y
                        else:
                            # recode as correct BY probe combination
                            new_evs[event + 1, 2] = 27

                        # go on to next trial
                        trial += 1
                        continue

                # if reaction was incorrect
                else:

                    # save response
                    reaction.append('Incorrect')

                    # if cue was an A
                    if new_evs[event, 2] == 5:
                        # recode as incorrect A-cue
                        new_evs[event, 2] = 28

                        # if probe was an X
                        if new_evs[event + 1, 2] == 11:
                            # recode as incorrect AX probe combination
                            new_evs[event + 1, 2] = 29

                        # if probe was a Y
                        else:
                            # recode as incorrect AY probe combination
                            new_evs[event + 1, 2] = 30

                        # go on to next trial
                        trial += 1
                        continue

                    # if cue was a B
                    else:
                        # recode as incorrect B-cue
                        new_evs[event, 2] = 31

                        # if probe was an X
                        if new_evs[event + 1, 2] == 11:
                            # recode as incorrect BX probe combination
                            new_evs[event + 1, 2] = 32

                        # if probe was a Y
                        else:
                            # recode as incorrect BY probe combination
                            new_evs[event + 1, 2] = 33

                        # go on to next trial
                        trial += 1
                        continue

            # if no reaction followed cue-probe combination
            elif new_evs[event + 2, 2] not in {1, 2, 3, 4}:

                # save reaction time as NaN
                rt.append(np.nan)
                reaction.append(np.nan)

                # if cue was an A
                if new_evs[event, 2] == 5:
                    # recode as missed A-cue
                    new_evs[event, 2] = 34

                    # if probe was an X
                    if new_evs[event + 1, 2] == 11:
                        # recode as missed AX probe combination
                        new_evs[event + 1, 2] = 35

                    # if probe was a Y
                    else:
                        # recode as missed AY probe combination
                        new_evs[event + 1, 2] = 36

                    # go on to next trial
                    trial += 1
                    continue

                # if cue was a B
                else:
                    # recode as missed B-cue
                    new_evs[event, 2] = 37

                    # if probe was an X
                    if new_evs[event + 1, 2] == 11:
                        # recode as missed BX probe combination
                        new_evs[event + 1, 2] = 38

                    # if probe was a Y
                    else:
                        # recode as missed BY probe combination
                        new_evs[event + 1, 2] = 39

                    # go on to next trial
                    trial += 1
                    continue

    # skip other events
    else:
        continue

###############################################################################
# 4) Set descriptive event names for extraction of epochs

# cue events
cue_event_id = {'Too_soon A': 18,
                'Too_soon B': 19,

                'Correct A': 22,
                'Correct B': 25,

                'Incorrect A': 28,
                'Incorrect B': 31,

                'Missed A': 34,
                'Missed B': 37}

# probe events
probe_event_id = {'Too_soon X': 20,
                  'Too_soon Y': 21,

                  'Correct AX': 23,
                  'Correct AY': 24,

                  'Correct BX': 26,
                  'Correct BY': 27,

                  'Incorrect AX': 29,
                  'Incorrect AY': 30,

                  'Incorrect BX': 32,
                  'Incorrect BY': 33,

                  'Missed AX': 35,
                  'Missed AY': 36,

                  'Missed BX': 38,
                  'Missed BY': 39}

###############################################################################
# 5) Create metadata structure to be added to the epochs

# only keep cue events
cue_events = new_evs[np.where((new_evs[:, 2] == 18) |
                              (new_evs[:, 2] == 19) |
                              (new_evs[:, 2] == 22) |
                              (new_evs[:, 2] == 25) |
                              (new_evs[:, 2] == 28) |
                              (new_evs[:, 2] == 31) |
                              (new_evs[:, 2] == 34) |
                              (new_evs[:, 2] == 37))]

# only keep probe events
probe_events = new_evs[np.where((new_evs[:, 2] == 20) |
                                (new_evs[:, 2] == 21) |
                                (new_evs[:, 2] == 23) |
                                (new_evs[:, 2] == 24) |
                                (new_evs[:, 2] == 26) |
                                (new_evs[:, 2] == 27) |
                                (new_evs[:, 2] == 29) |
                                (new_evs[:, 2] == 30) |
                                (new_evs[:, 2] == 32) |
                                (new_evs[:, 2] == 33) |
                                (new_evs[:, 2] == 35) |
                                (new_evs[:, 2] == 36) |
                                (new_evs[:, 2] == 38) |
                                (new_evs[:, 2] == 39))]

# reversed event_id dict
cue_event_id_rev = {val: key for key, val in cue_event_id.items()}
probe_event_id_rev = {val: key for key, val in probe_event_id.items()}

# check if events shape match
if cue_events.shape[0] != probe_events.shape[0]:
    cue_events = np.delete(cue_events, broken, 0)

# create list with reactions based on cue and probe event ids
same_stim, reaction_cues, reaction_probes, cues, probes = [], [], [], [], []
for cue, probe in zip(cue_events[:, 2], probe_events[:, 2]):
    response, cue = cue_event_id_rev[cue].split(' ')
    reaction_cues.append(response)
    # save cue
    cues.append(cue)

    # save response
    response, probe = probe_event_id_rev[probe].split(' ')
    reaction_probes.append(response)

    # check if same type of combination was shown in the previous trail
    if len(probes):
        stim = same_stim[-1]
        if probe == probes[-1] \
                and response == 'Correct' \
                and reaction_probes[-2] == 'Correct':
            stim += 1
            same_stim.append(stim)
        else:
            same_stim.append(0)
    else:
        stim = 0
        same_stim.append(0)

    # save probe
    probes.append(probe)

# subject's experimental group
if subject in [1, 4, 7, 8, 10, 13, 18, 19, 20, 21, 22, 23, 25]:
    group = 1
elif subject in [2, 3, 6, 9, 11, 12, 14, 15, 16, 17, 24, 26, 27]:
    group = 2
elif subject in [5]:
    group = np.nan

# create data frame with epochs metadata
metadata = {'group': np.delete(np.repeat(group, trial),  broken, 0),
            'block': np.delete(block, broken, 0),
            'trial': np.delete(np.arange(0, trial), broken, 0),
            'cue': cues,
            'probe': probes,
            'run': same_stim,
            'reaction_cues': reaction_cues,
            'reaction_probes': reaction_probes,
            'rt': np.delete(rt, broken, 0)}
metadata = pd.DataFrame(metadata)

# add subject id ands save RT measures for later analyses
rt_data = metadata.copy()
rt_data = rt_data.assign(subject=subject)

# export to .tsv
rt_data.to_csv(op.join(fname.rt, 'sub-%s-rt.tsv' % subject),
               sep='\t',
               index=False)

###############################################################################
# 6) Extract the epochs

# rejection threshold
reject = dict(eeg=250e-6)
decim = 1

if raw.info['sfreq'] == 256.0:
    decim = 2
elif raw.info['sfreq'] == 512.0:
    decim = 4
elif raw.info['sfreq'] == 1024.0:
    decim = 8

# extract cue epochs
cue_epochs = Epochs(raw, cue_events, cue_event_id,
                    metadata=metadata,
                    on_missing='ignore',
                    tmin=-2.0,
                    tmax=2.5,
                    baseline=None,
                    preload=True,
                    reject_by_annotation=True,
                    reject=reject,
                    decim=decim)

# extract probe epochs
probe_epochs = Epochs(raw, probe_events, probe_event_id,
                      metadata=metadata,
                      on_missing='ignore',
                      tmin=-3.0,
                      tmax=1.5,
                      baseline=None,
                      preload=True,
                      reject_by_annotation=True,
                      reject=reject,
                      decim=decim
                      )

###############################################################################
# 7) Save info about extracted and rejected epochs

# clean cue epochs
clean_cues = cue_epochs.selection
bad_cues = [x for x in set(list(range(0, trial)))
            if x not in set(cue_epochs.selection)]
# clean probe epochs
clean_probes = probe_epochs.selection
bad_probes = [x for x in set(list(range(0, trial)))
              if x not in set(probe_epochs.selection)]

###############################################################################
# 8) Save epochs

# output path for cues
cue_output_path = fname.output(processing_step='cue_epochs',
                               subject=subject,
                               file_type='epo.fif')
# output path for epochs
probe_output_path = fname.output(processing_step='probe_epochs',
                                 subject=subject,
                                 file_type='epo.fif')

# resample and save cue epochs to disk
cue_epochs.save(cue_output_path, overwrite=True)

# also save probe epochs to disk
probe_epochs.save(probe_output_path, overwrite=True)

###############################################################################
# 9) Create HTML report
epochs_summary = '<p>Cue epochs extracted: <br>' \
                 'A: %s <br>' \
                 'B: %s <br>' \
                 '<p>Probe epochs extracted: <br>' \
                 'AX: %s <br>' \
                 'AY: %s <br>' \
                 'BX: %s <br>' \
                 'BY: %s <br>' \
                 % (
                     len(cue_epochs['Correct A']),
                     len(cue_epochs['Correct B']),
                     len(probe_epochs['Correct AX']),
                     len(probe_epochs['Correct AY']),
                     len(probe_epochs['Correct BX']),
                     len(probe_epochs['Correct BY'])
                    )

with open_report(fname.report(subject=subject)[0]) as report:
    report.add_htmls_to_section(htmls=epochs_summary,
                                captions='Epochs summary',
                                section='Epochs',
                                replace=True)
    report.save(fname.report(subject=subject)[1], overwrite=True,
                open_browser=False)

###############################################################################
# 10) Save time bins for further analysis
#
# # extract corrects
# correct_epochs = cue_epochs['Correct A', 'Correct B']
# # apply baseline correction
# correct_epochs.apply_baseline((-0.300, -0.050))
#
# # epochs to pandas dataframe
# df = correct_epochs.to_data_frame(long_format=True)
#
# # get time roi N170
# n170 = df[((df["time"] >= 150) & (df["time"] <= 250))
#           & ((df["channel"] == 'PO8') |
#              (df["channel"] == 'PO7') |
#              (df["channel"] == 'FCz'))]
# n170 = n170.assign(subject=subject)
#
# # get time roi LPC
# LPC = df[((df["time"] >= 400) & (df["time"] <= 750))
#          & (df["channel"] == 'Pz')]
# LPC = LPC.assign(subject=subject)
#
# # get time roi CNV
# CNV = df[((df["time"] >= 1000) & (df["time"] <= 1750))
#          & (df["channel"] == 'C1')]
# CNV = CNV.assign(subject=subject)
#
# # export to .tsv files
# n170.to_csv(op.join(fname.rois,
#                     'sub-%s-n170.tsv' % subject),
#             sep='\t',
#             index=False)
# LPC.to_csv(op.join(fname.rois,
#                    'sub-%s-LPC.tsv' % subject),
#            sep='\t',
#            index=False)
#
# CNV.to_csv(op.join(fname.rois,
#                    'sub-%s-cnv.tsv' % subject),
#            sep='\t',
#            index=False)
