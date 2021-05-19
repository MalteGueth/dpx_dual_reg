"""
===============================================
Repair EEG artefacts caused by ocular movements
===============================================

Identify "bad" components in ICA solution (e.g., components which are highly
correlated the time course of the electrooculogram).

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import numpy as np
import matplotlib.pyplot as plt

from mne import open_report, events_from_annotations, Epochs
from mne.io import read_raw_fif
from mne.preprocessing import read_ica, corrmap

# All parameters are defined in config.py
from config import fname, parser, LoggingFormat

# Handle command line arguments
args = parser.parse_args()
subject = args.subject

print(LoggingFormat.PURPLE +
      LoggingFormat.BOLD +
      'Finding and removing bad components for subject %s' % subject +
      LoggingFormat.END)

###############################################################################
# 1) Import the output from previous processing step
input_file = fname.output(subject=subject,
                          processing_step='repair_bads',
                          file_type='raw.fif')
raw = read_raw_fif(input_file, preload=True)
raw.apply_proj()

###############################################################################
# 2) Import ICA weights from precious processing step
ica_file = fname.output(subject=subject,
                        processing_step='fit_ica',
                        file_type='ica.fif')
ica = read_ica(ica_file)

###############################################################################
# 3) Find bad components via correlation with template ICA
temp_subjs = [5]
temp_raws = []
temp_icas = []

# import template subjects
for subj in temp_subjs:
    temp_raws.append(read_raw_fif(fname.output(subject=subj,
                                               processing_step='repair_bads',
                                               file_type='raw.fif')))
    temp_icas.append(read_ica(fname.output(subject=subj,
                                           processing_step='fit_ica',
                                           file_type='ica.fif')))

# compute correlations with template ocular movements up/down and left/right
corrmap(icas=[temp_icas[0], ica],
        template=(0, 0), threshold=0.90, label='blink_up', plot=False)
corrmap(icas=[temp_icas[0], ica],
        template=(0, 1), threshold=0.85, label='blink_side', plot=False)

# # compute correlations with template ocular movements that look slightly
# # different
# corrmap(icas=[temp_icas[1], ica],
#         template=(0, 0), threshold=0.90, label='blink_misc', plot=False)
# corrmap(icas=[temp_icas[0], ica],
#         template=(0, 1), threshold=0.90, label='blink_misc', plot=False)

###############################################################################
# 4) Create summary plots to show signal correction on main experimental
# condition

# create a-cue epochs
a_evs = events_from_annotations(raw, regexp='^(70)')[0]
a_epo = Epochs(raw, a_evs,
               tmin=-1.0,
               tmax=2.0,
               reject_by_annotation=True,
               proj=False,
               preload=True)
a_epo.apply_baseline(baseline=(-0.3, -0.05))
a_evo = a_epo.average()

# loop over identified "bad" components
bad_components = []
for label in ica.labels_:
    bad_components.extend(ica.labels_[label])

for bad_comp in np.unique(bad_components):
    # show component frequency spectrum
    fig_comp = ica.plot_properties(a_epo,
                                   picks=bad_comp,
                                   psd_args={'fmax': 35.},
                                   show=False)[0]

    # show how the signal is affected by component rejection
    fig_evoked = ica.plot_overlay(a_evo, exclude=[bad_comp], show=False)
    plt.close(fig_evoked)

    # create HTML report
    with open_report(fname.report(subject=subject)[0]) as report:
        report.add_figs_to_section(fig_comp, 'Component %s identified ' 
                                             'by correlation with template'
                                   % bad_comp,
                                   section='ICA',
                                   replace=True)
        report.add_figs_to_section(fig_evoked, 'Component %s rejected'
                                   % bad_comp,
                                   section='ICA',
                                   replace=True)
        report.save(fname.report(subject=subject)[1], overwrite=True,
                    open_browser=False)

# add bad components  to exclusion list
ica.exclude = np.unique(bad_components)

# apply ica weights to data
ica.apply(raw)

###############################################################################
# 5) Save repaired data set
# output path
output_path = fname.output(processing_step='repaired_with_ica',
                           subject=subject,
                           file_type='raw.fif')

# save file
raw.save(output_path, overwrite=True)
