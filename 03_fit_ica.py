"""
==============================
Independent Component Analysis
==============================

Decompose EEG signal into independent components

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
from mne import open_report
from mne.io import read_raw_fif
from mne.preprocessing import ICA

# All parameters are defined in config.py
from config import fname, parser, LoggingFormat, n_jobs

# Handle command line arguments
args = parser.parse_args()
subject = args.subject

print(LoggingFormat.PURPLE +
      LoggingFormat.BOLD +
      'Fit ICA for subject %s' % subject +
      LoggingFormat.END)

###############################################################################
# 1) Import the output from previous processing step
input_file = fname.output(subject=subject,
                          processing_step='repair_bads',
                          file_type='raw.fif')
raw = read_raw_fif(input_file, preload=True)
# activate average reference
raw.apply_proj()

# filter data to remove drifts
raw_filt = raw.copy().filter(l_freq=1.0, h_freq=None, n_jobs=n_jobs)

###############################################################################
#  2) Set ICA parameters
n_components = 15
method = 'infomax'
reject = dict(eeg=250e-6)

###############################################################################
#  2) Fit ICA
ica = ICA(n_components=n_components,
          method=method,
          fit_params=dict(extended=True))

ica.fit(raw_filt,
        reject=reject,
        reject_by_annotation=True)

###############################################################################
# 3) Plot ICA components
ica_fig = ica.plot_components(picks=range(0, 20), show=False)

###############################################################################
# 4) Save ICA solution
# output path
output_path = fname.output(processing_step='fit_ica',
                           subject=subject,
                           file_type='ica.fif')
# save file
ica.save(output_path)

###############################################################################
# 5) Create HTML report
with open_report(fname.report(subject=subject)[0]) as report:
    report.add_figs_to_section(ica_fig, 'ICA solution',
                               section='ICA',
                               replace=True)
    report.save(fname.report(subject=subject)[1], overwrite=True,
                open_browser=False)
