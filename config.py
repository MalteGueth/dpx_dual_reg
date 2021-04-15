"""
========================
Study configuration file
========================

Configuration parameters and global variable values for the study.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import os
from os import path as op
import platform

import argparse
import numpy as np

from utils import FileNames

from mne.channels import make_standard_montage


###############################################################################
class LoggingFormat:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


###############################################################################
# User parser to handle command line arguments
parser = argparse.ArgumentParser(description='Parse command line argument for '
                                             'pre-processing of EEG data.')
parser.add_argument('subject',
                    metavar='sub###',
                    help='The subject to process',
                    type=int)

# Determine which user is running the scripts on which machine. Set the path to
# where the data is stored and determine how many CPUs to use for analysis.
node = platform.node()  # Maschine
system = platform.system()  # Os

# You want to add your machine to this list
if 'Jose' in node and 'n' in system:
    # iMac at work
    data_dir = '../data'
    n_jobs = 2  # This station has 4 cores (we'll use 2).
elif 'jose' in node and 'x' in system:
    # pc at home
    data_dir = '../data'
    n_jobs = 8  # This station has 16 cores (we'll use 8).
elif 'ma04' in node:
    data_dir = '../data'
    n_jobs = 2
else:
    # Defaults
    data_dir = '../data'
    n_jobs = 1

# For BLAS to use the right amount of cores
os.environ['OMP_NUM_THREADS'] = str(n_jobs)

###############################################################################
# Relevant parameters for the analysis.
sampling_rate = 256.0  # Hz
task_name = 'dpxtt'
task_description = 'DPX, effects of time on task'
# eeg channel names and locations
montage = make_standard_montage('standard_1020')
# channels to be exclude from import
exclude = ['EXG5', 'EXG6', 'EXG7', 'EXG8']

# subjects to use for analysis
subjects = np.arange(1, 28)

# relevant events in the paradigm
event_ids = {'correct_target_button': 13,
             'correct_non_target_button': 12,
             'incorrect_target_button': 113,
             'incorrect_non_target_button': 112,
             'cue_0': 70,
             'cue_1': 71,
             'cue_2': 72,
             'cue_3': 73,
             'cue_4': 74,
             'cue_5': 75,
             'probe_0': 76,
             'probe_1': 77,
             'probe_2': 78,
             'probe_3': 79,
             'probe_4': 80,
             'probe_5': 81,
             'start_record': 127,
             'pause_record': 245}

###############################################################################
# Templates for filenames
#
# This part of the config file uses the FileNames class. It provides a small
# wrapper around string.format() to keep track of a list of filenames.
# See fnames.py for details on how this class works.
fname = FileNames()

# directories to use for input and output
fname.add('data_dir', data_dir)
fname.add('bids_data', '{data_dir}/sub-{subject:03d}')
fname.add('subject_demographics',
          '{data_dir}/subject_data/subject_demographics.tsv')
fname.add('sourcedata_dir', '{data_dir}/sourcedata')
fname.add('derivatives_dir', '{data_dir}/derivatives')
fname.add('reports_dir', '{derivatives_dir}/reports')
fname.add('results', '{derivatives_dir}/results')
fname.add('rt', '{results}/rt')
fname.add('figures', '{results}/figures')
fname.add('tables', '{results}/tables')
fname.add('rois', '{results}/rois')


# The paths for data file input
# fname.add('source',
#           '{sourcedata_dir}/sub-{subject:02d}/eeg/sub-{subject:02d}_dpx_eeg.bdf')  # noqa
# alternative:
def source_file(files, source_type, subject):
    if source_type == 'eeg':
        return \
            op.join(files.sourcedata_dir,
                    'sub-%02d/%s/sub-%02d_dualreg_eeg.bdf' % (subject,
                                                              source_type,
                                                              subject))
    elif source_type == 'demo':
        return \
            op.join(files.sourcedata_dir,
                    'sub-%02d/%s/sub-%02d_dualreg_demographics.tsv' % (subject,
                                                                       source_type,
                                                                       subject))


# create full path for data file input
fname.add('source', source_file)


# The paths that are produced by the analysis steps
def output_path(path, processing_step, subject, file_type):
    path = op.join(path.derivatives_dir, processing_step, 'sub-%03d' % subject)
    os.makedirs(path, exist_ok=True)
    return op.join(path,
                   'sub-%03d-%s-%s' % (subject, processing_step, file_type))


# The full path for data file output
fname.add('output', output_path)


# The paths that are produced by the report step
def report_path(path, subject):
    h5_path = op.join(path.reports_dir, 'sub-%03d.h5' % subject)
    html_path = op.join(path.reports_dir, 'sub-%03d-report.html' % subject)
    return h5_path, html_path


# The full path for report file output
fname.add('report', report_path)

# Files produced by system check and validator
fname.add('system_check', './system_check.txt')
fname.add('validator', './validator.txt')
