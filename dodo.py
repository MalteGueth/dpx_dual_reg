"""
Do-it script to execute the entire pipeline using the doit tool:

All the filenames are defined in config.py

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)

Notes
-----
- for more on doit: http://pydoit.org
"""
from config import fname, subjects

# Configuration for the "doit" tool.
DOIT_CONFIG = dict(
    # While running scripts, output everything the script is printing to the
    # screen.
    verbosity=2,

    # When the user executes "doit list", list the tasks in the order they are
    # defined in this file, instead of alphabetically.
    sort='definition'
)


def task_check():
    """Check the system dependencies."""
    return dict(
        file_dep=['check_system.py'],
        targets=[fname.system_check],
        actions=['python check_system.py %s' % fname.data_dir]
    )


# This task executes a single analysis script for each subject, giving
# the subject as a command line parameter to the script.
def task_eeg_to_bids():
    """Step 00: Put EEG data into a BIDS-compliant directory structure."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `task_check`
            task_dep=['check'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=[fname.source(source_type='eeg', subject=subject),
                      '00_eeg_to_bids.py'],

            # The files produced by the script
            targets=[fname.bids_data(subject=subject)],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 00_eeg_to_bids.py %s' % subject]
        )


# This task executes a single analysis script for each subject, giving
# the subject as a command line parameter to the script.
def task_task_blocks():
    """Step 01: Extracts task segments (drop pauses in between)."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `eeg_to_bids`
            task_dep=['eeg_to_bids'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=['00_eeg_to_bids.py'],

            # The files produced by the script
            targets=[fname.output(processing_step='task_blocks',
                                  subject=subject,
                                  file_type='raw.fif')],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 01_task_blocks.py %s' % subject]
        )


# This task executes a single analysis script for each subject, giving
# the subject as a command line parameter to the script.
def task_repair_bad_channels():
    """Step 02: Identify and repair bad (i.e., noisy) EEG channels."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `task_blocks`
            task_dep=['task_blocks'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=['00_eeg_to_bids.py',
                      '01_task_blocks.py'],

            # The files produced by the script
            targets=[fname.output(processing_step='repair_bads',
                                  subject=subject,
                                  file_type='raw.fif')],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 02_repair_bad_eeg_channels.py %s' % subject]
        )


# This task executes a single analysis script for each subject, giving
# the subject as a command line parameter to the script.
def task_fit_ica():
    """Step 03: Decompose EEG signal into independent components."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `repair_bad_channels`
            task_dep=['repair_bad_channels'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=['00_eeg_to_bids.py',
                      '01_task_blocks.py',
                      '02_repair_bad_eeg_channels.py'],

            # The files produced by the script
            targets=[fname.output(processing_step='fit_ica',
                                  subject=subject,
                                  file_type='ica.fif')],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 03_fit_ica.py %s' % subject]
        )


def task_repair_artefacts():
    """Step 04: Repair EEG artefacts caused by ocular movements."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `fit_ica`
            task_dep=['fit_ica'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=['00_eeg_to_bids.py',
                      '01_task_blocks.py',
                      '02_repair_bad_eeg_channels.py',
                      '03_fit_ica.py'],

            # The files produced by the script
            targets=[fname.output(processing_step='repaired_with_ica',
                                  subject=subject,
                                  file_type='raw.fif')],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 04_repair_eeg_artefacts.py %s' % subject]
        )


def task_extract_epochs():
    """Step 05: Extract epochs from continuous EEG."""
    # Run the script for each subject in a sub-task.
    for subject in subjects:
        yield dict(
            # This task should come after `fit_ica`
            task_dep=['repair_artefacts'],

            # A name for the sub-task: set to the name of the subject
            name=subject,

            # If any of these files change, the script needs to be re-run. Make
            # sure that the script itself is part of this list!
            file_dep=['00_eeg_to_bids.py',
                      '01_task_blocks.py',
                      '02_repair_bad_eeg_channels.py',
                      '03_fit_ica.py',
                      '04_repair_eeg_artefacts.py'],

            # The files produced by the script
            targets=[fname.output(processing_step='cue_epochs',
                                  subject=subject,
                                  file_type='epo.fif'),
                     fname.output(processing_step='probe_epochs',
                                  subject=subject,
                                  file_type='epo.fif')],

            # How the script needs to be called. Here we indicate it should
            # have one command line parameter: the name of the subject.
            actions=['python 05_extract_epochs.py %s' % subject]
        )

# # # Here is another example task that averages across subjects.
# # def task_example_summary():
# #     """Step 01: Average across subjects."""
# #     return dict(
# #         task_dep=['example_step'],  # This task should come after
# #         `task_example_step`
# #         file_dep=[fname.output(subject=s) for s in subjects] + [
# #         '01_grand_average.py'],
# #         targets=[fname.grand_average],
# #         actions=['python 01_grand_average.py'],
# #     )
# #
# #
# # def task_figures():
# #     """Make all figures. Each figure is a sub-task."""
# #     # Make figure 1
# #     yield dict(
# #         name='figure_example1',
# #         file_dep=['figure_example1.py'],
# #         targets=[fname.figure1],
# #         actions=['python figure_example1.py'],
# #     )
# #
# #     # Make figure 2
# #     yield dict(
# #         name='figure_grand_average',
# #         file_dep=[fname.grand_average, 'figure_grand_average.py'],
# #         targets=[fname.figure_grand_average],
# #         actions=['python figure_grand_average.py'],
# #     )
