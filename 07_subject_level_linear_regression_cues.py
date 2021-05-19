"""
=============================================
Fit single subject linear model to cue epochs
=============================================

Mass-univariate analysis of cue evoked activity.

Authors: José C. García Alanis <alanis.jcg@gmail.com>

License: BSD (3-clause)
"""
import numpy as np
import patsy

from sklearn.utils.class_weight import compute_sample_weight
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from mne import read_epochs
from mne.decoding import Vectorizer, get_coef

# All parameters are defined in config.py
from config import subjects, fname, LoggingFormat, n_jobs

# dicts for storing individual sets of epochs/ERPs
cues = dict()

# baseline to be applied
baseline = (-0.300, -0.050)

###############################################################################
# 1) loop through subjects and compute ERPs for A and B cues
for subj in subjects:

    # remove subject 5 from analysis
    if subj == 5:
        continue

    # log progress
    print(LoggingFormat.PURPLE +
          LoggingFormat.BOLD +
          'Loading epochs for subject %s' % subj +
          LoggingFormat.END)

    # import the output from previous processing step
    input_file = fname.output(subject=subj,
                              processing_step='cue_epochs',
                              file_type='epo.fif')
    cue_epo = read_epochs(input_file, preload=True)

    # extract a and b epochs (only those with correct responses)
    # and apply baseline
    cues['subj_%s' % subj] = cue_epo['Correct A', 'Correct B']
    cues['subj_%s' % subj].apply_baseline(baseline).crop(tmin=-0.500)

###############################################################################
# 2) linear model parameters
# use first subject as generic information template for results
generic = cues['subj_%s' % subjects[0]].copy()

# save the generic info structure of cue epochs (i.e., channel names, number of
# channels, etc.).
epochs_info = generic.info

# only use times > -0.25
times_to_use = (generic.times >= -0.25) & (generic.times <= 2.45)
times = generic.times[times_to_use]
n_times = len(times)
n_channels = len(epochs_info['ch_names'])

# subjects
subjects = list(cues.keys())

# independent variables to be used in the analysis (i.e., predictors)
predictors = ['cue', 'block']

# compute number of predictors minus the intercept
# here: cue A before manipulation is set as reference level
# therefore design matrixs shows:
# 1) effect of cue B against cue A in baseline
# 2) effect of manipulation on cue A
# 3) effect of manipulation on cue B
generic_design = generic.metadata[predictors]
generic_design = patsy.dmatrix("cue + C(block, Treatment('Single')) + "
                               "cue:C(block, Treatment('Single'))",
                               generic_design,
                               return_type='dataframe')
n_predictors = len(generic_design.columns) - 1

###############################################################################
# 3) initialise place holders for the storage of results
betas = np.zeros((n_predictors,
                  len(cues.values()),
                  n_channels * n_times))

r_squared = np.zeros((len(cues.values()),
                      n_channels * n_times))

###############################################################################
# 4) Fit linear model for each subject
for subj_ind, subj in enumerate(cues):
    print(subj_ind, subj)

    # 4.1) create subject design matrix using epochs metadata
    metadata = cues[subj].metadata.copy()

    # only keep predictor columns
    design = metadata[predictors]
    # design['block'].cat.reorder_categories(['Single', 'Dual'], inplace=True)

    # create design matrix
    design = patsy.dmatrix("cue + C(block, Treatment('Single')) + "
                           "cue:C(block, Treatment('Single'))", design,
                           return_type='dataframe')
    design.columns = ['Intercept', 'cue[T.B]', 'block[T.Dual]',
                      'cue[T.B]:block[T.Dual]']
    design = design[['cue[T.B]', 'block[T.Dual]', 'cue[T.B]:block[T.Dual]']]

    # 4.2) vectorise channel data for linear regression
    # data to be analysed
    dat = cues[subj].get_data()
    dat = dat[:, :, times_to_use]
    Y = Vectorizer().fit_transform(dat)

    # 4.3) fit linear model with sklearn's LinearRegression
    weights = compute_sample_weight(class_weight='balanced',
                                    y=metadata.cue.to_numpy())
    linear_model = LinearRegression(n_jobs=n_jobs, fit_intercept=True)
    linear_model.fit(design, Y, sample_weight=weights)

    # 4.4) extract the resulting coefficients (i.e., betas)
    # extract betas
    coefs = get_coef(linear_model, 'coef_')

    # 4.5) extract model r_squared
    r2 = r2_score(Y, linear_model.predict(design),
                  multioutput='raw_values')
    # save model R-squared
    r_squared[subj_ind, :] = r2

    # save results
    for pred_i, predictor in enumerate(design.columns):
        print(pred_i, predictor)
        if 'cue' in predictor and 'block' not in predictor:
            # extract cue beats
            betas[pred_i, subj_ind, :] = coefs[:, pred_i]
        elif 'cue' not in predictor and 'block' in predictor:
            # extract cue beats
            betas[pred_i, subj_ind, :] = coefs[:, pred_i]
        elif 'cue' in predictor and 'block' in predictor:
            # extract cue beats
            betas[pred_i, subj_ind, :] = coefs[:, pred_i]

###############################################################################
# 5) Save subject-level results to disk
np.save(fname.results + '/subj_betas_cue_by_block_m250_robust.npy', betas)
np.save(fname.results + '/subj_r2_cue_m250_robust.npy', r_squared)
