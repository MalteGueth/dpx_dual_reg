##### generalized decoding across time for cue data

import glob
import os

import numpy as np

import mne

from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from mne.decoding import Vectorizer, SlidingEstimator, cross_val_multiscore, GeneralizingEstimator

path = './epochs/'

### start the decoding for baseline

for file in glob.glob(os.path.join(path, '*base-epo.fif')):
    
    # Read the epoched eeg files
    epochs_base = mne.read_epochs(file, preload=True)
    
    # To make the exmaple run faster, downsample the data
    epochs_base.resample(256)

    # Cut the epochs to the required length
    epochs_base.crop(tmin=-0.25, tmax=epochs_base.tmax)
    # Equalize the trial counts to balance conditions
    epochs_base_eq = epochs_base.copy().equalize_event_counts(['A', 'B'])[0]
    
    # Pick the EEG data for decoding (X) ...
    X=epochs_base_eq['A','B'].get_data()
    # ... and the categories / experimental conditions (y)
    y=epochs_base_eq['A','B'].events[:,2]
    
    # Create the pipeline with linear support vector classification
    # and balanced classes. Alternatively, instead of equlizing trial
    # counts, this argument can be set to unbalanced.
    clf = make_pipeline(Vectorizer(), StandardScaler(),
                        LinearSVC(class_weight='balanced')
                        )
    clf.fit(X, y)
    
    # Calculate scores for classification
    sl = SlidingEstimator(clf)
    scores_time_decoding = cross_val_multiscore(sl, X, y)
    
    # Append the results for each subject
    if file == './epochs/101_base-epo.fif':
        scores_td_base = scores_time_decoding
    else:
        scores_td_base = np.append(scores_td_base, scores_time_decoding, axis=0)
    
    # Again, calculate scores with a receiver operating curve
    gen = GeneralizingEstimator(clf, scoring='roc_auc')
    scores_gat = cross_val_multiscore(gen, X, y)
    
    if file == './epochs/101_base-epo.fif':
        scores_gat_base = scores_gat
    else:
        scores_gat_base = np.append(scores_gat_base, scores_gat, axis=0)
        
for file in glob.glob(os.path.join(path, '*reg-epo.fif')):
    
    epochs_reg = mne.read_epochs(file, preload=True)
    
    epochs_base.resample(256)
    
    epochs_base.crop(tmin=-0.25, tmax=epochs_base.tmax)
    epochs_reg_eq = epochs_reg.copy().equalize_event_counts(['A', 'B'])[0]

    X=epochs_reg_eq['A','B'].get_data()
    y=epochs_reg_eq['A','B'].events[:,2]
    
    clf = make_pipeline(Vectorizer(), StandardScaler(),
                        LinearSVC(class_weight='balanced')
                        )
    clf.fit(X, y)
    
    sl = SlidingEstimator(clf, scoring='roc_auc')
    scores_time_decoding = cross_val_multiscore(sl, X, y)
    
    if file == './epochs/101_reg-epo.fif':
        scores_td_reg = scores_time_decoding
    else:
        scores_td_reg = np.append(scores_td_reg, scores_time_decoding, axis=0)  
        
    gen = GeneralizingEstimator(clf)
    scores_gat = cross_val_multiscore(gen, X, y)
    
    if file == './epochs/101_reg-epo.fif':
        scores_gat_reg = scores_gat
    else:
        scores_gat_reg = np.append(scores_gat_reg, scores_gat, axis=0)  
        
###### Plot decoding results

import matplotlib.pyplot as plt
      
fig, ax = plt.subplots()
ax.plot(epochs_base.times, scores_td_base.T)
plt.show()

fig, ax = plt.subplots()
ax.plot(epochs_base.times, scores_td_base.mean(0))
plt.show()

tmin, tmax = epochs_base.times[[0, -1]]

fig, ax = plt.subplots()
im = ax.imshow(
    scores_gat_base.mean(0),
    origin="lower", cmap="RdBu_r",
    extent=(tmin, tmax, tmin, tmax),
    vmax=0.7, vmin=0.3);
ax.axhline(0., color='k')
ax.axvline(0., color='k')
plt.colorbar(im)



fig, ax = plt.subplots()
ax.plot(epochs_base.times, scores_td_reg.T)
plt.show()

fig, ax = plt.subplots()
ax.plot(epochs_base.times, scores_td_reg.mean(0))
plt.show()

tmin, tmax = epochs_base.times[[0, -1]]

fig, ax = plt.subplots()
im = ax.imshow(
    scores_gat_reg.mean(0),
    origin="lower", cmap="RdBu_r",
    extent=(tmin, tmax, tmin, tmax),
    vmax=0.7, vmin=0.3);
ax.axhline(0., color='k')
ax.axvline(0., color='k')
plt.colorbar(im)
