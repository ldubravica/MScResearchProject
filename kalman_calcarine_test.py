"""
Train and test minimal predictive processing model to reproduce LZ and TE changes in
the psychedelic state and schizophrenia patients. The model consists of a Kalman filter
trained on time series from the placebo condition of the LSD_MEG dataset, and LSD/SCH
are modelled as increased prior/likelihood variance, respectively.

Pedro Mediano and Hardik Rajpal, Apr 2021
"""

#from glob import glob
from itertools import product
# from lz76 import LZ76
from lz76.lz76 import LZ76
from pykalman import KalmanFilter
from scipy.io import loadmat
from scipy.signal import butter, filtfilt, decimate
from tqdm import tqdm

import h5py
import jpype
import matplotlib.pyplot as pl
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
import seaborn as sn
import sys

project_root = os.path.dirname(os.path.abspath(__file__))


def startJVM():
    """
    Start JVM using jpype, if not already running. Assumes the JIDT jar is in
    the current folder.
    """
    jarLocation = os.path.join(project_root, "infodynamics-dist-1.6.1", "infodynamics.jar")
    if not jpype.isJVMStarted():
        jpype.startJVM(jpype.getDefaultJVMPath(), "-ea","-Xmx1024m", "-Djava.class.path=" + jarLocation)


def TE(src, tgt, k=1, tau=1):
    """
    Computes transfer entropy using the JIDT Gaussian solver between two 1D
    time series.
    """
    startJVM()
    te_calc = jpype.JPackage("infodynamics.measures.continuous.gaussian").TransferEntropyCalculatorMultiVariateGaussian()
    te_calc.setProperty("DELAY",str(tau))
    te_calc.setProperty("K",str(k))
    te_calc.setProperty("L",str(k))
    te_calc.initialise(1, 1)
    te_calc.setObservations(src.tolist(), tgt.tolist())
    te_val = te_calc.computeAverageLocalOfObservations()
    return te_val


def simulate(x, f):
    """
    Given a KalmanFilter f and an input time series x, returns a tuple (s,o)
    where s are the filtered states and o are the model prediction errors.
    """
    s = f.filter(x)[0]
    o = (f.observation_matrices @ s.T).T 
    o = np.array(o)
    return (s.squeeze()[1:], o.squeeze()[1:] - x[1:])


def LZ(x):
    """
    Convenience wrapper function to compute normalised LZ of a time series
    after detrending and binarisation.
    """
    return LZ76(1*((x-np.mean(x)) > 0))*np.log(len(x))/len(x)


def copy_filter(kf, prior_factor=1, likelihood_factor=1):
    f = KalmanFilter(n_dim_state=kf.n_dim_state, n_dim_obs=kf.n_dim_obs,
                     initial_state_mean=kf.initial_state_mean,
                     initial_state_covariance = prior_factor*kf.initial_state_covariance,
                     transition_covariance = prior_factor*kf.transition_covariance,
                     observation_covariance = likelihood_factor*kf.observation_covariance,
                     em_vars = ['transition_matrices','observation_matrices'])
    return f


def butter_bandpass(lowcut, highcut, fs, order=1):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=1, axis=0):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data, axis)
    return y


def load_data(eeg_source, sub_id, random_trials=True):
    '''
    Loads data from primary visual cortex (calcarine) for a given subject from the loaded MatLab dictionary.
    Data is of the shape (Time x AAL x Trials), where AAL is the number of AAL regions (24 for left calcarine, 25 for right calcarine).
    The function returns a numpy array of shape (Trials x Time) for the calcarine region.
    If random_trials is True, it randomly selects 10 trials from the available data.

    Parameters:
    - eeg_source: dict, loaded MatLab dictionary containing all subjects' data
    - sub_id: str, subject ID to extract data for
    - random_trials: bool, whether to select random trials or not
    Returns:
    - np.ndarray: shape (Trials x Time) for the calcarine region
    '''
    calcarine_idx = 24  # Index for the left calcarine region in AAL atlas (24 for left, 25 for right)
    normalise = lambda x: x/x.std()
    
    # Extract subject data: [samples * sources * trials]
    subject_data = eeg_source[sub_id]
    # Extract calcarine region data and transpose to [trials * samples]
    calcarine_data = subject_data[:, calcarine_idx, :].T
    
    # Select random trials and process data
    nb_trials = 10
    total_trials = calcarine_data.shape[0]
    
    if random_trials:
        idx = np.random.randint(total_trials, size=min(nb_trials, total_trials))
    else:
        idx = np.arange(min(nb_trials, total_trials))
    
    D = calcarine_data[idx]
    D_proc = []
    for i in range(len(D)): #Bandpass filtering for the data. Comment this out if data is already filtered
        filtered_data = butter_bandpass_filter(D[i,:], lowcut=1, highcut=100, fs=600, order=3)
        D_proc.append(normalise(filtered_data))
    
    return np.array(D_proc)  # Returns the processed data as a numpy array of shape (Trials x Time) for a given subject


def model_run():
    # Read the Excel file
    file_name = 'Data_4_Import_REST.xlsx'
    excel_sheet_name = 'Depression Rest'
    file_path = os.path.join("Depression_Study", "depression_data", file_name)
    df_excel = pd.read_excel(file_path, sheet_name=excel_sheet_name)
    df_excel['depressed'] = df_excel['MDD'].apply(lambda x: 1 if x <= 2 else 0)

    # Load all subjects' source reconstructed data from a matlab file
    data_file = os.path.join(project_root, 'eeg_source.mat')
    eeg_source = h5py.File(data_file, 'r')
    eeg_source_open = eeg_source['eeg_source_open']
    eeg_source_closed = eeg_source['eeg_source_closed'] # TODO - figure out how to use this data
    
    # Get all subject IDs from the struct fields
    subject_ids = [key for key in eeg_source_open.keys() if key[1:] in df_excel[df_excel['depressed'] == 0]['id'].astype(str).values]
    subject_ids = subject_ids[:3]

    df = []
    
    # Loop over each subject
    for z, sub_id in enumerate(tqdm(subject_ids)):
        # print(f"{z+1}/{len(subject_ids)}")
        
        # Loads filtered and normalised data of V1 calcarine region for the subject of shape (Trials x Time)
        D = load_data(eeg_source_open, sub_id, random_trials=True)
        
        ## Train baseline filter
        kf = KalmanFilter()
        for X in D:
            kf.em(X, n_iter=1)

        # Second round of training of baseline filter
        kf2 = copy_filter(kf)
        for X in D:
            kf2.em(X, n_iter=1)

        # Simulate predictions
        sim = [simulate(x, kf2) for x in D] # Get predicted states (frontal regions signal) and prediction errors (sensory regions signal)
        baseline_lz = np.mean([LZ(s[1]) for s in sim]) # LZ estimated for prediction errors
        # baseline_te = np.mean([TE(s[0], s[1], k=1) for s in sim]) # TE estimated from filtered states to prediction errors - front-to-back / top-down
        baseline_te = np.mean([TE(s[1], s[0], k=1) for s in sim]) # TE estimated from prediction errors to filtered states - back-to-front / bottom-up

        factor_vec = 2.0**np.arange(-5,6,1)
        
        for eta_prior in factor_vec:
            for eta_likelihood in factor_vec:
                # print(np.round(eta,2),end="...")

                # Try different prior and likelihood factors to match the LZ and TE changes in the Depressed condition.
                # Train one filter and then copy it with different prior/likelihood factors.
            
                ## Define model filter for depression
                depressed_kf = copy_filter(kf, prior_factor=eta_prior, likelihood_factor=eta_likelihood)
                # prior_factor=eta - increases prior variance - mimics overactive top-down priors
                # likelihood_factor=eta - increases likelihood variance - mimics overactive sensory trust (or weakened frontal control)

                ## Retrain the model with the same data, but with different prior/likelihood factors
                for X in D:
                    depressed_kf.em(X, n_iter=1)
            
                ## Simulate and compute LZ+TE
                for t, x in enumerate(D):
                    # LSD
                    depressed_sim = simulate(x, depressed_kf)
                    lz = LZ(depressed_sim[1])
                    te = TE(depressed_sim[0], depressed_sim[1], k=1)
                    df.append(pd.DataFrame({'Subject': sub_id, 'Dataset': 'EEG', 'Prior Factor': eta_prior, 'Likelihood Factor': eta_likelihood, 
                                            'Trial': t, 'LZ': lz, 'TE': te, 'Model': 'Depression'}, index=[0]))
            
                    # Add baseline again, for convenience
                    df.append(pd.DataFrame({'Subject': sub_id, 'Dataset': 'EEG', 'Prior Factor': eta_prior, 'Likelihood Factor': eta_likelihood,
                                            'Trial': t, 'LZ': baseline_lz, 'TE': baseline_te, 'Model': 'Baseline'}, index=[0]))
            # print("\n")
    
    return pd.concat(df, ignore_index=True)



if __name__ == '__main__':
    #sub_list = np.unique([int(fname.split('_')[0]) for fname in os.listdir(data_folder)])
    #pool = mp.Pool(processes=4)
    #res = pool.map(model_run, iterable= sub_list, chunksize=len(sub_list)//4)
    df = model_run()#pd.concat(res,ignore_index=True)

    df = df.melt(id_vars=['Subject', 'Dataset', 'Prior Factor', 'Likelihood Factor', 'Model', 'Trial'],#, 'Run'],
                 value_name='Value', var_name='Measure')
    
    df.to_csv("model_results_latest_proc10.csv")
    
    # Plot varying Prior Factor (with fixed Likelihood Factor = 1.0)
    df_prior = df[df['Likelihood Factor'] == 1.0]
    
    pl.figure(1)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Prior Factor', ci=60,
                   kind='line', facet_kws={'sharey': False}, data=df_prior)
    g.axes[0,0].set_xscale('log', base=2)
    g.figure.suptitle('Effect of Prior Factor (Likelihood Factor = 1.0)', y=1.02)
    pl.tight_layout()
    
    # Plot varying Likelihood Factor (with fixed Prior Factor = 1.0)
    df_likelihood = df[df['Prior Factor'] == 1.0]
    
    pl.figure(2)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Likelihood Factor', ci=60,
                   kind='line', facet_kws={'sharey': False}, data=df_likelihood)
    g.axes[0,0].set_xscale('log', base=2)
    g.figure.suptitle('Effect of Likelihood Factor (Prior Factor = 1.0)', y=1.02)
    pl.tight_layout()
    
    df_depression = df.loc[df['Model']=='Depression']
    df_depression['Model'] = 'Depression'
    
    df_base = df.loc[df['Model']=='Baseline']
    df_all = pd.concat([df_depression, df_base], ignore_index=True)
    df_all.to_csv('model_results_final_proc10.csv')
    
    # Combined plot with both models - Prior Factor
    df_all_prior = df_all[df_all['Likelihood Factor'] == 1.0]
    
    pl.figure(3)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Prior Factor', ci=60,
                   kind='line', facet_kws={'sharey': False}, data=df_all_prior)
    g.axes[0,0].set_xscale('log', base=2)
    g.figure.suptitle('Comparison: Prior Factor Effect (Likelihood Factor = 1.0)', y=1.02)
    pl.tight_layout()
    
    # Combined plot with both models - Likelihood Factor
    df_all_likelihood = df_all[df_all['Prior Factor'] == 1.0]
    
    pl.figure(4)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Likelihood Factor', ci=60,
                   kind='line', facet_kws={'sharey': False}, data=df_all_likelihood)
    g.axes[0,0].set_xscale('log', base=2)
    g.figure.suptitle('Comparison: Likelihood Factor Effect (Prior Factor = 1.0)', y=1.02)
    pl.tight_layout()
    
    pl.show()

    ##############

    # Filter only the depression model (exclude baseline if needed)
    df_dep = df[df['Model'] == 'Depression']

    # Compute mean across trials and subjects
    agg = df_dep.groupby(['PriorFactor', 'LikelihoodFactor']).agg({'TE': 'mean', 'LZ': 'mean'}).reset_index()

    # Pivot for heatmaps
    te_pivot = agg.pivot(index='LikelihoodFactor', columns='PriorFactor', values='TE')
    lz_pivot = agg.pivot(index='LikelihoodFactor', columns='PriorFactor', values='LZ')

    pl.figure(figsize=(8,6))
    sn.heatmap(te_pivot, cmap='viridis', annot=True, fmt=".3f", cbar_kws={'label': 'Transfer Entropy'})
    pl.title("Transfer Entropy (Prediction Errors → Hidden States)")
    pl.xlabel("Prior Factor")
    pl.ylabel("Likelihood Factor")
    pl.tight_layout()
    pl.show()

    pl.figure(figsize=(8,6))
    sn.heatmap(lz_pivot, cmap='magma', annot=True, fmt=".3f", cbar_kws={'label': 'Lempel-Ziv Complexity'})
    pl.title("LZ Complexity of Prediction Errors")
    pl.xlabel("Prior Factor")
    pl.ylabel("Likelihood Factor")
    pl.tight_layout()
    pl.show()
