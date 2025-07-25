"""
Train and test minimal predictive processing model to reproduce LZ and TE changes in
the psychedelic state and schizophrenia patients. The model consists of a Kalman filter
trained on time series from the placebo condition of the LSD_MEG dataset, and LSD/SCH
are modelled as increased prior/likelihood variance, respectively.

Pedro Mediano and Hardik Rajpal, Apr 2021
"""
from scipy.io import loadmat
from pykalman import KalmanFilter
import numpy as np
import matplotlib.pyplot as pl
from lz76 import LZ76
#from glob import glob
from scipy.signal import decimate
import pandas as pd
import seaborn as sn
import multiprocessing as mp
import jpype
import os
from scipy.signal import butter, filtfilt
from tqdm import tqdm
from itertools import product
import sys

data_folder = 'F:/LSD_Data/Data1/'#'/path/to/LSD/data/'

def startJVM():
    """
    Start JVM using jpype, if not already running. Assumes the JIDT jar is in
    the current folder.
    """
    jarLocation = os.path.join(os.getcwd(),"infodynamics.jar")
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


def load_data(data_path,fname,random_trials=True):
    '''
    Loads data from primary visual cortex (calcarine) for a given subject and condition, specified by data_path and fname.
    Data is of the shape (AAL x Trials x Time), where AAL is the number of AAL regions (42 for left calcarine, 43 for right calcarine).
    The function returns a numpy array of shape (Trials x Time) for the calcarine region.
    If random_trials is True, it randomly selects 10 trials from the available data.

    Parameters:
    - data_path: str, path to the data folder
    - fname: str, name of the data file
    - random_trials: bool, whether to select random trials or not
    Returns:
    - np.ndarray: shape (Trials x Time) for the calcarine region
    '''
    # 42 for left calcarine
    # 43 for right calcarine

    calcarine_idx = 42 # Choose an index based on the AAL atlas used corresponding to calcarine or primary visual cortex

    normalise = lambda x: x/x.std()
    fpath = os.path.join(data_path,fname)
    #Keep in mind the shape of the loaded data
    #My version has AAL x Trials x Time
    if 'LSD' in data_path: # Replace with your actual data path for Depression/Healthy controls
        D = np.squeeze(loadmat(fpath)['timeseries'][0][calcarine_idx])
    if 'KET' in data_path:
        D = np.stack(loadmat(fpath)['dat'][0])[:,calcarine_idx]
    
    
    nb_trials = 10 #Sampling random trials from the available data 
    if random_trials:
        idx = np.random.randint(len(D), size=nb_trials)
    else:
        idx = np.arange(nb_trials)
    
    D = D[idx]
    D_proc = []
    for i in range(len(D)): #Bandpass filtering for the data. Comment this out if data is already filtered
       D_proc.append(normalise(butter_bandpass_filter(D[i,:], lowcut=1, highcut=100, fs=600,order=3)))
    
    
    return np.array(D_proc) #Returns the processed data as a numpy array of shape (Trials x Time) for a given subject
    


def model_run():

    #Collecting filenames from the data folders to train the baseline kalman filter on
    #Replace with your actual data path and select filenames for the Healthy Controls only.
    data_folder1 = 'F:/LSD_Data/Data1/'
    data_folder2 = 'F:/KET_Data/KET_Data/'
    
    #sub = arg[0]
    #eta = arg[1]
    #Use sub_list = [2] for the best example
    df = []
    flist1 = [i for i in os.listdir(data_folder1) if i.split('_')[1]=='PLA']
    flist2 = [i for i in os.listdir(data_folder2) if i.split('_')[1]=='PLA']
    datapath_list1 = [data_folder1 for i in range(len(flist1))]
    datapath_list2 = [data_folder2 for i in range(len(flist2))]
    flist = flist1 + flist2
    dfolder_list = datapath_list1 + datapath_list2
    #sub_list = np.unique([int(fname.split('_')[0]) for fname in os.listdir(data_folder)])
    #Looping over each subject
    z=0
    for fname,dfolder in tqdm(zip(flist,dfolder_list)):
        z+=1
        print(f"{z}/{len(flist)}")
        sub = fname.split('_')[0]
        if 'LSD' in dfolder:
            drug = 'LSD'
        else:
            drug = 'KET'
        
        #Load data for the subject
        D = load_data(dfolder, fname, random_trials= True)
        
        ## Train baseline filter
        kf = KalmanFilter()
        
        for X in D:
            kf.em(X, n_iter=1)
    
    
        # Second round of training of baseline filter, to match the LSD/SCH experimental conditions
        kf2 = copy_filter(kf)
        for X in D:
            kf2.em(X, n_iter=1)
    
        # Simulate predictions
        sim = [simulate(x, kf2) for x in D]
        baseline_lz = np.mean([LZ(s[1]) for s in sim]) #LZ estimated for prediction errors (stays same for Depression Dataset)
        baseline_te = np.mean([TE(s[0], s[1], k=1) for s in sim]) # TE estimated from filtered states to prediction errors (opposite direction for Depression dataset)
    
        factor_vec = 2.0**np.arange(-5,6,1)#
        
        for eta in factor_vec:
            print(np.round(eta,2),end="...")
            #sys.stdout.flush()
        
            ## Define model filters for LSD and SCH
            lsd_kf = copy_filter(kf, prior_factor=eta) # We increase the prior variance for LSD
            sch_kf = copy_filter(kf, likelihood_factor=eta) # We increase the likelihood variance for SCH
            # We will need to try different prior and likelihood factors to match the LZ and TE changes in the Depressed condition.
            # We only have 1 condition for the Depression dataset, so we just train one filter and then copy it with different prior/likelihood factors.
        
            ## Re-train models
            for X in D:
                # We now retrain the models with the same data, but with different prior/likelihood factors
                lsd_kf.em(X, n_iter=1)
                sch_kf.em(X, n_iter=1)
        
            ## Simulate and compute LZ+TE
            for t,x in enumerate(D):
                # LSD
                lsd_sim = simulate(x, lsd_kf)
                lz = LZ(lsd_sim[1])
                te = TE(lsd_sim[0], lsd_sim[1], k=1)
                #Renaming this as the prior model
                df.append(pd.DataFrame({'Subject': sub, 'Dataset':drug, 'Factor': eta, 'Trial': t,
                                        'LZ': lz, 'TE': te, 'Model': 'Prior'}, index=[0]))
        
                # SCH
                sch_sim = simulate(x, sch_kf)
                lz = LZ(sch_sim[1])
                te = TE(sch_sim[0], sch_sim[1], k=1)
                #Renaming this as the likelihood model
                df.append(pd.DataFrame({'Subject': sub, 'Dataset':drug, 'Factor': eta, 'Trial': t,
                                        'LZ': lz, 'TE': te, 'Model': 'Likelihood'}, index=[0]))
        
                # Add baseline again, for convenience
                df.append(pd.DataFrame({'Subject': sub, 'Dataset':drug, 'Factor': eta, 'Trial': t,
                                        'LZ': baseline_lz, 'TE': baseline_te, 'Model': 'Baseline'}, index=[0]))
        print("\n")
    
    return pd.concat(df, ignore_index=True)



if __name__ == '__main__':
    #sub_list = np.unique([int(fname.split('_')[0]) for fname in os.listdir(data_folder)])
    #pool = mp.Pool(processes=4)
    #res = pool.map(model_run, iterable= sub_list, chunksize=len(sub_list)//4)
    df = model_run()#pd.concat(res,ignore_index=True)

    #df = model_run()
    df = df.melt(id_vars=['Subject','Dataset','Factor', 'Model', 'Trial'],#, 'Run'],
                 value_name='Value', var_name='Measure')
    
    df.to_csv("model_results_latest_proc10.csv")
    
    pl.figure(1)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Factor', ci=60,
               kind='line', facet_kws={'sharey': False}, data=df)
    g.axes[0,0].set_xscale('log',basex=2)
    pl.tight_layout()
    
    df_schiz = df.loc[df['Model']=='Likelihood'].loc[df['Factor']<=1]
    df_schiz['Factor'] = 1/df_schiz['Factor'].values
    df_schiz['Model'] = 'Schizophrenia'
    
    df_lsd = df.loc[df['Model']=='Prior'].loc[df['Factor']>=1]
    df_lsd['Model'] = 'Psychedelics'
    
    df_base = df.loc[df['Model']=='Baseline'].loc[df['Factor']>=1]
    df_all = df_schiz.append(df_lsd,ignore_index=True)
    df_all = df_all.append(df_base,ignore_index=True)
    df_all.to_csv('model_results_final_proc10.csv')
    
    pl.figure(2)
    g = sn.relplot(col='Measure', y='Value', hue='Model', x='Factor', ci=60,
               kind='line', facet_kws={'sharey': False}, data=df_all)
    g.axes[0,0].set_xscale('log',basex=2)
    pl.tight_layout()
    pl.show()


