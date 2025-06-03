import mne
import matplotlib.pyplot as plt
import numpy as np
import os

from autoreject import AutoReject
from mne.preprocessing import ICA
from mne_icalabel import label_components
from pyprep.prep_pipeline import PrepPipeline
from scipy.io import loadmat, savemat

mne.set_log_level("WARNING")


### LOAD .MAT FILE > PROCESS CHANNELS & SIGNALS > EYE EPOCHS

def load_data(file_name):
    data = loadmat(file_name)["EEG"][0][0]
    signals = data[15]/1e6
    channels = data[21][0]
    
    # Extract channel information
    channel_names = np.array([channel[0][0] for channel in channels])
    channels_to_keep = [i for i, name in enumerate(channel_names) if name not in {"CB1","CB2","EKG"}]
    channel_names = channel_names[channels_to_keep]
    channel_types = ["eeg" if name not in {"VEOG","HEOG"} else "eog" for name in channel_names]
    channel_locs = [[channel[j][0][0] for j in range(4,7)] for channel in channels if channel[0][0] not in {"CB1", "CB2", "EKG"}]

    # Scale coordinates cm->m and rotate them
    channel_locs = np.array(channel_locs) / 1000
    channel_locs = np.stack([
        channel_locs[:, 1],
        channel_locs[:, 0],
        channel_locs[:, 2]
    ], axis=1)

    # Process events information
    def extract_event_value(event): # string and integer events are stored differently
        return event[0] if isinstance(event[0], np.str_) else event[0][0]

    events = np.array([[extract_event_value(event[j]) for j in range(2)] for event in data[25][0]])
    events = [[int(event[0]), int(event[1])] for event in events if str(event[0]).isdigit()] # drop string events
    events = events[1:len(events)-1] # remove first and last event 17    
    
    # Parse events into key timestamps
    event_labels = [event[0] % 2 for event in events] # open or closed eyes
    switching_indeces = [i+1 for i in range(len(event_labels)-1) if event_labels[i] != event_labels[i+1]]
    switching_tstamps = [events[0][1]] + [events[i][1] for i in switching_indeces] + [events[-1][1]]

    # Slice epochs based on initial labels
    signals = signals[channels_to_keep]
    epochs = [signals[:, switching_tstamps[i]:switching_tstamps[i+1]] for i in range(len(switching_tstamps) - 1)]

    if event_labels[0] == 0:
        epochs_open = epochs[::2]
        epochs_closed = epochs[1::2]
    else:
        epochs_closed = epochs[::2]
        epochs_open = epochs[1::2]

    return (epochs_open, epochs_closed, channel_names, channel_locs, channel_types)


### DEFINE INFO & MONTAGE

def setup_info(channel_names, channel_locs, channel_types):

    # Create montage
    montage = mne.channels.make_dig_montage(
        ch_pos      = dict(zip(channel_names, channel_locs)),
        coord_frame = 'head')

    # Create MNE Info object
    info = mne.create_info(
        ch_names = channel_names.tolist(),
        sfreq    = 500,
        ch_types = channel_types)
    info.set_montage(montage, match_alias={'VEOG':'eog','HEOG':'eog'})
    print(info)
    print()

    return info


### SCALING > HP FILTER > PREP

def initial_processing(epoch, info):

    # Apply high-pass filtering
    raw_eeg = mne.io.RawArray(epoch, info)
    hp_eeg = raw_eeg.copy().filter(l_freq=1, h_freq=100)

    # Apply PREP pipeline
    sfreq = info["sfreq"]
    prep_params = {
        "ref_chs": "eeg",
        "reref_chs": "eeg",
        "line_freqs": []
    }

    print("Starting PREP...\n")

    prep = PrepPipeline(hp_eeg, prep_params, info.get_montage(), random_state=0)
    prep.fit()
    hp_prep_eeg = prep.raw

    print("\nInitial processing complete")
    print(f" > channels fixed after interpolation: {prep.interpolated_channels}")
    print(f" > channels still noisy after interpolation: {prep.still_noisy_channels}")

    return hp_prep_eeg


### EPOCH CREATION & REJECTION

def epochize_and_filter(eeg):

    epochs = mne.make_fixed_length_epochs(eeg, duration=4.0, preload=True, reject_by_annotation=False, proj=False)
    ar = AutoReject(n_interpolate=[1,2,4,8], random_state=0) # up to 15% (high fidelity) of 64 channels
    ar_epochs, reject_log = ar.fit_transform(epochs, return_log=True)

    # ar_epochs contains all the epochs, with bad ones marked
    # reject_log contains the bad epochs

    return (epochs, ar_epochs, reject_log)


### ICA

def perform_ica(epochs, ica_method='fastica'):
    
    ica = ICA(n_components=0.9999999, random_state=0, method=ica_method, max_iter='auto', 
              fit_params=None if ica_method == 'fastica' else dict(extended=True))
    ica.fit(epochs, reject_by_annotation=True)
    
    print()
    print(ica)

    # TEST EOG
    eog_indices, eog_scores = ica.find_bads_eog(epochs, measure='correlation', threshold=0.6)
    print(f"\nExcluded EOG components: {eog_indices}")
    print(f"EOG scores: {eog_scores}\n")

    # # TEST ECG
    # ecg_indices, ecg_scores = ica.find_bads_ecg(epochs)
    # print(f"\nExcluded ECG components: {ecg_indices}")
    # print(f"ECG scores: {ecg_scores}\n")

    # TEST MUSCLE
    muscle_indices, muscle_scores = ica.find_bads_muscle(epochs)
    print(f"\nExcluded muscle components: {muscle_indices}")
    print(f"Muscle scores: {muscle_scores}\n")

    # Label ICA components
    icalab = label_components(epochs, ica, method='iclabel')

    # “Other” is a catch-all that for non-classifiable components. 
    # We will stay on the side of caution and assume we cannot blindly remove these.
    for idx, label in enumerate(icalab['labels']):
        if label not in {'brain', 'other'}:
            print(idx, icalab['y_pred_proba'][idx], "\t", label)
            ica.exclude.append(idx)

    total_count = len(icalab['labels'])
    brain_other_count = total_count - len(ica.exclude)
    print(f"\nBrain & other components ratio: {brain_other_count}/{total_count}")
    print(f"ICALabel exclusions: {ica.exclude}\n")

    # Remove ICA component labeled as non-brain
    ica_epochs = epochs.copy()
    ica.apply(ica_epochs) # modifies in-place

    return (ica, ica_epochs)


### SAVE DATA AS FIF

def save_fif(file_name, epochs_open, epochs_closed):

    # Create the directory if it doesn't exist
    save_dir = "export_test_fif"
    os.makedirs(save_dir, exist_ok=True)

    if epochs_open:
        epochs_open_concat = mne.concatenate_epochs(epochs_open)
        file_path = os.path.join(save_dir, f"{file_name}_open_concat.fif")
        epochs_open_concat.save(file_path, overwrite=True)

    if epochs_closed:
        epochs_closed_concat = mne.concatenate_epochs(epochs_closed)
        file_path = os.path.join(save_dir, f"{file_name}_closed_concat.fif")
        epochs_closed_concat.save(file_path, overwrite=True)


### SAVE DATA AS MAT

def save_mat(file_name, info, epochs_open, epochs_closed):

    # Create the directory if it doesn't exist
    save_dir = "export_mat"
    os.makedirs(save_dir, exist_ok=True)

    epochs_open_data = mne.concatenate_epochs(epochs_open).get_data() if epochs_open else []
    epochs_closed_data = mne.concatenate_epochs(epochs_closed).get_data() if epochs_closed else []
    channel_names_data = np.array(info["ch_names"], dtype=object).reshape(-1, 1)

    data = {
        'epochs_open': epochs_open_data,
        'epochs_closed': epochs_closed_data,
        'channel_names': channel_names_data
    }

    file_path = os.path.join(save_dir, f"{file_name}_processed.mat")
    # savemat(file_path, {'data': data})
    savemat(file_path, data)


# LOAD DATA, CREATE MONTAGE, DEVELOP

def process_file(file_name):

    # Load data
    folder_name = "data"
    file_path = os.path.join(folder_name, f"{file_name}.mat")
    epochs_open, epochs_closed, channel_names, channel_locs, channel_types = load_data(file_path)
    print(f"File {file_name}.mat loaded")

    open_count = len(epochs_open)
    closed_count = len(epochs_closed)
    total_count = open_count + closed_count
    print(f"{total_count} epochs ({len(epochs_open)} open & {len(epochs_closed)} closed)\n")

    info = setup_info(channel_names, channel_locs, channel_types)

    # IPNYB Access
    ipnyb_data = {}

    # Process EEG
    def process_eye_epochs(eye_epochs, type, count):
        processed_eye_epochs = []
        for idx, epoch in enumerate(eye_epochs):
            if len(epoch[0]) <= 303:
                print("\n***** SKIPPING EPOCH - shorter than 304 *****\n\n")
                continue
            
            try:
                hp_prep_eeg = initial_processing(epoch, info)
                (fixed_length_epochs, ar_epochs, reject_log) = epochize_and_filter(hp_prep_eeg)
                (ica, ica_epochs) = perform_ica(ar_epochs, ica_method='infomax')
                # Resetting bad channels is also important for concatenation
                ica_epochs.interpolate_bads(reset_bads=True, method=dict(eeg="spline"))
            except Exception as e:
                print(f"\n***** FAILED {type} epoch {idx + 1}/{count} (total: {total_count}) *****")
                print(f"{e}\n\n")
                continue
            
            processed_eye_epochs.append(ica_epochs)
            print(f"\n***** PROCESSED {type} epoch {idx + 1}/{count} (total: {total_count}) *****\n\n")

            ipnyb_local = {}

            ipnyb_local["hp_prep_eeg"] = hp_prep_eeg
            ipnyb_local["fixed_length_epochs"] = fixed_length_epochs
            ipnyb_local["ar_epochs"] = ar_epochs
            ipnyb_local["reject_log"] = reject_log
            ipnyb_local["ica"] = ica
            ipnyb_local["ica_epochs"] = ica_epochs

            ipnyb_data[f"{type}_{idx}"] = ipnyb_local
        
        return processed_eye_epochs
    
    processed_epochs_open = process_eye_epochs(epochs_open, "open", open_count)
    processed_epochs_closed = process_eye_epochs(epochs_closed, "closed", closed_count)

    ipnyb_data['file_name'] = file_name
    ipnyb_data['info'] = info
    ipnyb_data['processed_epochs_open'] = processed_epochs_open
    ipnyb_data['processed_epochs_closed'] = processed_epochs_closed

    # Save data
    try:
        save_fif(file_name, processed_epochs_open, processed_epochs_closed)
        print(f"FIF export SUCEEDED - {file_name}")
    except Exception as e:
        print(f"FIF export FAILED - {file_name}")
        print(e)

    try:
        save_mat(file_name, info, processed_epochs_open, processed_epochs_closed)
        print(f"MAT export SUCEEDED - {file_name}")
    except Exception as e:
        print(f"MAT export FAILED - {file_name}")
        print(e)

    return ipnyb_data


healthy_sample = ["509", "517", "519", "523", "533"]
depressed_sample = ["559", "561", "567", "587", "624"]
sample = healthy_sample + depressed_sample
# sample = ["567"]
sample_count = len(sample)
ipnyb = {}

for idx, unit in enumerate(sample):
    file_name = unit + "_Depression_REST"
    print(f"\n##### {file_name} - START PROCESSING - {idx + 1}/{sample_count} #####\n\n")
    
    try:
        ipnyb[unit] = process_file(file_name)
    except Exception as e:
        print(f"\n\n##### {file_name} - PROCESSING FAILED - {idx + 1}/{sample_count} #####")
        print(f"{e}\n")
        continue
    
    print(f"\n\n##### {file_name} - PROCESSING COMPLETED - {idx + 1}/{sample_count} #####\n")
