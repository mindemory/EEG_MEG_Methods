import numpy as np
import os, mne
from datetime import datetime
from mne_bids import write_raw_bids, BIDSPath
from mne.channels import make_dig_montage
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the data
rootPath = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/EEG'
bids_root = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/EEGBids'
groupName = 'GroupD'
recRoot = os.path.join(rootPath, groupName, 'Recording')
timeStampsRoot = os.path.join(rootPath, groupName, 'TimeStamps')

# Load electrode coordinates
elecTsv = os.path.join(recRoot, 'electrodeCoords.tsv')
electrodes_df = pd.read_csv(elecTsv, sep='\t')
chNames = electrodes_df['name'].tolist()
chCoords = electrodes_df[['x', 'y', 'z']].values
chPos = dict(zip(chNames, chCoords))

# Load fiducials and headshape from a .pos file
posFile = os.path.join(recRoot, 'fiducials_and_headshape.tsv')
fiducials_and_headshape = pd.read_csv(posFile, sep='\t', header=None, names=['Type', 'X', 'Y', 'Z'],  keep_default_na=False)

fiducials = {}
for _, row in fiducials_and_headshape.iterrows():
    if row['Type'] in ['LPA', 'RPA', 'NA']:
        fiducials[row['Type']] = [row['X'], row['Y'], row['Z']]
headshape_points = fiducials_and_headshape[fiducials_and_headshape['Type'] == 'EXTRA'][['X', 'Y', 'Z']].astype(float).values

dig_ch_pos = {
    'LPA': fiducials['LPA'], 
    'RPA': fiducials['RPA'],
    'NAS': fiducials['NA'],
}

montage = make_dig_montage(ch_pos=chPos,
                            coord_frame='head',
                            nasion=fiducials['NA'],
                            lpa=fiducials['LPA'],
                            rpa=fiducials['RPA'],
                            hsp=headshape_points)
# Visualize the montage
# montage.plot(kind='3d')
# plt.show()

# Get the list of all files in recRoot
recFiles = os.listdir(recRoot)

bdfFiles = [f for f in recFiles if f.endswith('.bdf')]
raw = mne.io.read_raw_bdf(os.path.join(recRoot, bdfFiles[0]), preload=False)
raw.info['line_freq'] = 60 # Set the line frequency to 60 Hz

# Remove the last 4 EXG channels
chansToRemove = ['EXG5', 'EXG6', 'EXG7', 'EXG8', 'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp']
raw.drop_channels(chansToRemove)
non_eeg_channels = {
    'EXG1': 'eog',
    'EXG2': 'eog',
    'EXG3': 'eog',
    'EXG4': 'eog',
}
raw.set_channel_types(non_eeg_channels)
raw.set_montage(montage, on_missing='ignore')


def fix_missing_flags(events_df):
    ''' 
    ome of the runs have missing flags for the start of the block. This code fixes the missing flags 
    by inserting a new row before the missing flag and updating the blockNumber and trialNumber accordingly.
    It also fixes the trialNumber for the next 30 trials after the missing flag.
    '''
    blkIdx = 1
    while True:
        max_block_number = events_df['blockNumber'].max()
        if blkIdx >= max_block_number:
            break
        trlCount = events_df[events_df['blockNumber'] == blkIdx]['trialNumber'].max()
        blkTypeThis = events_df[events_df['blockNumber'] == blkIdx]['blockType'].values[0]

        if blkTypeThis != 'storyStart':
            if trlCount > 30:
                # Find index of trial 31 in the block
                trl31Idx = events_df[(events_df['blockNumber'] == blkIdx) & (events_df['trialNumber'] == 31)].index[0]
                # Find index of first trial after where it is either 0 or 1
                trlNextIdx = events_df[(events_df['blockNumber'] == blkIdx+1) & ((events_df['trialNumber'] == 0) | (events_df['trialNumber'] == 1))].index[0]
                
                # Add blocknumber+1 for all blocks starting at trl31Idx
                events_df.loc[trl31Idx:, 'blockNumber'] += 1
                # Change trialNumber by subtracting 30 from trl31Idx to next 30 trials
                events_df.loc[trl31Idx:trlNextIdx-1, 'trialNumber'] -= 30

                # Insert a new row before trl31Idx
                sample_onset_avg = int((events_df.loc[trl31Idx, 'sampleOnset'] + events_df.loc[trl31Idx-1, 'sampleOnset']) / 2)
                time_onset_avg = (events_df.onset[trl31Idx] + events_df.onset[trl31Idx-1]) / 2
                last_trigger_code = events_df[(events_df['blockNumber'] <= blkIdx) & (events_df['trialNumber'] == 0)]['value'].values[-1]
                
                # Get trialOverBlocks from previous row and blockType from trlIdx row
                trial_over_blocks_prev = events_df.loc[trl31Idx-1, 'trialoverBlocks']
                # block_type_current = events_df.loc[trl31Idx, 'blockType']
                block_type_current = 'Unsure'

                new_row = {
                    'sampleOnset': sample_onset_avg,
                    'onset': time_onset_avg,
                    'duration': 0.0,
                    'value': last_trigger_code,
                    'blockNumber': blkIdx + 1,
                    'trialNumber': 0,
                    'trialoverBlocks': trial_over_blocks_prev,
                    'trial_type': None,
                    'blockType': block_type_current
                }
                events_df = pd.concat([events_df.iloc[:trl31Idx], pd.DataFrame([new_row]), events_df.iloc[trl31Idx:]]).reset_index(drop=True)
        else:
            if trlCount > 8:
                # Find index of trial 9 in the block
                trl9Idx = events_df[(events_df['blockNumber'] == blkIdx) & (events_df['trialNumber'] == 9)].index[0]
                # Find index of first trial after where it is either 0 or 1
                trlNextIdx = events_df[(events_df['blockNumber'] == blkIdx+1) & ((events_df['trialNumber'] == 0) | (events_df['trialNumber'] == 1))].index[0]
                
                # Add blocknumber+1 for all blocks starting at trl9Idx
                events_df.loc[trl9Idx:, 'blockNumber'] += 1
                # Change trialNumber by subtracting 8 from trl9Idx to next 8 trials
                events_df.loc[trl9Idx:trlNextIdx-1, 'trialNumber'] -= 8

                # Insert a new row before trl9Idx
                sample_onset_avg = int((events_df.loc[trl9Idx, 'sampleOnset'] + events_df.loc[trl9Idx-1, 'sampleOnset']) / 2)
                time_onset_avg = (events_df.onset[trl9Idx] + events_df.onset[trl9Idx-1]) / 2
                last_trigger_code = events_df[(events_df['blockNumber'] <= blkIdx) & (events_df['trialNumber'] == 0)]['value'].values[-1]
                
                # Get trialOverBlocks from previous row and blockType from trlIdx row
                trial_over_blocks_prev = events_df.loc[trl9Idx-1, 'trialoverBlocks']
                # block_type_current = events_df.loc[trl9Idx, 'blockType']
                block_type_current = 'Unsure'

                new_row = {
                    'sampleOnset': sample_onset_avg,
                    'onset': time_onset_avg,
                    'duration': 0.0,
                    'value': last_trigger_code,
                    'blockNumber': blkIdx + 1,
                    'trialNumber': 0,
                    'trialoverBlocks': trial_over_blocks_prev,
                    'trial_type': None,
                    'blockType': block_type_current
                }
                events_df = pd.concat([events_df.iloc[:trl9Idx], pd.DataFrame([new_row]), events_df.iloc[trl9Idx:]])

        blkIdx += 1
    return events_df

def addStimulusMeta(events_df, timeStampsRoot):
    timeStampsFiles = os.listdir(timeStampsRoot)
    timeStampsFiles = [f for f in timeStampsFiles if f.endswith('.mat')]
    timeStampsFiles.sort()  # Sort files in ascending order of time
    categories = [filename.split('_')[-1].split('.')[0] for filename in timeStampsFiles]
    blkTypeActual = [filename.split('_')[-2].split('_')[0] for filename in timeStampsFiles]
    events_df['stimulusCategory'] = None
    events_df['stimulusName'] = None
    events_df['stimulusDuration'] = None
    events_df['blockTypeActual'] = None
    events_df['value_fixed'] = events_df['value']

    maxBlocks = events_df['blockNumber'].max()
    for i, blkIdx in enumerate(range(1, maxBlocks+1)):
        print(f'Processing block {blkIdx}')
        trlIdx = events_df[events_df['blockNumber'] == blkIdx].index
        events_df.loc[trlIdx, 'stimulusCategory'] = categories[i]
        if blkTypeActual[i] == 'saud':
            events_df.loc[trlIdx, 'blockTypeActual'] = 'semanticAud'
            events_df.loc[trlIdx[0], 'value_fixed'] = 4
        elif blkTypeActual[i] == 'svis':
            events_df.loc[trlIdx, 'blockTypeActual'] = 'semanticVis'
            events_df.loc[trlIdx[0], 'value_fixed'] = 2
        elif blkTypeActual[i] == 'caud':
            events_df.loc[trlIdx, 'blockTypeActual'] = 'classicalAud'
            events_df.loc[trlIdx[0], 'value_fixed'] = 1
        elif blkTypeActual[i] == 'timingData':
            events_df.loc[trlIdx, 'blockTypeActual'] = 'story'
            events_df.loc[trlIdx[0], 'value_fixed'] = 8

        # Load the timeStamps file
        timeStampsFile = os.path.join(timeStampsRoot, timeStampsFiles[i])
        timeStamps = loadmat(timeStampsFile)
        if categories[i] != 'story2':
            timingData = timeStamps['timingData']

            # Get stimulus names
            sti_names = [None] + [entry['stiName'][0] for entry in timingData[0]]
            events_df.loc[trlIdx, 'stimulusName'] = sti_names

            # Get stimulus duration
            sti_durations = [None] + [round(entry['offsetTime'][0][0] - entry['onsetTime'][0][0], 4) for entry in timingData[0]]
            events_df.loc[trlIdx, 'stimulusDuration'] = sti_durations

        else:
            timingData = timeStamps['timingData3']

            # Get event names
            sti_names = [None] + [entry['event'][0] for entry in timingData[0]]
            events_df.loc[trlIdx, 'stimulusName'] = sti_names

            # Get event durations
            sti_durations = [None] + [round(entry['offsetTime'][0][0] - entry['onsetTime'][0][0], 4) for entry in timingData[0]]
            events_df.loc[trlIdx, 'stimulusDuration'] = sti_durations

    return events_df


# Extract events from STATUS channel
status_channel = 'Status'
if status_channel in raw.ch_names:
    status_data = raw.get_data(picks=[status_channel])[0]
    times = raw.times
    
    event_samples = []
    event_times = []
    event_values = []
    trl_types = []
    block_types = []
    block_numbers = []
    trial_numbers_within_block = []
    trial_numbers_over_blocks = []

    current_block_type = None
    current_block_number = 0
    current_trial_number = 0
    global_trial_number = 0

    for i in range(1, len(status_data)):
        if status_data[i] != status_data[i-1]:
            event_samples.append(i)
            event_times.append(times[i])
            event_values.append(status_data[i])

            if status_data[i] == 64:
                trl_types.append('even')
            elif status_data[i] == 128:
                trl_types.append('odd')
            elif status_data[i] == 16:
                trl_types.append('storySeg')
            else:
                trl_types.append(None)

            # Determine block type based on trigger code
            if status_data[i] in [1, 2, 4, 8]:
                if status_data[i] == 1:
                    current_block_type = 'classicalAud'
                elif status_data[i] == 2:
                    current_block_type = 'semanticVis'
                elif status_data[i] == 4:
                    current_block_type = 'semanticAud'
                elif status_data[i] == 8:
                    current_block_type = 'storyStart'

                current_block_number += 1
                current_trial_number = 0

            elif status_data[i] in [16, 64, 128]:
                current_trial_number += 1
                global_trial_number += 1

            # Append the current block type to the block_types list
            block_types.append(current_block_type)
            block_numbers.append(current_block_number)
            trial_numbers_within_block.append(current_trial_number)
            trial_numbers_over_blocks.append(global_trial_number)
            # checker.append(0)

    events_df = pd.DataFrame({
        'sampleOnset': event_samples,
        'onset': event_times,
        'duration': [0.0] * len(event_times),
        'value': event_values,
        'blockNumber': block_numbers,
        'trialNumber': trial_numbers_within_block,
        'trialoverBlocks': trial_numbers_over_blocks,
        'trial_type': trl_types,
        'blockType': block_types,
    })

    # Set column trigger_codes to be of type int
    events_df['value'] = events_df['value'].astype(int)
    valid_trigger_codes = [1, 2, 4, 8, 16, 64, 128]
    events_df = events_df[events_df['value'].isin(valid_trigger_codes)]

    # Reset index
    events_df.reset_index(drop=True, inplace=True)

    # Fix missing flags
    events_df = fix_missing_flags(events_df)
    print(events_df['blockNumber'].max()) 

    # Add stimulus meta data
    events_df = addStimulusMeta(events_df, timeStampsRoot)
    

# Set the subject ID and task name
if groupName[-1] == 'A':
    subject_id = '001'
elif groupName[-1] == 'B':
    subject_id = '002'
elif groupName[-1] == 'C':
    subject_id = '003'
elif groupName[-1] == 'D':
    subject_id = '004'

taskName = 'oddball'
 
# Write the raw data to BIDS
bids_path = BIDSPath(subject=subject_id, task=taskName, root=bids_root)
write_raw_bids(raw=raw, bids_path=bids_path, overwrite=True, verbose=True)
# Write the events to a TSV file
events_tsv_path = os.path.join(bids_root, 'sub-' + subject_id, 'eeg', 'sub-' + subject_id + '_task-' + taskName + '_events.tsv')
os.makedirs(os.path.dirname(events_tsv_path), exist_ok=True)
events_df.to_csv(events_tsv_path, sep='\t', index=False)