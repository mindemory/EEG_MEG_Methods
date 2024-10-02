'''
author: @amilleah
updated for the course: @mrugank.dake@nyu.edu (October 1st, 2024)
'''

import os, glob
import mne
from numpy import array
import subprocess
from datetime import datetime
# from mne_bids import write_raw_bids, BIDSPath
# some commented out parts could be modified to use with BIDS

#pip install mne-kit-gui

event_id = {'rsvp':164, 'rpvp':165} # use event ids from experiment

def reformat_date(date_str): #coordinate .sqd file naming convention with mrk file naming convention
    """Convert date from MM.DD.YY to YYMMDD."""
    date_obj = datetime.strptime(date_str, "%m.%d.%y")
    return date_obj.strftime("%y%m%d")

def process_folder(folder_path):
    folder_name = os.path.basename(folder_path)
    expt = 'PP' # file naming convention, your experiment name abbreviated: "prickly pear"
    if folder_name.startswith('R') and expt in folder_name: # and folder_name in current_left:
        parts = folder_name.split('_')
        subject_id = parts[0]
        date_str = parts[2]
        formatted_date = reformat_date(date_str)
    
        for filename in os.listdir(folder_path):
            if filename.endswith('NR.sqd') and 'emptyroom' not in filename.lower(): #process only noise reduced signal
                sqd_file = filename
                # beginning and end of recording markers, use if else statement for more blocks
                n1 = 1
                n2 = 2

                # Construct filenames based on the NR.sqd file
                base_name = os.path.splitext(sqd_file)[0]  # Remove extension
                mrk1 = os.path.join(folder_path, f'{formatted_date}-{n1}.mrk')
                mrk2 = os.path.join(folder_path, f'{formatted_date}-{n2}.mrk')
                elp = os.path.join(folder_path, f'{subject_id}_{date_str}_points.txt')
                hsp = os.path.join(folder_path, f'{subject_id}_{date_str}_HS.txt')

                kit_path = os.path.join(folder_path, sqd_file)
                # raw = mne.io.read_raw_kit(kit_path, elp=elp, hsp=hsp, mrk=[mrk1, mrk2], stim='164:165', preload=True)
                # print(len(raw.info['chs']))
                # print(raw.info['ch_names'])
                # events = mne.find_events(raw, stim_channel=['MISC 005','MISC 006'], min_duration=0.002)

                # Construct the mne kit2fiff command
                command = [
                    'mne', 'kit2fiff',
                    '--input', kit_path,
                    '--mrk', os.path.join(folder_path, mrk1),
                    '--mrk', os.path.join(folder_path, mrk2),
                    '--elp', os.path.join(folder_path, elp),
                    '--hsp', os.path.join(folder_path, hsp),
                    '--output', os.path.join(base_path, 'fif',f'{base_name}-raw.fif'), 
                    '--stim', '164:165', #stim channel naming convention, might break
                    '--slope','-',
                    '--stimthresh', '1',
                ]
                
                # Print command for debugging purposes
                print("Running command:", ' '.join(command))
                
                # Run the command
                try:
                    result = subprocess.run(command, check=True, capture_output=True, text=True)
                    print("Command output:", result.stdout)
                    raw = mne.io.read_raw_fif(os.path.join(bids_root,f'{base_name}-raw.fif'), preload=True)
                    events = mne.find_events(raw, min_duration=0.002)
                    print(events)
                    print(len(events))
                    print(len(raw.info['chs']))
                    # print(len(raw.info['ch_names']))
                except subprocess.CalledProcessError as e:
                    print("Command failed with error:")
                    print("Return code:", e.returncode)
                    print("Error output:", e.stderr)
                # filtered_events = events[np.isin(events[:, 2], list(event_id.values()))]  # Adjust `event_ids`
                # raw.plot(events=events, event_id=event_id)

                # bids_path = BIDSPath(subject=subject_id, task='relatedness', session=[i if i not in [subject_id, 'PP', date_str] else 'meg' for i in base_name.split('_')][0], root=bids_root)
                # input('press enter to continue...')
                # write_raw_bids(raw=raw, bids_path=bids_path, events=filtered_events, event_id=event_id, overwrite=True, format="FIF")

def process_all_folders(base_path):
    #Process all folders in the base path
    for folder_name in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder_name)
        if os.path.isdir(folder_path): #and folder_name in current_left:
            process_folder(folder_path)

# Set the base path to the directory containing all subject folders
base_path = '/raw' #experiment name here
bids_root = os.path.join(base_path, 'fif')
process_all_folders(base_path)