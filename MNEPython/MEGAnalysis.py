import os
import mne
import numpy as np
import matplotlib.pyplot as plt

# Load the data
rootPath = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/MEG'
groupName = 'GroupB'
recRoot = os.path.join(rootPath, groupName, 'Recording')

# Load the raw data
lstFiles = os.listdir(recRoot)
sqdFiles = [f for f in lstFiles if f.endswith('_NR.sqd')]
mrkFiles = [f for f in lstFiles if f.endswith('.mrk')]
elpFiles = [f for f in lstFiles if f.endswith('Points.txt')]
hsFiles = [f for f in lstFiles if f.endswith('HS.txt')]

# print(lstFiles)
# print(mrkFiles)
# print(elpFiles)
# print(hsFiles)

for ff in sqdFiles:
    mrk = []
    for mrkIdx in range(len(mrkFiles)):
        mrk.append(os.path.join(recRoot, mrkFiles[mrkIdx]))
    elp = os.path.join(recRoot, elpFiles[0])
    hsp = os.path.join(recRoot, hsFiles[0])
    kitpath = os.path.join(recRoot, ff)
    raw = mne.io.read_raw_kit(kitpath, mrk=mrk, elp=elp, hsp=hsp,
                              verbose=True)
    print(raw.ch_names)
    stim_channels = ['MISC 001', 'MISC 002', 'MISC 003', 'MISC 004', 'MISC 005', 'MISC 006', 'MISC 007', 'MISC 008']
    raw.set_channel_types({ch: 'stim' for ch in stim_channels})
    print(raw.get_channel_types())

    for stCh in stim_channels:
        thEvents = mne.find_events(raw, stim_channel=stCh, min_duration=0.002)
        print(f"Events for {stCh}: {len(thEvents)}")
    # Compute spectrogram per channel and plot
    # raw.compute_psd(fmax=70).plot()
    # plt.show()

    # # Screen data for quality check
    # raw.plot(scalings='auto')
    # plt.show()
