import numpy as np
import os, mne
from datetime import datetime
from mne_bids import write_raw_bids, BIDSPath
from mne.channels import make_dig_montage
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the data
rootPath = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/MEG'
bids_root = '/Users/mrugankdake/Library/CloudStorage/GoogleDrive-mdd9787@nyu.edu/My Drive/Coursework/EEG MEG methods/ClassData/MEGBids'
if not os.path.exists(bids_root):
    os.makedirs(bids_root)
groupName = 'GroupA'
recRoot = os.path.join(rootPath, groupName, 'Recording')
timeStampsRoot = os.path.join(rootPath, groupName, 'TimeStamps')

# Filter files based on their extensions
lstFiles = os.listdir(recRoot)
sqdFiles = [f for f in lstFiles if f.endswith('_NR.sqd')]
mrkFiles = [f for f in lstFiles if f.endswith('.mrk')]
elpFiles = [f for f in lstFiles if f.endswith('Points.txt')]
hsFiles = [f for f in lstFiles if f.endswith('HS.txt')]

# Load the data
for sqdFile in sqdFiles:
    # Load the sqd file
    raw = mne.io.read_raw_ctf(os.path.join(recRoot, sqdFile), preload=True)
    # Load the mrk file
    mrkFile = [f for f in mrkFiles if sqdFile.split('_')[0] in f][0]
    mrkData = pd.read_csv(os.path.join(recRoot, mrkFile), sep='\t')
    # Load the elp file
    elpFile = [f for f in elpFiles if sqdFile.split('_')[0] in f][0]
    elpData = pd.read_csv(os.path.join(recRoot, elpFile), sep='\t')
    # Load the hs file
    hsFile = [f for f in hsFiles if sqdFile.split('_')[0] in f][0]
    hsData = pd.read_csv(os.path.join(recRoot, hsFile), sep='\t')
    # Load the time stamps
    timeStampsFile = [f for f in os.listdir(timeStampsRoot) if sqdFile.split('_')[0] in f][0]
    timeStamps = loadmat(os.path.join(timeStampsRoot, timeStampsFile))
    # Create the BIDS path
    bids_path = BIDSPath(subject='01', session='01', task='rest', run='01', root=bids_root)
    # Create the montage
    montage = make_dig_montage(ch_pos=dict(zip(hsData['Channel'], hsData[['X', 'Y', 'Z']].values)),
                                elp=dict(zip(elpData['Label'], elpData[['X', 'Y', 'Z']].values)),
                                coord_frame='head')
    raw.set_montage(montage)
    # Write the BIDS data
    write_raw_bids(raw, bids_path, overwrite=True)