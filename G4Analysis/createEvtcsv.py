#!/usr/bin/env python

import numpy as np
import pandas as pd

if __name__ == "__main__":

    indir = '/Users/yuntse/data/coherent/preLArTPC/analysis/trackIDv2/CR2500K'
    inprefix = 'cosmic'
    outcsv = '/Users/yuntse/data/coherent/preLArTPC/analysis/evtv2/cosmicEvt.csv'

    # Concatenate the input data frames
    nInFiles = 125
    csvFiles = [ f'{indir}/{inprefix}_TrackID_{i:04d}.csv' for i in range(nInFiles) ]
    rawdfs = []

    for csvFile in csvFiles:
        rawdf = pd.read_csv( csvFile )
        rawdfs.append( rawdf ) 
    
    indf = pd.concat( rawdfs, ignore_index = True )

    savedParticleList = np.array([ 11, -11, 13, -13, 2212, 2112, 22, 211, -211 ])
    for p in indf.Pdg.unique():
        if not np.any(savedParticleList == p):
            print( f'Particle with PDG = {p} will not be saved!')

    # Define the output data frame
    columns = [ 'Run', 'Event', 'eTotalE', 'eMaxE', 'eMaxLength', 'muTotalE', 'muMaxE', 'muMaxLength',
                'pTotalE', 'pMaxE', 'pMaxLength', 'nTotalE', 'nMaxE', 'nMaxLength',
                'gTotalE', 'gMaxE', 'gMaxLength', 'piTotalE', 'piMaxE', 'piMaxLength' ]
    df = pd.DataFrame( columns = columns )

    # Find the unique event list
    evtList = []
    for row in indf.itertuples( index = False):
        evtID = (row.Run, row.Event)
        if evtID not in evtList:
            evtList.append( evtID )

    # Loop over the events
    for iRun, iEvt in evtList:
        event = indf[(indf.Run==iRun)&(indf.Event==iEvt)]

        # Electrons
        eTotalE = 0.
        eMaxE = 0.
        eMaxLength = 0.

        if any(np.abs(event.Pdg) == 11):
            electrons = event[np.abs(event.Pdg) == 11]
            eTotalE = electrons.dE.sum()
            eMaxE = electrons.dE.max()
            eMaxLength = electrons.StraightTrackLength.max()

        # Muons
        muTotalE = 0.
        muMaxE = 0.
        muMaxLength = 0.

        if any(np.abs(event.Pdg) == 13):
            muons = event[np.abs(event.Pdg) == 13]
            muTotalE = muons.dE.sum()
            muMaxE = muons.dE.max()
            muMaxLength = muons.StraightTrackLength.max()

        # Protons
        pTotalE = 0.
        pMaxE = 0.
        pMaxLength = 0.

        if any(event.Pdg == 2212):
            protons = event[event.Pdg == 2212]
            pTotalE = protons.dE.sum()
            pMaxE = protons.dE.max()
            pMaxLength = protons.StraightTrackLength.max()

        # Neutrons
        nTotalE = 0.
        nMaxE = 0.
        nMaxLength = 0.

        if any(event.Pdg == 2112):
            neutrons = event[event.Pdg == 2112]
            nTotalE = neutrons.dE.sum()
            nMaxE = neutrons.dE.max()
            nMaxLength = neutrons.StraightTrackLength.max()

        # gammas
        gTotalE = 0.
        gMaxE = 0.
        gMaxLength = 0.

        if any(event.Pdg == 22):
            gammas = event[event.Pdg == 22]
            gTotalE = gammas.dE.sum()
            gMaxE = gammas.dE.max()
            gMaxLength = gammas.StraightTrackLength.max()

        # Charged pions
        piTotalE = 0.
        piMaxE = 0.
        piMaxLength = 0.

        if any(np.abs(event.Pdg) == 211):
            pions = event[np.abs(event.Pdg) == 211]
            piTotalE = pions.dE.sum()
            piMaxE = pions.dE.max()
            piMaxLength = pions.StraightTrackLength.max()

        # Save the event information into the data frame
        evtDF = pd.DataFrame([{ 'Run': int(iRun), 'Event': int(iEvt), 'eTotalE': eTotalE, 'eMaxE': eMaxE, 
                                'eMaxLength': eMaxLength, 'muTotalE': muTotalE, 'muMaxE': muMaxE, 
                                'muMaxLength': muMaxLength, 'pTotalE': pTotalE, 'pMaxE': pMaxE, 
                                'pMaxLength': pMaxLength, 'nTotalE': nTotalE, 'nMaxE': nMaxE, 
                                'nMaxLength': nMaxLength, 'gTotalE': gTotalE, 'gMaxE': gMaxE, 
                                'gMaxLength': gMaxLength, 'piTotalE': piTotalE, 'piMaxE': piMaxE,
                                'piMaxLength': piMaxLength }])
        
        df = pd.concat([ df, evtDF ], ignore_index = True)
    
    # Save to csv
    df.to_csv( outcsv, index = False)
