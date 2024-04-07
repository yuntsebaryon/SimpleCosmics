#!/usr/bin/env python

import sys
import ROOT
import numpy as np
import pandas as pd

def particleLength( start, end ):
    return np.sqrt( np.sum([ dL**2 for dL in end-start ]) )

if __name__ == "__main__":

    # Default
    startFile = 0
    nFiles = 1

    startFile = int(sys.argv[1])
    nFiles = int(sys.argv[2])

    indir = '/Users/yuntse/data/coherent/preLArTPC/geant4/brn500Kv2'
    inPrefix = 'brn'
    outdir = '/Users/yuntse/data/coherent/preLArTPC/analysis/trackIDv2/brn500K'
    outPrefix = 'brn'

    # Settings for cosmic background
    # nInFilesInOne = 20

    # Settings for BRN
    nInFilesInOne = 10
    
    # Settings for nueArCC signal
    # nInFilesInOne = 1


    # unit in mm, half the dimension
    larX = 300.
    larY = 250.
    larZ = 300.

    columns = [ 'Run', 'Event', 'TrackID', 'MotherID', 'Pdg', 'StartE', 'dE', 'StartX', 'StartY', 'StartZ', 
                'EndX', 'EndY', 'EndZ', 'RealTrackLength', 'StraightTrackLength']
    df = pd.DataFrame( columns = columns )

    for iFile in range(startFile, startFile+nFiles):

        if iFile%nInFilesInOne == 0:
            df = pd.DataFrame( columns = columns )

        infilename = f'{indir}/{inPrefix}_g4_{iFile:04d}.root'
        print( f'Processing {infilename}....')
        infile = ROOT.TFile( infilename, 'READ')
        t = infile.Get("edep")

        Evt = -1
        TrackID = -1
        StartPos = np.array([0., 0., 7000.])

        for s in t:
    
            if np.abs(s.pdg) == 12 or np.abs(s.pdg) == 14 or np.abs(s.pdg) == 16 or s.pdg > 1000100000:
                continue
    
            mask = (( np.abs(s.startX) <= larX ) & ( np.abs(s.startY) <= larY ) & ( np.abs(s.startZ) <= larZ )) | \
                    (( np.abs(s.endX) <= larX ) & ( np.abs(s.endY) <= larY ) & ( np.abs(s.endZ) <= larZ ))
            
            if mask:
                if s.event == Evt and s.trackID == TrackID:
                    thisStartPos = np.array([ s.startX, s.startY, s.startZ ])
                    EndPos = np.array([ s.endX, s.endY, s.endZ ])
                    df.loc[df.index[-1], 'dE'] += s.dE
                    df.loc[df.index[-1], 'EndX'] = s.endX
                    df.loc[df.index[-1], 'EndY'] = s.endY
                    df.loc[df.index[-1], 'EndZ'] = s.endZ
                    df.loc[df.index[-1], 'RealTrackLength'] += particleLength( thisStartPos, EndPos )
                    df.loc[df.index[-1], 'StraightTrackLength'] = particleLength( StartPos, EndPos )
                
                else:
                    Evt = s.event
                    TrackID = s.trackID
                    StartPos[:] = [ s.startX, s.startY, s.startZ ]
                    EndPos = np.array([ s.endX, s.endY, s.endZ ])
                    RealTrackLength = particleLength( StartPos, EndPos )
                
                    df_temp = pd.DataFrame([ { 'Run': iFile, 'Event': s.event, 'TrackID': s.trackID, 
                                                'MotherID': s.motherID, 
                                                'Pdg': s.pdg, 'StartE': s.startE, 'dE': s.dE, 'StartX': s.startX, 
                                                'StartY': s.startY, 'StartZ': s.startZ, 'EndX': s.endX, 'EndY': s.endY,
                                                'EndZ': s.endZ, 'RealTrackLength': RealTrackLength, 
                                                'StraightTrackLength': RealTrackLength } ])
                    df = pd.concat([ df, df_temp ], ignore_index = True )

        infile.Close()

        if iFile%nInFilesInOne == (nInFilesInOne-1):
            iCSV = np.floor(iFile/nInFilesInOne).astype(int)
            outcsv = f'{outdir}/{outPrefix}_TrackID_{iCSV:04d}.csv'
            df.to_csv( outcsv, index = False)
