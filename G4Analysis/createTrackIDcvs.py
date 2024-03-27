#!/usr/bin/env python

import ROOT
import numpy as np
import pandas as pd

def particleLength( start, end ):
    return np.sqrt( np.sum([ dL**2 for dL in end-start ]) )

if __name__ == "__main__":

    # Settings for cosmic background
    # indir = '/Users/yuntse/data/coherent/preLArTPC/geant4/OneMEventsCR'
    # outdir = '/Users/yuntse/data/coherent/preLArTPC/analysis/TrackID'

    # nFiles = 1000
    # nInFilesInOne = 20

    # Settings for nueArCC signal
    indir = '/Users/yuntse/data/coherent/preLArTPC/geant4/nueArCCSignal'
    outdir = '/Users/yuntse/data/coherent/preLArTPC/analysis/nueArCCTrackID'

    nFiles = 10
    nInFilesInOne = 1


    # unit in mm, half the dimension
    larX = 300.
    larY = 250.
    larZ = 300.

    columns = [ 'Run', 'Event', 'TrackID', 'MotherID', 'Pdg', 'StartE', 'dE', 'StartX', 'StartY', 'StartZ', 
                'EndX', 'EndY', 'EndZ', 'RealTrackLength', 'StraightTrackLength']
    df = pd.DataFrame( columns = columns )

    for iFile in range(nFiles):

        if iFile%nInFilesInOne == 0:
            df = pd.DataFrame( columns = columns )

        infilename = f'{indir}/nueArCC_sns_yDir_g4_{iFile:02d}.root'
        print( f'Processing {infilename}....')
        infile = ROOT.TFile( infilename, 'READ')
        t = infile.Get("edep")

        Evt = -1
        TrackID = -1
        StartPos = np.array([0., 0., 7000.])

        for s in t:
    
            if s.pdg == 22 or np.abs(s.pdg) == 12 or np.abs(s.pdg) == 14 or np.abs(s.pdg) == 16 or np.abs(s.pdg) == 2112:
                continue
    
            mask = ( np.abs(s.startX) <= larX ) & ( np.abs(s.startY) <= larY ) & ( np.abs(s.startZ) <= larZ ) & \
                    ( np.abs(s.endX) <= larX ) & ( np.abs(s.endY) <= larY ) & ( np.abs(s.endZ) <= larZ )
            
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
            outcsv = f'{outdir}/nueArCCTrackID{iCSV:02d}.csv'
            df.to_csv( outcsv, index = False)