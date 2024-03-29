#!/usr/bin/env python

import ROOT
import numpy as np
import pandas as pd

def stepTotalLength(steps):
    return sum(np.hypot(steps.endX - steps.startX, steps.endY - steps.startY, steps.endZ - steps.startZ))

def particleLengths(particles): return particles.apply(stepTotalLength)

if __name__ == "__main__":

    indir = '/Users/yuntse/data/coherent/preLArTPC/geant4/nueArCCoutFiducial'
    inPrefix = 'nueArCC_sns_yDir'
    outcsv = '/Users/yuntse/data/coherent/preLArTPC/analysis/EFrac/nueArCCoutFiducialEFrac.csv'

    nFiles = 40

    columns = [ 'Run', 'Event', 'FiducialE', 'OutE' ]
    odf = pd.DataFrame( columns = columns )
    
    # unit in mm, half the dimension
    # Fiducial volume
    FX = 250.
    FY = 200.
    FZ = 250.
    # TPC volume
    TPCX = 300.
    TPCY = 250.
    TPCZ = 300.

    for iFile in range( nFiles ):

        infile = f'{indir}/{inPrefix}_g4_{iFile:02d}.root'
        df = pd.DataFrame(ROOT.RDataFrame("edep", infile).AsNumpy())

        # Find the events and steps in TPC
        isInTPC = ((df.startX.abs() <= TPCX)&(df.startY.abs() <= TPCY)&(df.startZ.abs() <= TPCZ))| \
                  ((df.endX.abs() <= TPCX)&(df.endY.abs() <= TPCY)&(df.endZ.abs() <= TPCZ))
        inTPC = df[isInTPC&( df.dE > 0. )]

        ByEvent = inTPC.groupby('event')

        # Find events with an identifiable muon
        muonicEvents = []
        for eventNo, evt in inTPC.groupby('event'):
            muonSteps = evt[evt.pdg.abs() == 13]
            if len(muonSteps) == 0: continue # there are no muons in TPC
            muons = muonSteps.groupby('trackID')
            if (particleLengths(muons) < 50.0).all(): continue # muons in TPC are all shorter than 5 cm
            # ... more conditions?
            muonicEvents.append(eventNo) # ok, it's really a nasty muon
        print(f"In Run {iFile}, found {len(muonicEvents)}/{len(inTPC.groupby('event'))} muonic events")

        # Sum over the energy deposition from electrons, muons, charged pions, and protons in different
        # detector portions
        for eventNo, evt in inTPC.groupby('event'):
            if eventNo in muonicEvents: continue
            cs = evt[(evt.pdg.abs()==11)|(evt.pdg.abs()==13)|(evt.pdg.abs()==211)|(evt.pdg==2212)]
            maskFi = (cs.startX.abs() <= FX)&(cs.startY.abs() <= FY)&(cs.startZ.abs() <= FZ)& \
                     (cs.endX.abs() <= FX)&(cs.endY.abs() <= FY)&(cs.endZ.abs() <= FZ)
            maskOut = ~maskFi & (cs.startX.abs() <= TPCX)&(cs.startY.abs() <= TPCY)&(cs.startZ.abs() <= TPCZ)& \
                     (cs.endX.abs() <= TPCX)&(cs.endY.abs() <= TPCY)&(cs.endZ.abs() <= TPCZ)
            FiducialE = cs[maskFi].dE.sum()
            OutE = cs[maskOut].dE.sum()
            outCS = pd.DataFrame([{ 'Run': iFile, 'Event': eventNo, 'FiducialE': FiducialE, 'OutE': OutE }])
            odf = pd.concat([ odf, outCS ], ignore_index = True)
    
    odf.to_csv( outcsv, index = False)
