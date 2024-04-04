#!/usr/bin/env python

import ROOT
import numpy as np
import pandas as pd

def stepTotalLength(steps):
    return sum(np.hypot(steps.endX - steps.startX, steps.endY - steps.startY, steps.endZ - steps.startZ))

def particleLengths(particles): return particles.apply(stepTotalLength)

if __name__ == "__main__":

    indir = '/Users/yuntse/data/coherent/preLArTPC/geant4/CR2500Kv2'
    inPrefix = 'cosmic'
    outcsv = '/Users/yuntse/data/coherent/preLArTPC/analysis/EFracv2/cosmicEFrac.csv'

    nFiles = 2500

    # In the output csv file, change InnerE -> FiducialE, FiducialE -> PartConE
    columns = [ 'Run', 'Event', 'InnerE', 'FiducialE', 'OutE', 'TotalE', 'MaxE', 'TopCRTE', 'BottomCRTE', 
               'FrontCRTE', 'BackCRTE', 'LeftCRTE', 'RightCRTE' ]
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
    # Inner volume
    InX = 200.
    InY = 150.
    InZ = 200.

    for iFile in range( nFiles ):

        infile = f'{indir}/{inPrefix}_g4_{iFile:04d}.root'
        df = pd.DataFrame(ROOT.RDataFrame("edep", infile).AsNumpy())

        # Find the events and steps in TPC
        isInTPC = ((df.startX.abs() <= TPCX)&(df.startY.abs() <= TPCY)&(df.startZ.abs() <= TPCZ))| \
                  ((df.endX.abs() <= TPCX)&(df.endY.abs() <= TPCY)&(df.endZ.abs() <= TPCZ))
        inTPC = df[isInTPC&( df.dE > 0. )]

        ByEvent = inTPC.groupby('event')

        # Find the events and steps in cosmic ray taggers (CRT)
        isInTopCRT = (((df.startZ<=635)&(df.endZ>=605))|((df.startZ>=605)&(df.endZ<=635)))&(df.startX.abs()<=500)& \
                     (df.startY.abs()<=500)&(df.endX.abs()<=500)&(df.endY.abs()<=500)
        inTopCRT = df[isInTopCRT&( df.dE > 0.)].groupby('event')

        isInBottomCRT = (((df.startZ>=-635)&(df.endZ<=-605))|((df.startZ<=-605)&(df.endZ>=-635)))& \
                        (df.startX.abs()<=500)&(df.startY.abs()<=500)&(df.endX.abs()<=500)&(df.endY.abs()<=500)
        inBottomCRT = df[isInBottomCRT&( df.dE > 0.)].groupby('event')

        isInFrontCRT = (((df.startY<=-505)&(df.endY>=-535))|((df.startY>=-535)&(df.endY<=-505)))& \
                       (df.startX.abs()<=500)&(df.startZ.abs()<=600)&(df.endX.abs()<=500)&(df.endZ.abs()<=600)
        inFrontCRT = df[isInFrontCRT&( df.dE > 0.)].groupby('event')

        isInBackCRT = (((df.startY>=505)&(df.endY<=535))|((df.startY<=535)&(df.endY>=505)))& \
                      (df.startX.abs()<=500)&(df.startZ.abs()<=600)&(df.endX.abs()<=500)&(df.endZ.abs()<=600)
        inBackCRT = df[isInBackCRT&( df.dE > 0.)].groupby('event')

        isInLeftCRT = (((df.startX<=-505)&(df.endX>=-535))|((df.startX>=-535)&(df.endX<=-505)))& \
                      (df.startY.abs()<=500)&(df.startZ.abs()<=600)&(df.endY.abs()<=500)&(df.endZ.abs()<=600)
        inLeftCRT = df[isInLeftCRT&( df.dE > 0.)].groupby('event')

        isInRightCRT = (((df.startX>=505)&(df.endX<=535))|((df.startX<=535)&(df.endX>=505)))& \
                       (df.startY.abs()<=500)&(df.startZ.abs()<=600)&(df.endY.abs()<=500)&(df.endZ.abs()<=600)
        inRightCRT = df[isInRightCRT&( df.dE > 0.)].groupby('event')

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
            maskIn = (cs.startX.abs() <= InX)&(cs.startY.abs() <= InY)&(cs.startZ.abs() <= InZ)& \
                     (cs.endX.abs() <= InX)&(cs.endY.abs() <= InY)&(cs.endZ.abs() <= InZ)
            maskFi = ~maskIn & (cs.startX.abs() <= FX)&(cs.startY.abs() <= FY)&(cs.startZ.abs() <= FZ)& \
                     (cs.endX.abs() <= FX)&(cs.endY.abs() <= FY)&(cs.endZ.abs() <= FZ)
            maskOut = ~maskIn & ~maskFi & (cs.startX.abs() <= TPCX)&(cs.startY.abs() <= TPCY)&(cs.startZ.abs() <= TPCZ)& \
                     (cs.endX.abs() <= TPCX)&(cs.endY.abs() <= TPCY)&(cs.endZ.abs() <= TPCZ)
            InnerE = cs[maskIn].dE.sum()
            FiducialE = cs[maskFi].dE.sum()
            OutE = cs[maskOut].dE.sum()
            TotalE = InnerE + FiducialE + OutE
            MaxE = cs.groupby('trackID').dE.sum().max()

            # Check CRT hits
            TopCRTE = 0.
            BottomCRTE = 0.
            FrontCRTE = 0.
            BackCRTE = 0.
            LeftCRTE = 0.
            RightCRTE = 0.

            if eventNo in inTopCRT.groups.keys():
                crtEvt = inTopCRT.get_group(eventNo)
                TopCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()
            
            if eventNo in inBottomCRT.groups.keys():
                crtEvt = inBottomCRT.get_group(eventNo)
                BottomCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()
            
            if eventNo in inFrontCRT.groups.keys():
                crtEvt = inFrontCRT.get_group(eventNo)
                FrontCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()

            if eventNo in inBackCRT.groups.keys():
                crtEvt = inBackCRT.get_group(eventNo)
                BackCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()

            if eventNo in inLeftCRT.groups.keys():
                crtEvt = inLeftCRT.get_group(eventNo)
                LeftCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()
            
            if eventNo in inRightCRT.groups.keys():
                crtEvt = inRightCRT.get_group(eventNo)
                RightCRTE = crtEvt[(crtEvt.pdg.abs()==11)|(crtEvt.pdg.abs()==13)|(crtEvt.pdg.abs()==211)|(crtEvt.pdg==2212)].dE.sum()

            outCS = pd.DataFrame([{ 'Run': iFile, 'Event': eventNo, 'InnerE': InnerE, 'FiducialE': FiducialE, 'OutE': OutE, 
                                   'TotalE': TotalE, 'MaxE': MaxE, 'TopCRTE': TopCRTE, 'BottomCRTE': BottomCRTE, 'FrontCRTE': FrontCRTE,
                                    'BackCRTE': BackCRTE, 'LeftCRTE': LeftCRTE, 'RightCRTE': RightCRTE }])
            odf = pd.concat([ odf, outCS ], ignore_index = True)
    
    odf.to_csv( outcsv, index = False)
