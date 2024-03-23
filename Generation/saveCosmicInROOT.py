#!/usr/bin/env python

import numpy as np
import ROOT
from array import array

if __name__ == "__main__":

    infile = '/Users/yuntse/data/coherent/preLArTPC/cosmic/Cosmic4Vector.npy'
    outdir = '/Users/yuntse/data/coherent/preLArTPC/cosmic'
    nEventPerFile = 1000
    nEvent = 1000000
    nCRperEvent = 3

    # muon's constant
    muPdg = 13
    # unit: GeV
    muMass = 0.105658

    # Load the input file
    vcosmic = np.load( infile )

    iFile = 0

    for i in range(nEvent):
        if i%nEventPerFile == 0:
            iFile = int(i/nEventPerFile)
            outFile = f'{outdir}/CosmicFlux{iFile:04d}.root'

            f = ROOT.TFile(outFile, "RECREATE")
            t = ROOT.TTree("kin", "Cosmic muons")

            # Declare the variables
            event = array( 'i', [0] )
            vpdg = ROOT.std.vector('int')()
            vE   = ROOT.std.vector('double')()
            vpx  = ROOT.std.vector('double')()
            vpy  = ROOT.std.vector('double')()
            vpz  = ROOT.std.vector('double')()
            vm   = ROOT.std.vector('double')()

            # Set the branches
            t.Branch('event', event, 'event/I')
            t.Branch('pdg', vpdg)
            t.Branch('E', vE)
            t.Branch('px', vpx)
            t.Branch('py', vpy)
            t.Branch('pz', vpz)
            t.Branch('m', vm)
    
        vpdg.clear()
        vE.clear()
        vpx.clear()
        vpy.clear()
        vpz.clear()
        vm.clear()
    
        event[0] = i
    
        for iCosmic in range(nCRperEvent):
        
            p = vcosmic[i*nCRperEvent+iCosmic][0]
            costh = vcosmic[i*nCRperEvent+iCosmic][1]
            phi = vcosmic[i*nCRperEvent+iCosmic][2]
            sinth = np.sqrt(1 - costh**2)
        
            E = np.sqrt( p**2 + muMass**2 )
            px = p*sinth*np.cos(phi)
            py = p*sinth*np.sin(phi)
            pz = -p*costh
        
            vpdg.push_back( muPdg )
            vE.push_back( E )
            vpx.push_back( px )
            vpy.push_back( py )
            vpz.push_back( pz )
            vm.push_back( muMass )
    
        t.Fill()

        if (i%nEventPerFile == nEventPerFile-1):
            print( f'saving event {i} in file {iFile}....')
            f.Write()
            f.Close()        