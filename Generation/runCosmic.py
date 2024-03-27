#!/usr/bin/env python

import os

if __name__ == "__main__":

    indir  = '/Users/yuntse/data/coherent/preLArTPC/cosmic/CR1M'
    outdir = '/Users/yuntse/data/coherent/preLArTPC/geant4/CR1M'
    nEventPerFile = 1000
    nFiles = 1000

    for iFile in range(nFiles):

        script = f'{outdir}/mac/run_cosmic_{iFile:04d}.mac'
        infile = f'{indir}/CosmicFlux{iFile:04d}.root'
        outfile = f'{outdir}/cosmic_g4_{iFile:04d}.root'

        with open( script, 'w') as f:
            f.write(f"""
/LArG4/generator/infile {infile}
/LArG4/output/outfile {outfile}
/run/beamOn {nEventPerFile}
""")
            
        cmd = f'./SimpleLArG4 {script}'
        os.system( cmd )
