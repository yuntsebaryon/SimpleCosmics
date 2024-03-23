#!/usr/bin/env python

import numpy as np
import pandas as pd

if __name__ == "__main__":

    incsv = '/Users/yuntse/data/coherent/preLArTPC/analysis/TrackID.csv'
    outcsv = '/Users/yuntse/data/coherent/preLArTPC/analysis/Evt.csv'

    indf = pd.read_csv( incsv )

    columns = [ 'Run', 'Event', 'eTotalE', 'eMaxE', 'eMaxLength', 'muLength', 'muTotalE' ]
    df = pd.DataFrame( columns = columns )

    evtID = np.array([ -1, -1 ])
    eTotalE = 0.
    eMaxE = 0.

    for row in indf.itertuples( index = False):
        if ((row.Run, row.Event) != evtID).any():
            evtID = np.array([ row.Run, row.Event ])
            if np.abs(row.Pdg) == 13:
                df_temp = pd.DataFrame([{ 'Run': row.Run, 'Event': row.Event, 'eTotalE': 0., 'eMaxE': 0.,
                                            'eMaxLength': 0., 'muLength': row.StraightTrackLength,
                                            'muTotalE': row.dE }])
            
                df = pd.concat([ df, df_temp ], ignore_index = True)
            else:
                df_temp = pd.DataFrame([{ 'Run': row.Run, 'Event': row.Event, 'eTotalE': row.dE, 'eMaxE': row.dE,
                                            'eMaxLength': row.StraightTrackLength, 'muLength': 0.,
                                            'muTotalE': 0. }])
            
                df = pd.concat([ df, df_temp ], ignore_index = True)
        else:
            if np.abs(row.Pdg) == 13:
                if row.StraightTrackLength > df.loc[df.index[-1], 'muLength']:
                    df.loc[df.index[-1], 'muLength'] = row.StraightTrackLength
                df.loc[df.index[-1], 'muTotalE'] += row.dE
            else:
                df.loc[df.index[-1], 'eTotalE'] += row.dE
                if row.dE > df.loc[df.index[-1], 'eMaxE']:
                    df.loc[df.index[-1], 'eMaxE'] = row.dE
                    df.loc[df.index[-1], 'eMaxLength'] = row.StraightTrackLength
    
    df.to_csv( outcsv, index = False)