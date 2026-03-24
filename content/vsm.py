#!/usr/bin/env python3
########################################################################
# vsm.py
# Script to analyse DMS Model 10 VHD data files
# Leon Abelmann
########################################################################
import sys # Enable reading parameters from command line
from pathlib import Path # To handle filenames from arg list

import VSMPlot

class VSMSettings:
    "Settings for analysis and plotting"
    
    # Moment is fitted on a range determined by a maximum drop
    MomentDrop = 0.018 # Say 1%

    # Background signal of holder plus empty substrate
    # E.g. -4.4e-8 emu/Oe. 1 emu/Oe = 1e-3 Am2/1e-4 T = 10 Am2/T
    BackgroundHolder      = -2.45e-7 # Am2/T

    BackgroundHolderError =  0.2e-7
    # If you estimate background from curve, define number of points
    BGNumPoints = 6 # Number of points to use for diamagnetic background
    # estimate

    # Number of points to calculate susceptibility
    SlopeNumPoints = 50

    # Number of points to calculate remanance
    RemNumPoints = 50

    # Sometimes first points are bad, skip them
    SkipPoints = 0

    # Extra's to display in graph
    ShowSlopes = True


def plot(filename,settings):
    return VSMPlot.plotVSM(filename, settings)

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        print("Usage: vsm.py <datafile.VHD>")
        sys.exit(2)

    for arg in sys.argv[1:]:
        # Get filename with path and extension
        filename=Path(str(arg))
        print("Processing ", filename)
        settings = VSMSettings()
        plot(filename, settings)
 
