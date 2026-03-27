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

    # Restrict the field axis (e.g. to show coercivity more clearly
    # Leave commented for autoscale on maximum range
    # FieldLimit = 0.5 # T. 
    
    # Background signal of holder plus empty substrate
    # E.g. -4.4e-8 emu/Oe. 1 emu/Oe = 1e-3 Am2/1e-4 T = 10 Am2/T
    # BackgroundHolder      = -1.67e-7 # Am2/T
    # BackgroundHolderError =  0.2e-7
    
    # If you estimate background from curve, define number of points
    BGNumPoints = 6 # Number of points to use for background estimate

    # Moment is fitted on a range determined by a maximum drop
    MomentDrop = 0.01 # Say 1%
    ShowSaturation = True # Show points used for saturation magnetisation
    
    # Number of points to calculate susceptibility
    SlopeNumPoints = 5
    ShowSlopes = True     # Show points used for susceptibility estimate

    # Number of points to calculate remanence
    RemNumPoints = 5
    ShowRemanence = False # Show remanence estimate

    # Show lines at coercivity
    ShowCoercivity = False
    
    # Sometimes first points are bad, skip them
    SkipPoints = 0

    # Generate outputfiles
    PrintSVG = False # SVG file, easy for editing in InkScape e.g
    PrintPDF = True  # PDF file, for sharing or presentations
    SaveCSV  = True  # Save curve as csv file

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
 
