########################################################################
# vsm.py
# Script to analyse DMS Model 10 VHD data files
# Leon Abelmann
########################################################################
import VSMPlot

class VSMSettings():
    "Settings for analysis and plotting"
    
    # Moment is fitted on a range determined by a maximum drop
    MomentDrop = 0.018 # Say 1%

    # Background signal of holder plus empty substrate
    # E.g. -4.4e-8 emu/Oe. 1 emu/Oe = 1e-3 Am2/1e-4 T = 10 Am2/T
    BackgroundHolder      = -2.45e-7 # Am2/T

    BackgroundHolderError =  0.2e-7
    # If you estimate background from curve, define number of points
    # BGNumPoints = 6 # Number of points to use for diamagnetic background
    # estimate

    # Number of points to calculate susceptibility
    SlopeNumPoints = 50

    # Number of points to calculate remanance
    RemNumPoints = 50

    # Sometimes first points are bad, skip them
    SkipPoints = 0

def plot(filename,settings):
    return VSMPlot.PlotVSM(filename, settings)


 
