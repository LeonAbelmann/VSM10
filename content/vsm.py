#!/usr/bin/env python3
########################################################################
# vsm.py
# Script to analyse DMS Model 10 VHD data files
# Leon Abelmann
########################################################################
import sys # Enable reading parameters from command line
from pathlib import Path # To handle filenames from arg list
import os
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from sistr import sistr
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Moment is fitted on a range determined by a maximum drop
MomentDrop = 0.018 # Say 1%

# Background signal of holder plus empty substrate
# E.g. -4.4e-8 emu/Oe. 1 emu/Oe = 1e-3 Am2/1e-4 T = 10 Am2/T
BackgroundHolder      = -2.45e-7 # Am2/T

BackgroundHolderError =  0.2e-7
# If you estimate background from curve, define number of points
# BGNumPoints = 6 # Number of points to use for diamagnetic background estimate

# Number of points to calculate susceptibility
SlopeNumPoints = 50

# Number of points to calculate remanance
RemNumPoints = 50

# Sometimes first points are bad, skip them
SkipPoints = 0

# ######################################################################
# (Field, Angle, MomentPar, MomentPer) = LoadVHD(filename)
# Load file into Field [T], Angle[deg] and Moment [Am2] arrays
######################################################################
def LoadVHD(filename):
    Field     = np.array([]) # Field values
    Angle     = np.array([]) # Angle values
    MomentPar = np.array([]) # Moment parallel to field
    MomentPer = np.array([]) # Moment parallel to field

    fn = open(filename)
    
    # Data starts at @@Data
    count = 0
    line =""
    while not line.startswith('@@Data'):
        count += 1
        line = fn.readline()
    fn.readline() # Drop first line (Section...)
    count += 1
    
    # Load lines of numbers into array data
    # Data ends with @@END Data
    count= 0
    data = []
    line = fn.readline()
    while not line.startswith('@@END Data'):
        count +=1
        a_list = line.split()
        map_object = map(float, a_list)
        data.append(list(map_object))
        line = fn.readline()

    for row in data:
        FieldValue  = row[8]*1e-4 # Convert Oe to T 
        FieldAngle = row[6]
        MParValue  = row[11]
        MPerValue  = row[12]
        #print(FieldValue,MParValue,MPerValue,FieldAngle)
        Field     = np.append(Field,FieldValue)
        Angle     = np.append(Angle,FieldAngle)
        MomentPar = np.append(MomentPar,MParValue)
        MomentPer = np.append(MomentPer,MPerValue)

    return (Field, Angle, MomentPar, MomentPer)

# ######################################################################
# MomentPar = CorrectBackground(Field,MomentPar, Background)
# Correct moment for background signal. Background in Am2/T
######################################################################
def CorrectBackground(Field,MomentPar,Background):
    MomentPar = MomentPar - Background*Field
    return MomentPar


# ######################################################################
# FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,Moment)
# Return op and down curves seperately
######################################################################
def CutCurves(Field,Moment):
    # Find index of largest value in Field
    halfway = np.argmax(Field)
    #print(halfway)
    # Cut curves
    FieldUp = Field[1:halfway]
    MomentUp = Moment[1:halfway]
    FieldDown = Field[halfway:]
    MomentDown = Moment[halfway:]
    return FieldUp, MomentUp, FieldDown, MomentDown

# ######################################################################
# (Moment, Error, Count) = CalcMoment(Field,MomentPar)
# Estimate saturation moment over "count" points from end of curve
######################################################################
def CalcMoment(Field,MomentPar):
    FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,MomentPar)
    # Average moment from lowest field. Find index where deviation
    # with average becomes larger than 0.5%
    average = MomentUp[0]
    newM = MomentUp[1]
    count = 1
    while ((abs((average - newM)/average) < MomentDrop) and
               (count < len(Field)-1)) :
        count = count + 1
        average = (average*(count-1) + newM)/count
        #print(count, average,newM,abs((average - newM)/average))
        newM = MomentUp[count]
    
    MomentUpStart   = np.average(MomentUp[0:count])
    MomentUpEnd     = np.average(MomentUp[len(MomentUp)-count:])
    MomentDownStart = np.average(MomentDown[0:count])
    MomentDownEnd   = np.average(MomentDown[len(MomentDown)-count:])
    
    # print("Count ", count)
    # print("MomentUpStart   ",MomentUpStart)
    # print("MomentUpEnd     ",MomentUpEnd)
    # print("MomentDownStart ",MomentDownStart)
    # print("MomentDownEnd   ",MomentDownEnd)
    moments = [abs(MomentUpStart), abs(MomentUpEnd),
                   abs(MomentDownStart), abs(MomentDownEnd)]
    Moment = np.average(moments)
    Error  = (np.max(moments)-np.min(moments))/2
    print("Moment : ", Moment, " Error: ", Error)
    return Moment, Error, count

# ######################################################################
# Background, Error = CalcBackground(Field, MomentPar, NumPoints)
# Estimate slope of curve over first or last NumPoints of the curves
######################################################################
def CalcBackground(Field, MomentPar, NumPoints):
    halfway = np.argmax(Field)
    length = len(Field)
    sections = [[0,NumPoints],
                [halfway-NumPoints,halfway],
                [halfway, halfway+NumPoints],
                [length-NumPoints,length]]
    i = 0;
    slopes = []
    for sec in sections:
        slope = np.polyfit(Field[sec[0]:sec[1]],
                            MomentPar[sec[0]:sec[1]], 1)[0]
        #print("Slope: ",slope)
        slopes = np.append(slopes,slope)
        i = i + 1
    #print("Slopes: ", slopes)
    Background = np.average(slopes)
    Error  = (np.max(slopes)-np.min(slopes))/2
    print("Background : ", Background, " Error: ", Error)
    return Background, Error

#######################################################################
# Bc = findZeroCrossing(Field, Moment)
# Find Field at which moment crosses zero by least square fit of
# straight line
######################################################################
def findZeroCrossing(Field,Moment):
# https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([Field,np.ones(len(Field))]).T
    # Fit y = mx + c
    m,c = np.linalg.lstsq(A,Moment,rcond=None)[0]
    # Zero crossing for y=0, so x=-c/m
    Bc = -c/m
    return Bc

#######################################################################
# Bc = findRemanence(Field, Moment)
# Find moment for B=0
######################################################################
def findRemanence(Field,Moment):
# https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
    A = np.vstack([Field,np.ones(len(Field))]).T
    # Fit y = mx + c
    m,c = np.linalg.lstsq(A,Moment,rcond=None)[0]
    # Remancence at x=0, so return c
    return c

# Basis functions for linear model:  func = a0*X0 + a1*X1 + a2*X2 + ...
def basis(x):
    ''' Returns basis functions for linear model
    
    The function to be fit is assumed to be of the form
    f(x) = a0*X0 + a1*X1 + a2*X2 + a3*X3 + ...
    where a0, a1, a2, ... are constants, and X0, X1, X2, ... are defined below.
    '''
    #X2 = x**2   
    X1 = x
    X0 = 0.*X1 + 1. # Need array of len(x), thus the 0.*X1
    return np.array([X0,X1])

#######################################################################
# Slope, Error = findSlope(Field, Moment)
# Fit to straight line, return slope
######################################################################
def findSlope(Field,Moment):
#http://www.eg.bucknell.edu/~phys310/jupyter/linear_fit_example.html
    X = basis(Field).T    # Basis functions evaluated at all x (the X_j(x_i)) of N.R.)
    u = np.ones(len(Moment))*3e-9 # Uncertainty in Moment
    W = np.diag(1/u)  # Matrix with uncertainties on diagonal
    Xw = np.dot(W,X)  # A_ij of Eq. (14.3.4)
    Yw = np.dot(Moment,W)  # b_i of Eq. (14.3.5)
    fit = np.linalg.lstsq(Xw,Yw,rcond=1)  # lstq returns: best values, chi2, ..
    covariance = np.linalg.inv(np.dot(Xw.T,Xw))
    uncertainty = np.sqrt(np.diag(covariance))
    slope = fit[0][1]
    error = np.sqrt(covariance[1,1])
    #print("Slope",slope)
    #print("Error", error)
    return slope, error

#######################################################################
# BcUp, BcUpError, BcDown, BcDownError  = CalcBc(Field, MomentPar)
# Interpolate the zero crossing to calculate coercivity
######################################################################
def CalcBc(Field, MomentPar):
    # Split into up and down part
    FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,MomentPar)

    # Find first index where moment crosses zero
    crossing = np.where(np.diff(np.sign(MomentUp)))[0][0]
    #print(crossing)

    # Find zero crossing with two sets of points
    f1 = findZeroCrossing(FieldUp[crossing-2:crossing+3],
                         MomentUp[crossing-2:crossing+3])
    f2 = findZeroCrossing(FieldUp[crossing-3:crossing+2],
                         MomentUp[crossing-3:crossing+2])
    # print(crossing,f1,f2)
    # for i in range(6):
    #     print(crossing-2+i,\
    #               FieldUp[crossing-2+i],MomentUp[crossing-2+i])
    BcUp     = (f1+f2)/2 # Bc is average of both interpolations
    BcUpError= abs(f1-f2) # Error is difference
    
    # Same for down curve
    crossing = np.where(np.diff(np.sign(MomentDown)))[0][0]
    #print(crossing)

    # Find zero crossing with two sets of points
    f1 = findZeroCrossing(FieldDown[crossing-2:crossing+3],
                         MomentDown[crossing-2:crossing+3])
    f2 = findZeroCrossing(FieldDown[crossing-3:crossing+2],
                         MomentDown[crossing-3:crossing+2])

    # print(crossing,f1(0),f2(0))
    # for i in range(6):
    #     print(crossing-2+i,\
    #               FieldDown[crossing-2+i],MomentDown[crossing-2+i])
    BcDown     = (f1+f2)/2 # Bc is average of both interpolations
    BcDownError= abs(f1-f2) # Error is difference
    
    print("BcUp  : ", BcUp, " error: ", BcUpError)
    print("BcDown: ", BcDown, " error: ", BcDownError)
    return BcUp, BcUpError, BcDown, BcDownError

#######################################################################
# BrUp, BrUpError, BrDown, BrDownError  = CalcBr(Field, MomentPar)
# Interpolate around B=0 to calculate remanence
######################################################################
def CalcBr(Field, MomentPar):
    # Split into up and down part
    FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,MomentPar)
    
    # Find first index where field crosses zero
    crossing = np.where(np.diff(np.sign(FieldUp)))[0][0]
    #print("Crossing Up: ", crossing, "Field : ", FieldUp[crossing])

    # Find moment at B=0, do it with two different sets of points:
    m0 = findRemanence(FieldUp[crossing-2:crossing+3],
                         MomentUp[crossing-2:crossing+3])
    m1 = findRemanence(FieldUp[crossing-1:crossing+4],
                         MomentUp[crossing-1:crossing+4])
    BrUp     = (m0+m1)/2  # Take average
    BrUpError= abs(m0-m1) # Error is difference
    
    # Same for down curve
    m0 = findRemanence(FieldDown[crossing-2:crossing+3],
                         MomentDown[crossing-2:crossing+3])
    m1 = findRemanence(FieldDown[crossing-1:crossing+4],
                         MomentDown[crossing-1:crossing+4])
    BrDown     = (m0+m1)/2  # Take average
    BrDownError= abs(m0-m1) # Error is difference
    
    print("BrUp  : ", BrUp, " error: ", BrUpError)
    print("BrDown: ", BrDown, " error: ", BrDownError)
    return BrUp, BrUpError, BrDown, BrDownError

#######################################################################
# Up, Uerror, Down, DError, indexUp, indexDown =
# CalcSlope(Field, Momentpar, BcUp, BcDown, NumPoints)
# Find slope of curve around coercivity
######################################################################
def CalcSlope(Field, Moment, BcUp, BcDown, NumPoints):
    # Split into up and down part
    FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,Moment)
    
    indexUp    = np.where(FieldUp>BcUp)[0][0]
    indexDown  = np.where(FieldDown<BcDown)[0][0]
    print("MomentUp   : ", MomentUp[indexUp],
          "MomentDown : ", MomentDown[indexDown])

    fitrange = math.ceil(NumPoints/2)
    
    Up, UError =  findSlope(FieldUp[indexUp-fitrange:indexUp+fitrange],
                    MomentUp[indexUp-fitrange:indexUp+fitrange])
    print("Slope Up  : ", Up, " error: ", UError)

    Down, DError = findSlope(FieldDown[indexDown-fitrange:indexDown+fitrange],
                    MomentDown[indexDown-fitrange:indexDown+fitrange])
    print("Slope Down: ", -Down, " error: ", DError)

    # Return negative of down slope, more logical
    return Up, UError, -Down, DError, indexUp, indexDown



class Parameters():
    "Parameters extracted from hysteresis loop"
    class par():
        "Value plus error"
        val = 0
        err = 0
        def all(self):
            return "%g (%g)"%(self.val,self.err)
    FileName   = ""    # VHD filename
    Moment     = par() # Saturation moment, Am2
    Remanence  = par() # Remanent moment, Am2
    Coercivity = par() # Coercive field, T
    Slope      = par() # dm/dB at B=Bc
    Angle      = par() # Field angle (set by VSM)
    def print(self):
        print("FileName   : ",self.FileName)
        print("Angle      : ",self.Angle.all())
        print("Moment     : ",self.Moment.all())
        print("Remanence  : ",self.Remanence.all())
        print("Coercivity : ",self.Coercivity.all())
        print("Slope      : ",self.Slope.all())
    def header(self):
        txt = "# Filename,"
        txt = txt + "Angle (deg),AngleError (deg),"
        txt = txt + "Moment (Am2),MomentError (Am2),"
        txt = txt + "Remanence (Am2),RemanenceError (Am2),"
        txt = txt + "Coercivity (T),CoercivityError (T),"
        txt = txt + "Slope (Am2/T),SlopeError (Am2/T)"
        return txt
    def println(self):
        txt = self.FileName+", "
        txt = txt + "%g, %g, "%(self.Angle.val, self.Angle.err)
        txt = txt + "%g, %g, "%(self.Moment.val, self.Moment.err)
        txt = txt + "%g, %g, "%(self.Remanence.val, self.Remanence.err)
        txt = txt + "%g, %g, "%(self.Coercivity.val, self.Coercivity.err)
        txt = txt + "%g, %g  "%(self.Slope.val, self.Slope.err)
        return txt


#######################################################################
# plotVSM(filename)
# Extract parameters. Plot curve to pdf file and
# Parameters:
# Result.Moment(.val/.err) : Saturation moment
# 
#######################################################################
def plotVSM(filename):
    print("Analyzing ",filename)
    
    # Prepare results object
    Results = Parameters()
    Results.FileName = filename
    
    # Load data from file
    (Field,Angle,MomentPar,MomentPer)=LoadVHD(filename)
    print("data loaded")
    
    # Skip first measurement points
    Field     = Field[SkipPoints:] 
    Angle     = Angle[SkipPoints:]
    MomentPar = MomentPar[SkipPoints:]
    MomentPer = MomentPer[SkipPoints:]
 
    # Either define background from holder
    Background = BackgroundHolder
    BGError    = BackgroundHolderError
    # Or calculate diamagnetic background signal
    # Background, BGError = CalcBackground(Field, MomentPar, BGNumPoints)
    
    # Correct parallel curve for diamagnetic background signal
    MomentPar = CorrectBackground(Field,MomentPar,Background)

    # Assume the VSM loop is for one angle:
    Results.Angle.val = Angle[0]
    Results.Angle.err = 0
    
    # Calculate moment with error and section of curve that was used.
    Moment, MError, NumPoints = CalcMoment(Field, MomentPar)

    Results.Moment.val = Moment
    Results.Moment.err = MError
    
    # Calculate Remanence
    BrUp, BrUpError, BrDown, BrDownError  = CalcBr(Field, MomentPar)
    Br = (BrDown-BrUp)/2 # Remanence is avarage
    # Error is sum of fit errors or difference between values,
    # whichever is bigger
    BrError = max((BrUpError + BrDownError)/math.sqrt(2),\
                     abs(BrUp + BrDown))

    Results.Remanence.val = Br
    Results.Remanence.err = BrError
    
    # Calculate Coercivity
    BcUp, BcUpError, BcDown, BcDownError  = CalcBc(Field, MomentPar)
    Bc = (BcUp-BcDown)/2 # Average
    BcError = (BcUpError+BcDownError)/math.sqrt(2)
    Results.Coercivity.val = Bc
    Results.Coercivity.err = BcError
 
    # Calculate slope around Coercivity
    SlopeUp, SUError, SlopeDown, SDError, indexUp, indexDown = CalcSlope(Field, MomentPar, BcUp, BcDown, SlopeNumPoints)    
    Slope     = (SlopeUp + SlopeDown)/2
    # Error is sum of fit errors or difference between values,
    # whichever is bigger
    SlopeError = max((SUError + SDError)/math.sqrt(2),
                     abs(SlopeUp - SlopeDown))
    Results.Slope.val = Slope
    Results.Slope.err = SlopeError
                    
    # Choose field label
    if (abs(Field[0]) < 1) :
        # Field in mT
        FieldText = "mT"
        FieldUnit = 1e-3
    else:
        # Field in T
        FieldText = "T"
        FieldUnit = 1

    # Choose magnetic moment label
    if (abs(MomentPar[0]) > 1e-2) :
        #Moment in mAm2
        MomentText = "mAm$^2$"
        MomentUnit = 1e-3
    else:
        if  (abs(MomentPar[0]) > 1e-5) :
            # Moment in uAm2
            MomentText = u'\u03bc' + "Am$^2$" # Upright mu
            MomentUnit = 1e-6
        else:
            MomentText = "nAm$^2$"
            MomentUnit = 1e-9

    # Plot
    fig, ax = plt.subplots(figsize=(6,4), dpi = 100)
    
    plt.title(filename,fontsize=10)
    
    # Labels on axes
    ax.set_ylabel('Moment / '+MomentText, fontsize = 14)
    ax.set_xlabel('Field / '+FieldText, fontsize = 14)
    # Horizontal and vertical line at (0,0)
    plt.axhline(y=0, color='k', linewidth=1)
    plt.axvline(x=0, color='k', linewidth=1)
    # Tics on axis
    # ax.xaxis.set_major_locator(MultipleLocator(0.5))
    # ax.xaxis.set_major_formatter('{x:.1f}')
    # ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # Tics on all axes
    ax.tick_params(which='both',top=True,right=True)

    # Limit field axis
    #plt.xlim(-xlimit,xlimit)
    plt.ylim(-1.1*Moment/MomentUnit,1.1*Moment/MomentUnit)

    # Plot MomentPer, below all other graphs
    ax.plot(Field/FieldUnit,MomentPer/MomentUnit,
                'ro',markersize="3",label=r"m$_\perp$")
    ax.plot(Field/FieldUnit, MomentPer/MomentUnit,
                'r',linewidth="1")
    
    # Show coercivity lines:
    # plt.axvline(x=BcUp/FieldUnit)
    # plt.axvline(x=BcDown/FieldUnit)
    
    # Plot MomentPar
    ax.plot(Field/FieldUnit, MomentPar/MomentUnit,'ko',markersize="3",
                label=r"m$_\parallel$")
    ax.plot(Field/FieldUnit, MomentPar/MomentUnit,'k',linewidth="1")

    # Show slope lines used for susceptilibity
    # slopeRange=math.ceil(SlopeNumPoints/2)
    # SlopeField=Field[indexUp-slopeRange:indexUp+slopeRange]
    # ax.plot(SlopeField/FieldUnit,\
    #         (SlopeUp*(SlopeField-BcUp))/MomentUnit,'b')
    # halfway = np.argmax(Field)
    # SlopeField=Field[halfway+indexDown-slopeRange:\
    #                  halfway+indexDown+slopeRange]
    # ax.plot(SlopeField/FieldUnit,\
    #         (SlopeDown*(SlopeField-BcDown))/MomentUnit,'b')
    
    # Show which points were used to calculate moment
    # length  = len(Field)
    # index = [[0,NumPoints],
    #              [halfway-NumPoints,halfway+NumPoints],
    #              [length-NumPoints,length]]
    # for sec in index:
    #     ax.plot(Field[sec[0]:sec[1]]/FieldUnit,\
    #                 MomentPar[sec[0]:sec[1]]/MomentUnit,\
    #                 'bo',markersize="3")
                    
    # Show remanence points
    # print("Br : ", [BrUp/MomentUnit,BrDown/MomentUnit])
    # ax.plot([0,0],[BrUp/MomentUnit,BrDown/MomentUnit],'bo',markersize="5")

    # Info block
    # textblock = ""
    # textblock =  "Moment     : %s \n"%\
    #   sistr(Moment/MomentUnit, MError/MomentUnit, MomentText)
    # # Background in uAm2/T
    # textblock = textblock + "Background : %s \n"%\
    #   sistr(Background/(MomentUnit/FieldUnit),\
    #   BGError/(MomentUnit/FieldUnit), \
    #    MomentText+"/"+FieldText)
    # textblock = textblock + "Remanence  : %s \n"%\
    #   sistr(Br/MomentUnit,BrError/MomentUnit,MomentText)
    # textblock = textblock + "Bc         : %s \n"%\
    #   sistr(Bc/FieldUnit,BcError/FieldUnit,FieldText)
    # textblock = textblock + "Rel. Slope : %s \n"%\
    #   sistr(Slope/Moment,SlopeError/Moment,"1/"+FieldText)
    textblock = "Angle : %.1f deg"%\
                (Angle[0])

    plt.annotate(textblock, xy=(0.01,0.95), xycoords='axes fraction',\
                 fontfamily='monospace', fontsize=9)

    ax.legend(loc='lower right')
    # Text instead of legend
    # plt.text(10,240,r"m$_\parallel$", fontsize=12)
    # plt.text(10, 50,r"m$_\perp$", color = "r", fontsize=12)
        
    # ax[0].grid(which='major', axis='x')
    # ax[0].legend()
    # ax[0].set_ylim(ymin=0, ymax=160)

    # fig.tight_layout()
    # fig.show()
    display(fig)
    
    print("plotting done. preparing pdf")
    
    #pdffile = filename.with_suffix('.pdf')
    root, _ = os.path.splitext(filename)
    #pdffile = root + ".pdf"
    #fig.savefig(pdffile)

    root, _ = os.path.splitext(filename)    
    svgfile = root + ".svg"
    try:
        fig.savefig(svgfile)
    except:
        print("saving file failed")
    return Results
 


# if __name__ == "__main__":
#     if (len(sys.argv) < 2):
#         print("Usage: vsm.py <datafile.VHD>")
#         sys.exit(2)

#     directory = Path(sys.argv[1]).parent
#     with open(str(directory)+'/vsm_zoom.csv', 'w') as csvfile:
#         writer = csv.writer(csvfile)
#         Results = Parameters()
#         csvfile.write(Results.header()+"\n")
#         for arg in sys.argv[1:]:
#             # Get filename with path and extension
#             filename=Path(str(arg))
#             print("Processing ", filename)
#             # Extract parameters and plot loops:
#             Results = plotVSM(filename)
#             # Print results to stdout
#             # Results.print()
#             # Append results to output file
#             csvfile.write(Results.println()+"\n")

 
