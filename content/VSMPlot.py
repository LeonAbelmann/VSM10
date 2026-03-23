########################################################################
# VSMPlot.py
# Script to analyse and plot DMS Model 10 VHD data files
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

from VSMCalc import Parameters, CalcMoment, CalcBr, CalcBc, CalcSlope
from VSMLoad import LoadVHD, CorrectBackground, CutCurves

#######################################################################
# plotVSM(filename, settings)
# Extract parameters. Plot curve to pdf file and
# Parameters:
# Result.Moment(.val/.err) : Saturation moment
# 
#######################################################################
def plotVSM(filename, settings):
    print("Analyzing ",filename)
    
    # Prepare results object
    Results = Parameters()
    Results.FileName = filename
    
    # Load data from file
    (Field,Angle,MomentPar,MomentPer)=LoadVHD(filename)
    print("data loaded")
    
    # Skip first measurement points
    Field     = Field[settings.SkipPoints:] 
    Angle     = Angle[settings.SkipPoints:]
    MomentPar = MomentPar[settings.SkipPoints:]
    MomentPer = MomentPer[settings.SkipPoints:]
 
    # Either define background from holder
    Background = settings.BackgroundHolder
    BGError    = settings.BackgroundHolderError
    # Or calculate diamagnetic background signal
    # Background, BGError = CalcBackground(Field, MomentPar, BGNumPoints)
    
    # Correct parallel curve for diamagnetic background signal
    MomentPar = CorrectBackground(Field,MomentPar,Background)

    # Assume the VSM loop is for one angle:
    Results.Angle.val = Angle[0]
    Results.Angle.err = 0
    
    # Calculate moment with error and section of curve that was used.
    Moment, MError, NumPoints = CalcMoment(Field, MomentPar, settings)

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
    SlopeUp, SUError, SlopeDown, SDError, indexUp, indexDown = CalcSlope(Field, MomentPar, BcUp, BcDown, settings.SlopeNumPoints)    
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
 


 
