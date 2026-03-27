########################################################################
# VSMCalc.py
# Scripts to analyse DMS Model 10 VHD data files
# Leon Abelmann
########################################################################
import csv
import numpy as np
import math
from sistr import sistr

from VSMLoad import CutCurves

class Parameters():
    "Parameters extracted from hysteresis loop"
    class par():
        "Value plus error"
        val = 0
        err = 0
        unit = 1
        text = ""
        def all(self):
            return "%g (%g) %s"%(self.val,self.err,self.text)
        
    FileName   = ""    # VHD filename
    Background = par() # Linear background signal Am2/T
    Moment     = par() # Saturation moment, Am2
    Remanence  = par() # Remanent moment, Am2
    Coercivity = par() # Coercive field, T
    Slope      = par() # dm/dB at B=Bc
    Angle      = par() # Field angle (set by VSM)

    # Machine settings for Angle
    Angle.err  = 1
    Angle.text = "deg"
    
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
    def textblock(self):
        obj = self.Moment
        txt =  "Moment     : %s \n"%\
               sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        obj = self.Background
        txt = txt + "Background : %s \n"%\
              sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        obj = self.Remanence
        txt = txt + "Remanence  : %s \n"%\
              sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        obj = self.Coercivity
        txt = txt + "Coercivity : %s \n"%\
              sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        obj = self.Slope
        txt = txt + "Slope      : %s \n"%\
              sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        obj = self.Angle
        txt = txt + "Angle      : %s"%\
              sistr(obj.val/obj.unit, obj.err/obj.unit, obj.text)
        return txt


#######################################################################
# (Moment, Error, Count) = CalcMoment(Field,MomentPar)
# Estimate saturation moment over "count" points from end of curve
######################################################################
def CalcMoment(Field,MomentPar,Settings):
    FieldUp,MomentUp,FieldDown,MomentDown = CutCurves(Field,MomentPar)
    # Average moment from lowest field. Find index where deviation
    # with average becomes larger than 0.5%
    average = MomentUp[0]
    newM = MomentUp[1]
    count = 1
    while ((abs((average - newM)/average) < Settings.MomentDrop) and
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
    print("Slope Down: ", Down, " error: ", DError)

    return Up, UError, Down, DError, indexUp, indexDown
 
