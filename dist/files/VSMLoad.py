########################################################################
# VSMLoad.py
# Script load and manipulate DMS Model 10 VHD data files
# Leon Abelmann
########################################################################

#######################################################################
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

#######################################################################
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

#######################################################################
# Moment = CorrectBackground(Field,Moment, Background)
# Correct moment for background signal. Background in Am2/T
######################################################################
def CorrectBackground(Field,Moment,Background):
    MomentPar = Momen - Background*Field
    return Momen
