#######################################################################
# sistr(value, error, str, sig)
# Returns string in correct format. If
# value = 1.235, error = 0.32, unit = "m",
# sig=1 results in "1.2+-0.3 m"
# sig=2 results in "1.24+-0.32 m"
#######################################################################
def sistr(value, error, unit="", sig=1):
    # Determine significant numbers
    scale = 0
    err = error
    if (error == 0):
        err = 0
        scale = 0
    else:
        if (error < 1):
            scale = -1
            while (err < 0.1):
                err = err*10
                scale = scale-1
        while (err > 10):
            err = err/10
            scale = scale+1
        
    # Format error to sig significant numbers
    n = 10**scale
    errRound = round(error/n)*n
    valRound = round(value/n)*n
    #print(value, error, " : ", valRound, errRound, " scale: ", n)
    returnstr = r"%g$\pm$%g %s"%(valRound,errRound,unit)
    #print(value, error, " : ", returnstr, " scale: ", n)
    return returnstr

sistr(1.234,0.032,'m')
sistr(1.234,0.32,'m')
sistr(1.234,3.2,'m')
sistr(1.234,32,'m')
