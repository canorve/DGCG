
import os.path
import os
import stat
from astropy.io import fits
import numpy as np
import sys
import scipy.special
import scipy
import subprocess as sp

from dgcg.lib import core 
from dgcg.lib import mgetools



def PrintVar(parvar,FileHandle):
    "Print DGCG user parameters in log file"

# parameters file

    line = "Img = {} \n".format(parvar.Img)
    FileHandle.write(line)

    line = "SexCat = {} \n".format(parvar.SexCat)
    FileHandle.write(line)

    line = "SigImg = {} \n".format(parvar.SigImg)
    FileHandle.write(line)

    line = "PsfDir = {} \n".format(parvar.PsfDir)
    FileHandle.write(line)

    line = "MagZpt = {} \n".format(parvar.MagZpt)
    FileHandle.write(line)

    line = "PlateScale = {} \n".format(parvar.PlateScale)
    FileHandle.write(line)

    line = "FitFunc = {} \n".format(parvar.FitFunc)
    FileHandle.write(line)

#    line = "GalBot = {} \n".format(parvar.GalBot)
#    FileHandle.write(line)

    line = "GalClas = {} \n".format(parvar.GalClas)
    FileHandle.write(line)

    line = "ConvBox = {} \n".format(parvar.ConvBox)
    FileHandle.write(line)

    line = "FitBox = {} \n".format(parvar.FitBox)
    FileHandle.write(line)

    line = "MagDiff = {} \n".format(parvar.MagDiff)
    FileHandle.write(line)

    line = "KronScale = {} \n".format(parvar.KronScale)
    FileHandle.write(line)

    line = "SkyScale = {} \n".format(parvar.SkyScale)
    FileHandle.write(line)

    line = "Offset = {} \n".format(parvar.Offset)
    FileHandle.write(line)

    line = "SkyAnnuli = {} \n".format(parvar.SkyWidth)
    FileHandle.write(line)

    line = "NSer = {} \n".format(parvar.NSer)
    FileHandle.write(line)

    line = "MaxFit = {} \n".format(parvar.MaxFit)
    FileHandle.write(line)

#    line = "MagMin = {} \n".format(parvar.MagMin)
#    FileHandle.write(line)

    line = "MagMax = {} \n".format(parvar.MagMax)
    FileHandle.write(line)

    line = "FlagSex = {} \n".format(parvar.FlagSex)
    FileHandle.write(line)

    line = "ConsFile = {} \n".format(parvar.ConsFile)
    FileHandle.write(line)

    line = "Region = {} \n".format(parvar.Region)
    FileHandle.write(line)

    line = "Bxmin = {} ".format(parvar.Bxmin)
    FileHandle.write(line)

    line = "Bymin = {} ".format(parvar.Bymin)
    FileHandle.write(line)

    line = "Bxmax = {} ".format(parvar.Bxmax)
    FileHandle.write(line)

    line = "Bymax = {} \n".format(parvar.Bymax)
    FileHandle.write(line)

    line = "Split = {} \n".format(parvar.Split)
    FileHandle.write(line)

#    line = "AutoSatRegion = {} \n".format(parvar.AutoSatRegion)
#    FileHandle.write(line)

    line = "SatRegionScale = {} \n".format(parvar.SatRegionScale)
    FileHandle.write(line)

    line = "Ds9SatReg = {} \n".format(parvar.Ds9SatReg)
    FileHandle.write(line)

    line = "FileOut = {} \n".format(parvar.FileOut)
    FileHandle.write(line)

    line = "SegFile = {} \n".format(parvar.SegFile)
    FileHandle.write(line)

    line = "SkyFile = {} \n".format(parvar.SkyFile)
    FileHandle.write(line)

    line = "PixPrefix = {} \n".format(parvar.PixPrefix)
    FileHandle.write(line)

    line = "Ds9OutName = {} \n".format(parvar.Ds9OutName)
    FileHandle.write(line)

    line = "Ds9OutNum = {} \n".format(parvar.Ds9OutNum)
    FileHandle.write(line)

    line = "Ds9FitReg = {} \n".format(parvar.Ds9FitReg)
    FileHandle.write(line)

    line = "BoxOut = {} \n".format(parvar.BoxOut)
    FileHandle.write(line)

    line = "BoxSkyOut = {} \n".format(parvar.BoxSkyOut)
    FileHandle.write(line)

    line = "SexSort = {} \n".format(parvar.SexSort)
    FileHandle.write(line)

    line = "SexArSort = {} \n".format(parvar.SexArSort)
    FileHandle.write(line)

    line = "Erase = {} \n".format(parvar.Erase)
    FileHandle.write(line)

    line = "Nice = {} \n".format(parvar.Nice)
    FileHandle.write(line)

    line = "Overwrite = {} \n".format(parvar.Overwrite)
    FileHandle.write(line)

    line = "Execute = {} \n".format(parvar.Execute)
    FileHandle.write(line)

    line = "TempDir = {} \n".format(parvar.TempDir)
    FileHandle.write(line)

#    line = "PixDir = {} \n".format(parvar.PixDir)
#    FileHandle.write(line)

#    line = "MaskPixDir = {} \n".format(parvar.MaskPixDir)
#    FileHandle.write(line)

    line = "SkyDir = {} \n".format(parvar.SkyDir)
    FileHandle.write(line)

    line = "OutputDir = {} \n".format(parvar.OutputDir)
    FileHandle.write(line)

    line = "InputDir = {} \n".format(parvar.InputDir)
    FileHandle.write(line)

    line = "NCol = {} \n".format(parvar.NCol)
    FileHandle.write(line)

    line = "NRow = {} \n".format(parvar.NRow)
    FileHandle.write(line)

    line = "Buffer = {} \n".format(parvar.Buffer)
    FileHandle.write(line)

    line = "Total = {} \n".format(parvar.Total)
    FileHandle.write(line)

    line = "LogFile = {} \n".format(parvar.LogFile)
    FileHandle.write(line)

    line = "BtFile = {} \n".format(parvar.BtFile)
    FileHandle.write(line)

    line = "FitsFile = {} \n".format(parvar.FitsFile)
    FileHandle.write(line)

    line = "SimulFit = {} \n".format(parvar.SimulFit)
    FileHandle.write(line)

    line = "Contrast = {} \n".format(parvar.Contrast)
    FileHandle.write(line)

    line = "Bias = {} \n".format(parvar.Bias)
    FileHandle.write(line)

# ellipsectgalfit

    line = "ranx = {} \n".format(parvar.ranx)
    FileHandle.write(line)

    line = "rany = {} \n".format(parvar.rany)
    FileHandle.write(line)

    line = "plotsub = {} \n".format(parvar.flagsub)
    FileHandle.write(line)

    line = "dpi = {} \n".format(parvar.dpi)
    FileHandle.write(line)

    line = "outsb = {} \n".format(parvar.flagout)
    FileHandle.write(line)

    line = "plotpix = {} \n".format(parvar.flagpix)
    FileHandle.write(line)






    line = "  \n\n"
    FileHandle.write(line)


    return True


####
##  functions not checked, caution!
####

####
# Deprecated
###
def ReadFitlog(output):
    "Read fit.log file of GALFIT output and create a new file in columns format"
    # Read fit.log file of GALFIT output
    # and create a new file in columns format

    outflag = True

    inputf = "fit.log"

    FOUT = open(output, "w")

    with open(inputf) as FITLOG:
        try:
            # All lines including the blank ones
            lines = (line.rstrip() for line in FITLOG)
            lines = (line.split('#', 1)[0] for line in lines)  # remove comments
            # remove lines containing only comments
            lines = (line.rstrip() for line in lines)
            lines = (line for line in lines if line)  # Non-blank lines

            lines = list(lines)
            index = 0

            function="sersic"

            while index < len(lines):

                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "Input":     # input image

#Input image     : tempfits/temp-2-1.fits[1:762,1:474]
#Init. par. file : obj2099
#Restart file    : galfit.48
#Output image    : A1213.fits-111701.79-p290557.3-2099-out.fits


                    inputimage=tmp[3]

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    initfile=tmp[4]

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    restartfile=tmp[3]


                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    outim = tmp[3]


                    while (tmp[0] != "Chi^2"):

                        index += 1
                        line = lines[index]
                        (tmp) = line.split()


                        if (tmp[0] == "sersic"):

# sersic    : (235.08, 132.86)   14.75   6.09    4.84    0.86   -57.85
#               ( 0.00,   0.00)   0.01     0.09    0.09   0.00   1.21

                            function=tmp[0]

                            (posxser,posyser,magser,reser,nser,axiser,paser)=GetSerVars(line)

                            index += 1
                            line = lines[index]

                            (erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser)=GetSerErrVars(line)


#### print line to Output
                            statout = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage,initfile,restartfile,outim,function,posxser,posyser, magser, reser, nser, axiser, paser, erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser)
                            FOUT.write(statout)
##################



                        if (tmp[0] == "expdisk"):   ## add x y position
# expdisk   : ( {482.78}, {338.62})  16.52      2.40    0.97    22.60
#             (   {0.00},   {0.00})   0.00      0.02    0.00     6.34


                            function=tmp[0]
                            #(tmp) = line.split()

                            (posxexp,posyexp,magexp,rsexp,axisexp,paexp)=GetExpVars(line)

                        #errs
                            index += 1
                            line = lines[index]
#                       (tmp) = line.split()

                            (erposxexp,erposyexp,ermagexp,errsexp,eraxisexp,erpaexp)=GetExpErrVars(line)

#### print line to Output
                            line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, posxexp, posyexp, magexp, rsexp, 1, axisexp, paexp, erposxexp, erposyexp, ermagexp, errsexp, 0, eraxisexp, erpaexp)
                            FOUT.write(line)
##################


                        if (tmp[0] == "gaussian"):


                            function=tmp[0]

                            (posxbar,posybar,magbar, fwhm, arbar,pabar)=GetGaussVars(line)

                            index += 1
                            line = lines[index]

                            (erposxbar,erposybar,ermagbar,erfwhm,erarbar,erpabar)=GetGaussErrVars(line)

#### print line to Output
                            statout = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, posxbar, posybar, magbar, fwhm, 0.5, arbar, pabar, erposxbar ,erposybar , ermagbar,erfwhm, 0, erarbar, erpabar)
                            FOUT.write(statout)
##################


                        if( tmp[0] == "sky"):

                            function=tmp[0]
                            sky=tmp[4]

                            sky = sky.strip("[")
                            sky = sky.strip("]")
                            sky = sky.strip("{")
                            sky = sky.strip("}")

                            index += 1
                            line = lines[index]
                            (tmp) = line.split()

                            ersky=tmp[0]

                            ersky = ersky.strip("[")
                            ersky = ersky.strip("]")
                            ersky = ersky.strip("{")
                            ersky = ersky.strip("}")

                            line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim,function, 0, 0, sky, 0, 0, 0, 0, 0, 0, ersky, 0, 0, 0, 0)
                            FOUT.write(line)


                        if(tmp[0] == "Chi^2"):

                            chi = tmp[2]
                            chi = chi.strip(",")

                            ndof= tmp[5]

                            index += 1
                            line  = lines[index]
                            (tmp) = line.split()

                            chinu=tmp[2]

                            function=tmp[0]


                            line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, 0, 0, chinu, chi, ndof, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
                            FOUT.write(line)


                index += 1
        except:
            print("Can't open fit.log \n")
            outflag = False


    FOUT.close()

    return outflag


def ReadFitlog2(output):
    "Read fit.log file of GALFIT output and create a new file in columns format"
    # Read fit.log file of GALFIT output
    # and create a new file in columns format

    outflag = True

    inputf = "fit.log"

    FOUT = open(output, "w")

    FITLOG = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in FITLOG)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    function="sersic"

    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()
        if tmp[0] == "Input":     # input image

#Input image     : tempfits/temp-2-1.fits[1:762,1:474]
#Init. par. file : obj2099
#Restart file    : galfit.48
#Output image    : A1213.fits-111701.79-p290557.3-2099-out.fits


            inputimage=tmp[3]

            index += 1
            line = lines[index]
            (tmp) = line.split()

            initfile=tmp[4]

            index += 1
            line = lines[index]
            (tmp) = line.split()

            restartfile=tmp[3]


            index += 1
            line = lines[index]
            (tmp) = line.split()

            outim = tmp[3]


            while (tmp[0] != "Chi^2/nu"):

                index += 1
                line = lines[index]
                (tmp) = line.split()


                if (tmp[0] == "sersic"):

# sersic    : (235.08, 132.86)   14.75   6.09    4.84    0.86   -57.85
#               ( 0.00,   0.00)   0.01     0.09    0.09   0.00   1.21

                    function=tmp[0]

                    (posxser,posyser,magser,reser,nser,axiser,paser)=GetSerVars(line)

                    index += 1
                    line = lines[index]

                    (erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser)=GetSerErrVars(line)


#### print line to Output
                    statout = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage,initfile,restartfile,outim,function,posxser,posyser, magser, reser, nser, axiser, paser, erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser)
                    FOUT.write(statout)
##################

                if (tmp[0] == "expdisk"):   ## add x y position
# expdisk   : ( {482.78}, {338.62})  16.52      2.40    0.97    22.60
#             (   {0.00},   {0.00})   0.00      0.02    0.00     6.34

                    function=tmp[0]
                            #(tmp) = line.split()

                    (posxexp,posyexp,magexp,rsexp,axisexp,paexp)=GetExpVars(line)

                    #errs
                    index += 1
                    line = lines[index]
#                   (tmp) = line.split()

                    (erposxexp,erposyexp,ermagexp,errsexp,eraxisexp,erpaexp)=GetExpErrVars(line)

#### print line to Output
                    line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, posxexp, posyexp, magexp, rsexp, 1, axisexp, paexp, erposxexp, erposyexp, ermagexp, errsexp, 0, eraxisexp, erpaexp)
                    FOUT.write(line)
##################


                if (tmp[0] == "gaussian"):

                    function=tmp[0]

                    (posxbar,posybar,magbar, fwhm, arbar,pabar)=GetGaussVars(line)

                    index += 1
                    line = lines[index]

                    (erposxbar,erposybar,ermagbar,erfwhm,erarbar,erpabar)=GetGaussErrVars(line)

#### print line to Output
                    statout = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, posxbar, posybar, magbar, fwhm, 0.5, arbar, pabar, erposxbar ,erposybar , ermagbar,erfwhm, 0, erarbar, erpabar)
                    FOUT.write(statout)
##################

                if(tmp[0] == "sky"):

                    function=tmp[0]
                    sky=tmp[4]

                    sky = sky.strip("[")
                    sky = sky.strip("]")
                    sky = sky.strip("{")
                    sky = sky.strip("}")

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    ersky=tmp[0]

                    ersky = ersky.strip("[")
                    ersky = ersky.strip("]")
                    ersky = ersky.strip("{")
                    ersky = ersky.strip("}")

                    line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim,function, 0, 0, sky, 0, 0, 0, 0, 0, 0, ersky, 0, 0, 0, 0)
                    FOUT.write(line)


                if(tmp[0] == "Chi^2"):

                    chi = tmp[2]
                    chi = chi.strip(",")

                    ndof= tmp[5]

                    index += 1
                    line  = lines[index]
                    (tmp) = line.split()

                    chinu=tmp[2]

                    function=tmp[0]


                    line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(inputimage, initfile, restartfile, outim, function, 0, 0, chinu, chi, ndof, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
                    FOUT.write(line)


        index += 1

    FOUT.close()
    FITLOG.close()


    return outflag



def GetSerVars(statement):

    (left,right) = statement.split(")")
    (trash,pos)  = left.split("(")
    (posxser,posyser) = pos.split()

    posxser = posxser.strip(",")  # equivalent in python

    (words) = right.split()

    magser  = words[0]
    reser   = words[1]
    nser    = words[2]
    axiser  = words[3]
    paser   = words[4]


# cleaning variables from [] {}

    posxser = posxser.strip("[")
    posxser = posxser.strip("]")
    posxser = posxser.strip("{")
    posxser = posxser.strip("}")
    posyser = posyser.strip("[")
    posyser = posyser.strip("]")
    posyser = posyser.strip("{")
    posyser = posyser.strip("}")

    magser = magser.strip("[")
    magser = magser.strip("]")
    magser = magser.strip("{")
    magser = magser.strip("}")

    reser = reser.strip("[")
    reser = reser.strip("]")
    reser = reser.strip("{")
    reser = reser.strip("}")

    nser = nser.strip("[")
    nser = nser.strip("]")
    nser = nser.strip("{")
    nser = nser.strip("}")

    axiser = axiser.strip("[")
    axiser = axiser.strip("]")
    axiser = axiser.strip("{")
    axiser = axiser.strip("}")

    paser = paser.strip("[")
    paser = paser.strip("]")
    paser = paser.strip("{")
    paser = paser.strip("}")


    return (posxser,posyser,magser,reser,nser,axiser,paser)


def GetSerErrVars(statement):


    (left,right) = statement.split(")")
    (trash,pos)  = left.split("(")
    (erposxser,erposyser) = pos.split()
    erposxser = erposxser.strip(",")  # equivalent in python
    (words) = right.split()


    ermagser   = words[0]
    erreser    = words[1]
    ernser     = words[2]
    eraxiser   = words[3]
    erpaser    = words[4]



#    erposxser = erposxser.strip(",")  # equivalent in python
#    erposxser = erposxser.strip("(")  # equivalent in python

#    erposyser = erposyser.strip(")")

    # cleaning err variables

    erposxser = erposxser.strip("[")
    erposxser = erposxser.strip("]")
    erposxser = erposxser.strip("{")
    erposxser = erposxser.strip("}")
    erposyser = erposyser.strip("[")
    erposyser = erposyser.strip("]")
    erposyser = erposyser.strip("{")
    erposyser = erposyser.strip("}")

    ermagser = ermagser.strip("[")
    ermagser = ermagser.strip("]")
    ermagser = ermagser.strip("{")
    ermagser = ermagser.strip("}")

    erreser = erreser.strip("[")
    erreser = erreser.strip("]")
    erreser = erreser.strip("{")
    erreser = erreser.strip("}")

    ernser = ernser.strip("[")
    ernser = ernser.strip("]")
    ernser = ernser.strip("{")
    ernser = ernser.strip("}")

    eraxiser = eraxiser.strip("[")
    eraxiser = eraxiser.strip("]")
    eraxiser = eraxiser.strip("{")
    eraxiser = eraxiser.strip("}")

    erpaser = erpaser.strip("[")
    erpaser = erpaser.strip("]")
    erpaser = erpaser.strip("{")
    erpaser = erpaser.strip("}")


    return (erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser)


def GetExpVars(statement):

    words = statement.split()

    posxexp = words[3]
    posyexp = words[4]
    magexp  = words[5]
    rsexp   = words[6]
    axisexp = words[7]
    paexp   = words[8]

    posxexp = posxexp.strip(",")  # equivalent in python
    posyexp = posyexp.strip(")")

    posxexp = posxexp.strip("[")
    posxexp = posxexp.strip("]")
    posxexp = posxexp.strip("{")
    posxexp = posxexp.strip("}")
    posyexp = posyexp.strip("[")
    posyexp = posyexp.strip("]")
    posyexp = posyexp.strip("{")
    posyexp = posyexp.strip("}")

    magexp = magexp.strip("[")
    magexp = magexp.strip("]")
    magexp = magexp.strip("{")
    magexp = magexp.strip("}")

    rsexp = rsexp.strip("[")
    rsexp = rsexp.strip("]")
    rsexp = rsexp.strip("{")
    rsexp = rsexp.strip("}")

    axisexp = axisexp.strip("[")
    axisexp = axisexp.strip("]")
    axisexp = axisexp.strip("{")
    axisexp = axisexp.strip("}")

    paexp = paexp.strip("[")
    paexp = paexp.strip("]")
    paexp = paexp.strip("{")
    paexp = paexp.strip("}")

    return (posxexp,posyexp,magexp,rsexp,axisexp,paexp)


def GetExpErrVars(statement):

#    words = statement.split()
    partrash,stats= statement.split("(")
    words = stats.split()


    erposxexp     = words[0]
    erposyexp     = words[1]
    ermagexp      = words[2]
    errsexp       = words[3]
    eraxisexp     = words[4]
    erpaexp       = words[5]

    erposxexp = erposxexp.strip(",")  # equivalent in python
    erposyexp = erposyexp.strip(")")

    erposxexp = erposxexp.strip("[")
    erposxexp = erposxexp.strip("]")
    erposxexp = erposxexp.strip("{")
    erposxexp = erposxexp.strip("}")
    erposyexp = erposyexp.strip("[")
    erposyexp = erposyexp.strip("]")
    erposyexp = erposyexp.strip("{")
    erposyexp = erposyexp.strip("}")

    ermagexp = ermagexp.strip("[")
    ermagexp = ermagexp.strip("]")
    ermagexp = ermagexp.strip("{")
    ermagexp = ermagexp.strip("}")

    errsexp = errsexp.strip("[")
    errsexp = errsexp.strip("]")
    errsexp = errsexp.strip("{")
    errsexp = errsexp.strip("}")


    eraxisexp = eraxisexp.strip("[")
    eraxisexp = eraxisexp.strip("]")
    eraxisexp = eraxisexp.strip("{")
    eraxisexp = eraxisexp.strip("}")

    erpaexp = erpaexp.strip("[")
    erpaexp = erpaexp.strip("]")
    erpaexp = erpaexp.strip("{")
    erpaexp = erpaexp.strip("}")


    return (erposxexp,erposyexp,ermagexp,errsexp,eraxisexp,erpaexp)



def GetGaussVars(statement):

#   gaussian  : ( 261.93, 350.00)   15.49    12.56    0.47    36.77
#              ( 0.01, 0.01)   0.00     0.04     0.00    0.12

    words = statement.split()

    posxbar = words[3]
    posybar = words[4]

    posxbar = posxbar.strip(",")
#    posxbar = posxbar.strip("(")

    posybar = posybar.strip(")")

    magbar = words[5]
    fwhm   = words[6]
    arbar  = words[7]
    pabar  = words[8]

    posxbar  =  posxbar.strip("[")
    posxbar  =  posxbar.strip("]")
    posxbar  =  posxbar.strip("{")
    posxbar  =  posxbar.strip("}")

    posybar  =  posybar.strip("[")
    posybar  =  posybar.strip("]")
    posybar  =  posybar.strip("{")
    posybar  =  posybar.strip("}")

    magbar = magbar.strip("[")
    magbar = magbar.strip("]")
    magbar = magbar.strip("{")
    magbar = magbar.strip("}")

    fwhm = fwhm.strip("[")
    fwhm = fwhm.strip("]")
    fwhm = fwhm.strip("{")
    fwhm = fwhm.strip("}")

    arbar = arbar.strip("[")
    arbar = arbar.strip("]")
    arbar = arbar.strip("{")
    arbar = arbar.strip("}")

    pabar = pabar.strip("[")
    pabar = pabar.strip("]")
    pabar = pabar.strip("{")
    pabar = pabar.strip("}")

    return (posxbar,posybar,magbar,fwhm,arbar,pabar)



def GetGaussErrVars(statement):

#   gaussian  : (261.93, 350.00)   15.49    12.56    0.47    36.77
#              (0.01, 0.01)   0.00     0.04     0.00    0.12

    words = statement.split()

    erposxbar=words[0]
    erposybar=words[1]

    erposxbar = erposxbar.strip(",")  # equivalent in python
    erposybar = erposxbar.strip(")")

    ermagbar = words[2]
    erfwhm   = words[3]
    erarbar  = words[4]
    erpabar  = words[5]


    erposxbar  =  erposxbar.strip("[")
    erposxbar  =  erposxbar.strip("]")
    erposxbar  =  erposxbar.strip("{")
    erposxbar  =  erposxbar.strip("}")

    erposybar  =  erposybar.strip("[")
    erposybar  =  erposybar.strip("]")
    erposybar  =  erposybar.strip("{")
    erposybar  =  erposybar.strip("}")

    ermagbar = ermagbar.strip("[")
    ermagbar = ermagbar.strip("]")
    ermagbar = ermagbar.strip("{")
    ermagbar = ermagbar.strip("}")

    erfwhm = erfwhm.strip("[")
    erfwhm = erfwhm.strip("]")
    erfwhm = erfwhm.strip("{")
    erfwhm = erfwhm.strip("}")

    erarbar = erarbar.strip("[")
    erarbar = erarbar.strip("]")
    erarbar = erarbar.strip("{")
    erarbar = erarbar.strip("}")

    erpabar = erpabar.strip("[")
    erpabar = erpabar.strip("]")
    erpabar = erpabar.strip("{")
    erpabar = erpabar.strip("}")


    return (erposxbar,erposybar,ermagbar,erfwhm,erarbar,erpabar)



def MakeOutput(inputf, parvar, obj):
    "Read output file of ReadFitlog and computes B/T and other parameters and creates final file"

#    consflag=1
#    inflag=True
    countimg=1
    countgal=0
    logbflag=True
    logdflag=True
    skyflag=True

    errflag=False

#    mainflag=False

    sigma=40


## putting some limit to errors


#    maglim =5
#    relim=100
#    nlim=5
#    palim=50
#    axiserlim=0.3
####################

#    magexplim=5
#    rslim=100
#    paexplim=50
#    axisexplim=0.3


# check if  the error files exist in actual directory

# REMOVED

#    if(!(-e "bdsigmas.txt") )
#    {
#	print "can't found bdsigmas \n";
#    }
#    if(!(-e "bsigmas.txt") )
#    {
#	print "can't found bsigmas \n";
#    }
############################################
#############################################
# NEW value of these flags
# Meaning of the fitflag:

#	0 = good fit
#	1 = ran into constraints during the fit
#	2 = GALFIT found some error in one (or several) parameters
#	4 = nan errors in fit.log
#   8 = galfit crashed

# refit = True # it needs to rerun GALFIT on it

###############################################
###############################################

### Constraints file


    (consmagmin1,consmagmax1,consremin1,consremax1,consnmin1,consnmax1) = ReadConstraints(parvar.ConsFile,1)
    (consmagmin2,consmagmax2,consremin2,consremax2,consnmin2,consnmax2) = ReadConstraints(parvar.ConsFile,2)


    with open(inputf) as INCAT:
#        try:

        ds9images = parvar.Ds9OutName+"-"+ str(countimg)

        OUTBT    = open(parvar.BtFile,"w")
        OUTDS9   = open(ds9images,"w")
        OUTINFO  = open("inforegion","w")
        OUTFLAGS = open("objflags","w")


        lineds9 = "#!/bin/bash \n"
        OUTDS9.write(lineds9)

        lineds9 = "ds9 -tile -cmap grey -cmap value {} {} -invert -log -zmax ".format(parvar.Contrast,parvar.Bias)
        OUTDS9.write(lineds9)

        lines = (line.rstrip() for line in INCAT) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        lines = (line.rstrip() for line in lines)           # remove lines containing only comments
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0
        while index < len(lines):

            line = lines[index]
            tmp=line.split()

            mainflag=False
            btflag = False
            barflag=False
            refit=False

            fitflag=0
            distmax=5
            distmaxdisk=5

            magbar=99
            posxbar=0
            posybar=0

            pixrebar=0
            nbar=0
            axisbar=0
            pabar=0
            erposxbar=0
            erposybar=0
            ermagbar=0
            erpixrebar=0
            ernbar=0
            eraxisbar=0
            erpabar=0

            skyflag=True
            logbflag=True
            logdflag=True
            errflag=False

            while(tmp[4] != "Chi^2/nu"):

                if (mainflag == False):

                    (inputimage,initfile,restartfile,outim,function,posxser,posyser,magser,pixreser,nser,axiser,paser,erposxser,erposyser,ermagser,erreser,ernser,eraxiser,erpaser) = line.split()

                    mainflag=True



## checking  errors in fitting
                    if posxser.find("*") != -1:
                        errflag = True
                        posxser = posxser.strip("*")
                        erposxser = erposxser.strip("*")

                    if posyser.find("*") != -1:
                        errflag = True
                        posyser = posyser.strip("*")
                        erposyser = erposyser.strip("*")


                    if magser.find("*") != -1:
                        errflag = True
                        magser = magser.strip("*")
                        ermagser = ermagser.strip("*")

                    if pixreser.find("*") != -1:
                        errflag = True
                        pixreser = pixreser.strip("*")
                        erpixreser = erpixreser.strip("*")

                    if nser.find("*") != -1:
                        errflag = True
                        nser = nser.strip("*")
                        ernser = ernser.strip("*")


                    if axiser.find("*") != -1:
                        errflag = True
                        axiser = axiser.strip("*")
                        eraxiser = eraxiser.strip("*")

                    if paser.find("*") != -1:
                        errflag = True
                        paser = paser.strip("*")
                        erpaser = erpaser.strip("*")



                    (inputimage,coorde)=inputimage.split("[")
                    (coorx,coory)=coorde.split(",")
                    (xmin,xmax)=coorx.split(":")
                    (ymin,ymax)=coory.split(":")
                    ymax=ymax.strip("]")


                    objx=initfile
                    (trash,pppnum)=objx.split("-")

                    pppnum=int(pppnum)


                    restart=restartfile


                    (cluster,cordra,corddect,pppnum2,trash) = outim.split("-")

                    (dirtrash,cluster) = cluster.split("/")  #filter removing directory

                    corddec = corddect[1:]

                    sign = corddect[0:1]

                    rasec = cordra[4:]

                    decsec = corddec[4:]

                    ramin=cordra[2:4]

                    decmin=corddec[2:4]

                    rax=cordra[0:2]
                    decx=corddec[0:2]


                    ra = rax+":"+ramin+":"+rasec

                    if(sign == "p"):
                        dec = decx+":"+decmin+":"+decsec
                    elif (sign == "m"):
                        sdecx = (-1)*int(decx)  #change to int instead of float
                        dec = str(sdecx)+":"+decmin+":"+decsec

                    ra2 = float(rax)+ ( float(ramin) + float(rasec) / 60) / 60

                    if(sign == "p"):
                        dec2 = float(decx) +  (float(decmin)  + float(decsec) / 60 ) / 60

                    elif(sign == 'm'):
                        dec2 = float(decx) +  (float(decmin)  + float(decsec) / 60 ) / 60
                        dec2 = (-1)*dec2

                    nser=float(nser)

                    xmin = int(xmin)
                    xmax = int(xmax)
                    ymin = int(ymin)
                    ymax = int(ymax)

                    posxser  = float(posxser)
                    posyser  = float(posyser)
                    magser   = float(magser)
                    pixreser = float(pixreser)
                    nser = float(nser)
                    axiser = float(axiser)
                    paser = float(paser)
                    erposxser = float(erposxser)
                    erposyser = float(erposyser)
                    ermagser = float(ermagser)
                    erreser = float(erreser)
                    ernser = float(ernser)
                    eraxiser = float(eraxiser)
                    erpaser = float(erpaser)


                    kser = GetK(nser)

                    xser  = posxser - xmin
                    yser  = posyser - ymin


####begin
                    if (parvar.FitFunc == "BD"):
                        index += 1
                        line = lines[index]


                        (inputimageexp,initfileexp,restartfileexp,outimexp,functionexp,posxexp,posyexp,magexp,pixrsexp,nexp,axisexp,paexp,erposxexp,erposyexp,ermagexp,errsexp,ernexp,eraxisexp,erpaexp)  = line.split()

## removing *s

                        posxexp   = posxexp.strip("*")
                        erposxexp = erposxexp.strip("*")

                        posyexp   = posyexp.strip("*")
                        erposyexp = erposyexp.strip("*")

                        posxexp  = float(posxexp)
                        posyexp  = float(posyexp)
                        erposxexp= float(erposxexp)
                        erposyexp = float(erposyexp)


                        dist = np.sqrt ( ( posxser - posxexp )**2 + ( posyser - posyexp )**2 )


                        if (dist < distmaxdisk):

                            btflag=True
                            if magexp.find("*") != -1:
                                errflag = True
                                magexp = magexp.strip("*")
                                ermagexp = ermagexp.strip("*")

                            if pixrsexp.find("*") != -1:
                                errflag = True
                                pixrsexp = pixrsexp.strip("*")
                                errsexp  = errsexp.strip("*")


                            if axisexp.find("*") != -1:
                                errflag = True
                                axisexp = axisexp.strip("*")
                                eraxisexp  = eraxisexp.strip("*")


                            if paexp.find("*") != -1:
                                errflag = True
                                paexp = paexp.strip("*")
                                erpaexp  = erpaexp.strip("*")


                            magexp   = float(magexp)
                            pixrsexp = float(pixrsexp)
                            axisexp = float(axisexp)
                            paexp = float(paexp)
                            ermagexp = float(ermagexp)
                            errsexp = float(errsexp)
                            eraxisexp = float(eraxisexp)
                            erpaexp = float(erpaexp)



                            pixreexp= (1.678347) * pixrsexp
                            kexp=1.678347
### end

                    else:
                        magexp=99
                        reexp=0
                        rsexp=0
                        pixreexp=0
                        pixrsexp=0
                        errsexp=0
                        posxexp=0
                        posyexp=0
                        axisexp=0
                        eraxisexp=0
                        paexp=0
                        erpaexp=0
                        meanmeexp=99
                        meexp=99
                        msexp=99
                        ermeanmeexp=99
                        ermeexp=99
                        ermsexp=99


                        ermagexp=99
                        erposxexp=0
                        erposyexp=0

                        kexp=1.678347


                    if(btflag == True):
                        fluxbulge =  10 **((-magser)/2.5 )
                        fluxdisk  =  10 **((-magexp)/2.5 )
                        magtotal  = (-2.5) * np.log10 (fluxbulge + fluxdisk)
                        bulgetotal= BulgeTotal(magser,magexp)

                    elif(btflag == False):
                        magtotal=magser
                        bulgetotal=1
                        index -= 1


## begin
                    if (parvar.FitFunc == "BD"):
                        index += 1
                        line = lines[index]
                        (tmp) = line.split()


                        tmpx = tmp[5]
                        tmpy = tmp[6]

                        tmpx   = tmpx.strip("*")
                        tmpy   = tmpy.strip("*")


                        dist = np.sqrt ( ( posxser - float(tmpx) )**2 + ( posyser - float(tmpy) )**2 )

                        if (dist < distmax and barflag == False and btflag == True):

                            barflag = True

                            posxbar=float(tmp[5])
                            posybar=float(tmp[6])
                            magbar=float(tmp[7])
                            pixrebar=float(tmp[8])
                            nbar=float(tmp[9])
                            axisbar=float(tmp[10])
                            pabar=float(tmp[11])
                            erposxbar=float(tmp[12])
                            erposybar=float(tmp[13])
                            ermagbar=float(tmp[14])
                            erpixrebar=float(tmp[15])
                            ernbar=float(tmp[16])
                            eraxisbar=float(tmp[17])
                            erpabar=float(tmp[18])
## end

                if (tmp[4] == "sky"):
                    sky   = float(tmp[7])
                    ersky = float(tmp[14])

                index += 1
                line  = lines[index]
                (tmp) = line.split()

            chinu=float(tmp[7])
            if (barflag == True):
                fluxbar =  10 **((-magbar)/2.5)
                magserbar = (-2.5) * np.log10 (fluxbulge + fluxbar)
                magtotal  = (-2.5) * np.log10 (fluxbulge + fluxdisk + fluxbar)

                bulgetotal= BulgeTotal(magserbar,magexp)


            reser=pixreser*parvar.PlateScale
            erreser=erreser*parvar.PlateScale

            rsexp=pixrsexp*parvar.PlateScale
            errsexp=errsexp*parvar.PlateScale

            reexp=pixreexp*parvar.PlateScale


            if (reser > 0):

                meanmeser = magser + 2.5 * np.log10(2 * np.pi * reser * reser)

                fn = (( axiser * nser * np.exp( kser)) / (kser ** (2 * nser )) ) * ( np.exp(scipy.special.gammaln(2*nser)) )

                meser = meanmeser +  2.5 * np.log10( fn )

                ermeanmeser=np.sqrt(ermagser**2 + (erreser**2)*(5*np.log10(np.e)/reser)**2 )

# ATENTION check how to include axiser in devfn
###   computing ermeser
                dudn  =  np.e*scipy.special.gamma(2*nser) + nser * np.e**kser * scipy.special.gamma(2*nser) * scipy.special.polygamma(0, 2*nser)
                dvdn  =  2*np.log(kser)*(kser ** (2 * nser ))
                devfn = (kser**(2 * nser) * dudn - (nser*np.e**kser)*scipy.special.gamma(2*nser)*dvdn) / ((kser ** (2 * nser ))**2)


                pardevfn=(2.5/np.log(10))*(devfn/fn)

                ermeser= np.sqrt(ermeanmeser**2 + (pardevfn**2 )*(ernser**2) )

###
            else:
                logbflag=False

            if (btflag == True):
                if (reexp > 0):

                    meanmeexp = magexp  + 2.5 * np.log10(2 * np.pi * reexp * reexp)

                    fnexp= (axisexp * np.exp( kexp) / (kexp ** (2 )) ) * ( np.exp(scipy.special.gammaln(2)) )
                    meexp = meanmeexp +  2.5 * np.log10( fnexp )

                    msexp = meexp - 1.822

                    ermeanmeexp=np.sqrt(ermagexp**2 + (erreexp**2)*(5*np.log10(np.e)/reexp)**2 )

                    ermeexp=ermeanmeexp
                    ermsexp=ermeexp


                else:
                    logdflag=False

            else:

                magexp=99
                ermagexp=0
                reexp=0
                rsexp=0
                pixreexp=0
                pixrsexp=0
                errsexp=0
                posxexp=0
                posyexp=0
                axisexp=0
                eraxisexp=0
                paexp=0
                erpaexp=0
                meanmeexp=99
                meexp=99




#            filepix= "{}/{}-{}".format(parvar.PixDir,parvar.PixPrefix,pppnum) #delete, not needed

#            imres = "{}/{}".format(parvar.OutputDir,outim)
            imres = "{}".format(outim)



            imsig = "{}/sigma-{}.fits".format(parvar.MaskDir,pppnum)

#            immask = "{}/mask-{}".format(parvar.MaskDir,pppnum)

            immask = parvar.SegFile

## averaged weight radius

            disktotal = 1 - bulgetotal

            bdradius = pixreser * bulgetotal + pixreexp * disktotal

            rmin=2


            maskind = obj.Num == pppnum
            ind=np.where(maskind == True)
            indx=ind[0][0]

            sexxser = obj.XPos[indx]
            sexyser = obj.YPos[indx]

            sexxmin = obj.XMin[indx]
            sexxmax = obj.XMax[indx]
            sexymin = obj.YMin[indx]
            sexymax = obj.YMax[indx]

            goffx =  obj.gOFFX[indx]
            goffy =  obj.gOFFY[indx]


            (tidal,objchinu,bump,snr,ndof) = Tidal(imres, imsig, immask, pppnum, sexxser, sexyser, sexxmin, sexxmax, sexymin, sexymax, goffx, goffy, rmin, sky, btflag)


#               def Tidal(imgout,imsig,immask,num,xser,yser,xlo,xhi,ylo,yhi,goffx,goffy,rmin,sky,btflag):
#               imgout = galfit output image
#               imsig = sigma image (The one created by GALFIT)
#               immask = Mask image: The original
#               num = object number
#               xser, yser = Sextractor object position on the big image
#               xlo, xhi, ylo, yhi = coordinate of object of the big image
#               goffx, goffy = coordinate tranformation from big to galfit image
#               rmin = minimum radius to compute Tidal and Bumpiness (to avoid PSF Mismatch)
#               sky = Background value
#               btflag = is bulge/disk or sersic component?


            if (parvar.FitFunc == "BD"):
                qarg= axisexp
                parg= paexp
            else:
                qarg= axiser
                parg= paser

   
            mgetools.EllipSec(restart,qarg,parg,parvar.flagsub,parvar.flagpix,parvar.ranx,parvar.rany,parvar.dpi,parvar.flagout)


            if ( magser == consmagmin1  or magser == consmagmax1 ):

                fitflag=fitflag+1
                refit=True

            elif ( pixreser == consremin1 or pixreser == consremax1 ):
                fitflag=fitflag+1
                refit=True

            elif ( nser == consnmin1  or nser == consnmax1 ):
                fitflag=fitflag+1
                refit=True


            elif (btflag == True):
                if ( magexp == consmagmin2 or magexp == consmagmax2):
                    fitflag=fitflag+1
                    refit=True


                elif (pixrsexp == consremin2 or pixrsexp == consremax2):
                    fitflag=fitflag+1
                    refit=True


            if (errflag==True):
                fitflag=fitflag+2
                refit=True




# big errors? if so, then flag it  only nan errors

            if ((ermagser == "nan") or (erreser == "nan") or (ernser == "nan") or (erpaser == "nan") or (eraxiser == "nan") or (ermagser == "-nan") or (erreser == "-nan") or (ernser == "-nan") or (erpaser == "-nan") or (eraxiser == "-nan")):
                fitflag=fitflag+4
                refit=True

            elif (btflag == True):
                if((ermagexp == "nan") or (errsexp == "nan") or (erpaexp == "nan") or (eraxisexp == "nan") or (ermagexp == "-nan") or (errsexp == "-nan") or (erpaexp == "-nan") or (eraxisexp == "-nan")):
                    fitflag=fitflag+4
                    refit=True


            if (sky > 0):
                magsky =  (-2.5) * np.log10 ( sky / parvar.ExpTime ) + 2.5 * np.log10 (parvar.PlateScale**2) + parvar.MagZpt

                meanmesky = meanmeser - magsky

            else:
                skyflag=False


            fileppp = "{}/Reg-{}.reg".format(parvar.OutputDir,pppnum)

            if (btflag == True):
                dx=posxser-posxexp
                dy=posyser-posyexp
            else:
                dx=0
                dy=0


            countgal+=1


            lineds9 = "{}[1] -regions {} {}[2] -regions {} {}[3] -regions {} ".format(imres,fileppp,imres,fileppp,imres,fileppp)
            OUTDS9.write(lineds9)


            lineds9="{} {} {} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {} \n".format(fileppp,pppnum,xser,yser,dx,dy,magtotal,pixreser,nser,bulgetotal,chinu,snr,axiser,axisexp,pixrsexp,tidal,objchinu,bump,paser,paexp,fitflag)
            OUTINFO.write(lineds9)

            lineds9="{} {}\n".format(pppnum,fitflag)
            OUTFLAGS.write(lineds9)


            if (countgal >= parvar.Ds9OutNum):
                countimg+=1
                OUTDS9.close()

                st = os.stat(ds9images)
                os.chmod(ds9images, st.st_mode | stat.S_IEXEC)


                ds9images=parvar.Ds9OutName+"-"+str(countimg)

                OUTDS9   = open(ds9images,"w")


                lineds9 = "#!/bin/bash \n"
                OUTDS9.write(lineds9)

                lineds9 = "ds9 -tile -cmap grey -cmap value {} {} -invert -log -zmax ".format(parvar.Contrast,parvar.Bias)
                OUTDS9.write(lineds9)

                countgal=1

            ermagtotal = 0
            erbt = 0
            ermagsersim = 0
            erresersim = 0
            ernsersim = 0
            ermagexpsim = 0
            errsexpsim = 0


            maskind = obj.Num == pppnum

            ind=np.where(maskind == True)

            indx=ind[0][0]



            obj.ra[indx] = ra
            obj.dec[indx] = dec
            obj.InputImage[indx] = inputimage
            obj.OutIm[indx] = outim
            obj.Objx[indx] = objx
            obj.PPPnum[indx] = pppnum
            obj.Restart[indx] = restart
            obj.Cluster[indx] =  cluster


            posxser = posxser  +  obj.OFFX[indx]
            posyser = posyser  +  obj.OFFY[indx]

#            posxexp = posxexp  +  obj.OFFX[maskind][0]
#            posyexp = posyexp  +  obj.OFFY[maskind][0]

            posxexp = posxexp  +  obj.OFFX[indx]
            posyexp = posyexp  +  obj.OFFY[indx]



            obj.PosXSer[indx] = posxser
            if np.isnan(erposxser):
                obj.ErPosXSer[indx]=99999
            else:
                obj.ErPosXSer[indx] = float(erposxser)

            obj.PosYSer[indx] =  posyser
            if np.isnan(erposyser):
                obj.ErPosYSer[indx]=99999
            else:
                obj.ErPosYSer[indx] = float(erposyser)

            obj.MagSer[indx] =  magser
#            obj.ErMagSer[indx] =  ermagser
            if np.isnan(ermagser):
                obj.ErMagSer[indx]=99999
            else:
                obj.ErMagSer[indx] = ermagser


            obj.ReSer[indx] =  reser
#            obj.ErReSer[indx] = erreser
            if np.isnan(erreser):
                obj.ErReSer[indx]=99999
            else:
                obj.ErReSer[indx] = erreser



            obj.NSer[indx] = nser
#            obj.ErNSer[indx] = ernser
            if np.isnan(ernser):
                obj.ErNSer[indx]=99999
            else:
                obj.ErNSer[indx] = ernser

            obj.AxisSer[indx] = axiser
#            obj.ErAxisSer[indx] = eraxiser
            if np.isnan(eraxiser):
                obj.ErAxisSer[indx]=99999
            else:
                obj.ErAxisSer[indx] = eraxiser

            obj.PaSer[indx] = paser
#            obj.ErPaSer[indx] = erpaser
            if np.isnan(erpaser):
                obj.ErPaSer[indx]=99999
            else:
                obj.ErPaSer[indx] = erpaser


            obj.KSer[indx] = kser

            obj.MeanMeSer[indx] = meanmeser
            obj.ErMeanMeSer[indx] = ermeanmeser

            obj.MeSer[indx] =  meser
            obj.ErMeSer[indx] =  ermeser



            obj.PosXExp[indx] = posxexp
#            obj.ErPosXExp[indx] = erposxexp
            if np.isnan(erposxexp):
                obj.ErPosXExp[indx]=99999
            else:
                obj.ErPosXExp[indx] = float(erposxexp)



            obj.PosYExp[indx] = posyexp
#            obj.ErPosYExp[indx] = erposyexp
            if np.isnan(erposyexp):
                obj.ErPosYExp[indx]=99999
            else:
                obj.ErPosYExp[indx] = float(erposyexp)



            obj.MagExp[indx] = magexp
#            obj.ErMagExp[indx] = ermagexp
            if np.isnan(ermagexp):
                obj.ErMagExp[indx]=99999
            else:
                obj.ErMagExp[indx] = float(ermagexp)



            obj.RsExp[indx] = rsexp
#            obj.ErRsExp[indx] = errsexp
            if np.isnan(errsexp):
                obj.ErRsExp[indx]=99999
            else:
                obj.ErRsExp[indx] = float(errsexp)

            obj.AxisExp[indx] = axisexp
#            obj.ErAxisExp[indx] = eraxisexp
            if np.isnan(eraxisexp):
                obj.ErAxisExp[indx]=99999
            else:
                obj.ErAxisExp[indx] = float(eraxisexp)


            obj.PaExp[indx] = paexp
#            obj.ErPaExp[indx] = erpaexp
            if np.isnan(erpaexp):
                obj.ErPaExp[indx]=99999
            else:
                obj.ErPaExp[indx] = float(erpaexp)


            obj.MeanMeExp[indx] = meanmeexp
            obj.ErMeanMeExp[indx] = ermeanmeexp

            obj.MeExp[indx] = meexp
            obj.MsExp[indx] = msexp

            obj.ErMeExp[indx] = ermeexp
            obj.ErMsExp[indx] = ermsexp



            obj.Sky[indx] =  sky
#            obj.ErSky[indx] = ersky
            if np.isnan(ersky):
                obj.ErSky[indx]=99999
            else:
                obj.ErSky[indx] = float(ersky)


            obj.MagTotal[indx] = magtotal
            obj.ErMagTotal[indx] = ermagtotal

            obj.BulgeTotal[indx] = bulgetotal
            obj.ErBt[indx] = erbt
            obj.ChiNu[indx] = chinu
            obj.Tidal[indx] = tidal
            obj.ObjChiNu[indx] = objchinu
            obj.Bump[indx] = bump
#            obj.MeanMesky[indx] = meanmesky
            obj.SNR[indx] = snr
            obj.NDof[indx] = ndof
            obj.FitFlag[indx]= fitflag

            obj.ReFit[indx] = refit



            if (parvar.FitFunc == "BD"):
                linefits = "{} {} {} {} {} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} \n".format(inputimage,outim,objx,pppnum,restart,cluster,ra,dec,posxser,erposxser,posyser,erposyser,magser,ermagser,reser,erreser,nser,ernser,axiser,eraxiser,paser,erpaser,kser,meanmeser,ermeanmeser,meser,ermeser,posxexp,erposxexp,posyexp,erposyexp,magexp,ermagexp,rsexp,errsexp,axisexp,eraxisexp,paexp,erpaexp,meanmeexp,ermeanmeexp,meexp,ermeexp,msexp,ermsexp,sky,ersky,magtotal,ermagtotal,bulgetotal,erbt,chinu,tidal,objchinu,bump,snr,ndof,fitflag,refit)
                OUTBT.write(linefits)
                    # 56 col.
            else:
                linefits = "{} {} {} {} {} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {} {} {} \n".format(inputimage,outim,objx,pppnum,restart,cluster,ra,dec,posxser,erposxser,posyser,erposyser,magser,ermagser,reser,erreser,nser,ernser,axiser,eraxiser,paser,erpaser,kser,meanmeser,ermeanmeser,meser,ermeser,sky,ersky,chinu,tidal,objchinu,bump,snr,ndof,fitflag,refit)
                OUTBT.write(linefits)
                   # 35 col.
#                    print ("Skipping object {}: {} <= 0, {} <=0 or {} = 0 \n".format(pppnum,reser,rsexp,sky))


            index+=1


        OUTDS9.close()
        OUTINFO.close()
        OUTFLAGS.close()
        OUTBT.close()


        st = os.stat(ds9images)
        os.chmod(ds9images, st.st_mode | stat.S_IEXEC)


        MakeRegs(parvar.FitFunc)


    return True


def MakeRegs(fitfunc):
    "Make DS9 Regions files for each object "

    INFO = open("inforegion")

    lines = (line.rstrip() for line in INFO) # All lines including the blank ones
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    lines = (line.rstrip() for line in lines)           # remove lines containing only comments
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0
    while index < len(lines):

        line = lines[index]
        (fileppp,pppnum,xser,yser,dx,dy,magtotal,pixreser,nser,bulgetotal,chinu,snr,axiser,axisexp,pixrsexp,tidal,objchinu,bump,paser,paexp,flag)=line.split()

        xup= float(xser)
        yup= float(yser) + 50

        xbot=float(xser)
        ybot=float(yser) - 30


        OUTREG = open (fileppp,"w")


        if (fitfunc == "BD"):

            linereg='point({},{}) # point=boxcircle color=blue font="times 14 bold" text={{ID: {} Flag: {} SNR: {} dx: {} dy: {} MAG: {} B/T: {} }}\n'.format(xser,yser,pppnum,flag,snr,dx,dy,magtotal,bulgetotal)
            OUTREG.write(linereg)


            linereg='# text({},{}) color=blue font="times 14 bold" text={{Tidal={} LocalChinu={} Bump={} chinu: {}  }} \n'.format(xup,yup,tidal,objchinu,bump,chinu)
            OUTREG.write(linereg)


            linereg= '# text({},{}) color=blue font="times 14 bold" text={{N :{} Re: {}  Rs: {}  qB: {}  qD: {} PaB: {} PaD: {} }} \n'.format(xbot,ybot,nser,pixreser,pixrsexp,axiser,axisexp,paser,paexp)
            OUTREG.write(linereg)

        else:

            linereg='point({},{}) # point=boxcircle color=blue font="times 14 bold" text={{ID: {} Flag: {} SNR: {} MAG: {}  }}\n'.format(xser,yser,pppnum,flag,snr,magtotal)
            OUTREG.write(linereg)


            linereg='# text({},{}) color=blue font="times 14 bold" text={{Tidal={} LocalChinu={} Bump={} chinu: {}  }} \n'.format(xup,yup,tidal,objchinu,bump,chinu)
            OUTREG.write(linereg)


            linereg= '# text({},{}) color=blue font="times 14 bold" text={{N :{} Re: {}  q: {}  PaB: {} }} \n'.format(xbot,ybot,nser,pixreser,axiser,paser)
            OUTREG.write(linereg)

        OUTREG.close()

        index+=1

    INFO.close()

    return True


def ReadConstraints(file, numobj):
    "Read constraints file and get constraints for parameters"

# Read Constraints File
# and get contraints for
# parameters

    flag = True

# Constraints file

    magmin = -99
    magmax = 99

    remin = -10000
    remax = 10000

    nmin = -100
    nmax = 100

    with open(file) as INCONS:
        try:

            # All lines including the blank ones
            lines = (line.rstrip() for line in INCONS)
            lines = (line.split('#', 1)[0]
                     for line in lines)  # remove comments
            # remove lines containing only comments
            lines = (line.rstrip() for line in lines)
            lines = (line for line in lines if line)  # Non-blank lines

            lines = list(lines)
            index = 0
            while index < len(lines):

                line = lines[index]
                (params) = line.split()

                if len(params) > 3:

                    if params[3] == "to":     # input image
                        if params[0] == numobj:

                            if(params[1] == "mag" or params[1] == 3):

                                magmin = params[2]
                                magmax = params[4]

                            if(params[1] == "re" or params[1] == 4 or params[1] == "rs"):

                                remin = params[2]
                                remax = params[4]

                            if(params[1] == "n" or params[1] == 5):

                                nmin = params[2]
                                nmax = params[4]

                index += 1

        except:
            print("Can't open constraints file or it was not provided \n")

    return (magmin, magmax, remin, remax, nmin, nmax)


def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

## solve the Sersic equation
# to get the dependence of K over
# Sersic index

    count = 1

    #limits
    lima=0
    limb=100

#fx is the function to solve
    fxa = fx(n,lima)
    fxb = fx(n,limb)

    resk= (lima + limb)/2

    fxres=fx(n,resk)


    if(fxa * fxb < 0):

        while(np.abs(fxres) > 0.00000001):

            if(fxa * fxres > 0):
                lima=resk
            elif(fxa * fxres < 0):
                limb=resk
            elif(fxres==0):
                break
            resk= (lima + limb)/2
            fxres=fx(n,resk)

            count+=1

            if (count >= 10000):
                break

    else:
        print("no solution in the range: ({},{})\n".format(lima,limb))

    return (resk)


def fx(n,k):
    "function to solve to get the relation between Sersic index and K"


    func = np.exp(scipy.special.gammaln(2*n)) - 2 * np.exp(scipy.special.gammaln(2*n)) * scipy.special.gammainc(2*n,k)


    return(func)


def BulgeTotal(mb,md):
    "Return Bulge to Total luminosity ratio"
# Given bulge mag and disk mag
# returns bulge to total luminosity
# ratio


    fluxb = 10 **((-mb)/2.5)

    fluxd = 10 **((-md)/2.5)

    BT= fluxb /(fluxb + fluxd)

    return (BT)


def Tidal(imgout,imsig,immask,num,xser,yser,xlo,xhi,ylo,yhi,goffx,goffy,rmin,sky,btflag):
    "Computes Tidal  values as defined in Tal et al. 2009 AJ. It algo computes Bumpiness"
    "(Blakeslee 2006 ApJ) value defined between rmin and rkron (see manual)"

# imgout = galfit output image
# imsig = sigma image (The one created by GALFIT)
# immask = Mask image: The original
# num = object number
# xser, yser = Sextractor object position on the big image
# xlo, xhi, ylo, yhi = coordinate of object of the big image
# goffx, goffy = coordinate tranformation from big to galfit image
# rmin = minimum radius to compute Tidal and Bumpiness (to avoid PSF Mismatch)
# sky = Background value
# btflag = is bulge/disk or sersic component?

    pixcount=0
    pixcountsnr=0
    pixcountchi=0
    sumtidal=0
    sumsig=0
    snr=0
    sigma=1
    flux=0
    sumflux=0
    meanflux=0
    var=0
    varb=0
    res=0
    resbump=0
    sumres=0
    sumchinu=0
    chinu=0
    varchi=1
    sflux=0
    sumsflux=0
    meanres=0
    sumvarres=0
    varres=0
    numbump=0
    ndof=0
    xglo=yglo=xghi=yghi=0
    xgser=ygser=0

    tidal=0
    objchinu=0
    bump=0
    snr=0


    hdug = fits.open(imgout) # galaxy
    datg = hdug[1].data

    datm = hdug[2].data

    hdus = fits.open(imsig) #sigma
    dats = hdus[0].data

    hdu = fits.open(immask) #mask
    dat = hdu[0].data


############################################################################
################ Code to compute rmin #######################################

#    for objchinu, Tidal and SNR
    maskm = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates

#   mask including rmin for Bumpiness only
    maskbum = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates

#   transforming coordinates to galfit output image
    xglo = xlo - goffx
    yglo = ylo - goffy

    xghi = xhi - goffx
    yghi = yhi - goffy

    xgser = xser - goffx
    ygser = yser - goffy

#############

    theta = 0

    ypos, xpos = np.mgrid[yglo - 1:yghi, xglo - 1:xghi]

    dx = xpos - xgser
    dy = ypos - ygser

    dist = np.sqrt(dx**2 + dy**2)

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / rmin, np.cos(landa) / rmin)

    xell = xgser + rmin * np.cos(angle)
    yell = ygser + rmin * np.sin(angle)

    dell = np.sqrt((xell - xgser)**2 + (yell - ygser)**2)

    mask = dist < dell

#  correcting for rmin
    maskbum[mask] = False

    if maskbum.any():

# for Bumpiness

        galfluxbum  = datg[yglo - 1:yghi, xglo - 1:xghi][maskbum]
        modfluxbum  = datm[yglo - 1:yghi, xglo - 1:xghi][maskbum]
        sigfluxbum  = dats[yglo - 1:yghi, xglo - 1:xghi][maskbum]
####

# for Tidal, SNR, and objchinu

        sumflux  = np.sum(datg[yglo - 1:yghi, xglo - 1:xghi][maskm] - sky)
        sumsig   = np.sum(dats[yglo - 1:yghi, xglo - 1:xghi][maskm])

        galflux  = datg[yglo - 1:yghi, xglo - 1:xghi][maskm]
        modflux  = datm[yglo - 1:yghi, xglo - 1:xghi][maskm]
        sigflux  = dats[yglo - 1:yghi, xglo - 1:xghi][maskm]

        resflux = (galflux - modflux)**2

#  local chinu

        varchi = sigflux**2
        chinu  = np.sum(resflux/varchi)

        pixcountchi = np.size(datm[yglo - 1:yghi, xglo - 1:xghi][maskm])

#  free parameters: sersic = 7 ;   Sersic + disk = 11

        if(pixcountchi > 11):

            if (btflag==False):
                ndof=pixcountchi-7

            elif(btflag==True):
                ndof=pixcountchi-11

            objchinu= chinu / ndof
        else:
            objchinu=-1
            ndof=-1


        if(np.size(dats[yglo - 1:yghi, xglo - 1:xghi][maskm]) > 0):
# snr
            meanflux = sumflux  / np.size(datg[yglo - 1:yghi, xglo - 1:xghi][maskm])
            sigma = sumsig / np.size(dats[yglo - 1:yghi, xglo - 1:xghi][maskm])
            snr = meanflux / sigma

        else:
            snr=-1

# Tidal parameter

        tgal = np.abs((galflux)/(modflux) - 1)
        sumtidal=np.sum(tgal)
        pixcountid=np.size(datm[yglo - 1:yghi, xglo - 1:xghi][maskm])

        if pixcountid > 0:
            tidal = (sumtidal / pixcountid)
        else:
            tidal=-1

# bumpiness

        resbump  = (galfluxbum - modfluxbum)**2
        sflux    = np.sum(np.abs(modfluxbum - sky))
        varres   = sigfluxbum * sigfluxbum
        numbump  = np.sum(resbump - varres)

        pixcountbum=np.size(datm[yglo - 1:yghi, xglo - 1:xghi][maskbum])


# Bumpiness
        if pixcountbum > 0:

            meansflux = sflux / pixcountbum
            meanres   = numbump / pixcountbum

            if (meanres < 0):
                meanres=0

            bump      = (np.sqrt(meanres)) / meansflux

        else:
            bump=-1



    hdug.close()
    hdus.close()
    hdu.close()


    return (tidal,objchinu,bump,snr,ndof)




def PosCor(offsetfile, btfile):
    "Correct x,y position from subpanel to the x,y of the actual image"

#    my ($flag,$line1,$line2);

    flag=True

    # copy file btfile
    runcmd = "cp {} btemp".format(btfile)
    errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

    IN1 = open("btemp")
    IN2 = open(offsetfile)

    lines = (line.rstrip() for line in IN1) # All lines including the blank ones
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    lines = (line.rstrip() for line in lines)     # remove lines containing only comments
    lines = (line for line in lines if line)  # Non-blank lines


    lines2 = (line2.rstrip() for line2 in IN2) # All lines including the blank ones
    lines2 = (line2.split('#', 1)[0] for line2 in lines2)  # remove comments
    lines2 = (line2.rstrip() for line2 in lines2)     # remove lines containing only comments
    lines2 = (line2 for line2 in lines2 if line2)  # Non-blank lines


    lines = list(lines)
    index = 0

    lines2 = list(lines2)
    index2 = 0

    POSOUT = open (btfile,w)

    while index < len(lines):
        line = lines[index]
        (inputimage,outim,objt,pppnum,restart,cluster,ra,dec,posxser,erposxser,posyser,erposyser,magser,ermagser,reser,erreser,nser,ernser,axisser,eraxisser,paser,erpaser,kser,meanmeser,meser,posxexp,erposxexp,posyexp,erposyexp,magexp,ermagexp,rsexp,errsexp,axisexp,eraxisexp,paexp,erpaexp,meanmeexp,meexp,sky,ersky,magtotal,ermagtotal,bulgetotal,erbt,chinu,tidal,objchinu,bump,meanmesky,snr,ndof,fitflag)=line.split()

        while index2 < len(lines2):
            line2 = lines2[index2]
            (num,output,xoff,yoff,intpos,xini,yini)=line2.split()

            if (pppnum == num ):
                posxser = posxser + xoff
                posyser = posyser + yoff

                lineout= "{1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22} {23} {24} {25} {26} {27} {28} {29} {30} {31} {32} {33} {34} {35} {36} {37} {38} {39} {40} {41} {42} {43} {44} {45} {46} {47} {49} {50} {51} {52} {53}".format(inputimage, outim, objt, pppnum, restart, cluster, ra, dec, posxser, erposxser, posyser, erposyser, magser, ermagser, reser, erreser, nser, ernser, axisser, eraxisser, paser, erpaser, kser, meanmeser, meser, posxexp, erposxexp, posyexp, erposyexp, magexp, ermagexp, rsexp, errsexp, axisexp, eraxisexp, paexp, erpaexp, meanmeexp, meexp, sky, ersky, magtotal, ermagtotal, bulgetotal, erbt, chinu, tidal, objchinu, bump, meanmesky, snr, ndof, fitflag)
                POSOUT.write(lineout)

                break

            index2+=1
        index+=1
        index2=0


    os.remove("btemp")

    IN.close()
    IN2.close()
    POSOUT.close()

    return (flag)


def JoinSexOut(parvar, obj):
    "Join Sextractor with output catalog and add a flag indicating if the fitting"
    "crashed or it was successfully fitted. fitflag=8"


    flag=True

    OBJIN  = open (parvar.ListObjs)
#    BTIN   = open (parvar.BtFile)
    SEXOUT = open (parvar.SexOut,"w")

    lines = (line.rstrip() for line in OBJIN) # All lines including the blank ones
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    lines = (line.rstrip() for line in lines)     # remove lines containing only comments
    lines = (line for line in lines if line)  # Non-blank lines


#    lines2 = (line2.rstrip() for line2 in BTIN) # All lines including the blank ones
#    lines2 = (line2.split('#', 1)[0] for line2 in lines2)  # remove comments
#    lines2 = (line2.rstrip() for line2 in lines2)     # remove lines containing only comments
#    lines2 = (line2 for line2 in lines2 if line2)  # Non-blank lines

    lines = list(lines)
    index = 0

#    lines2 = list(lines2)
#    index2 = 0


    while index < len(lines):

        line = lines[index]
        btflag = False
        (idx,i)=line.split()
        idx=int(idx)
        i=int(i)

############
        maskind = obj.Num == idx

        ind=np.where(maskind == True)

        indx=ind[0][0]
############

        fitflag=core.SearchObj(idx)

        if fitflag == False:
            obj.FitFlag[i] = obj.FitFlag[i] + 8
            obj.ReFit[i]=True

            print("obj.Num",obj.Num[i])
            print("obj.FitFlag ",obj.FitFlag[i])
            print("obj.ReFit ",obj.ReFit[i])


        if maskind.any():


            if (parvar.FitFunc == "BD"):

# 84 cols

                lineout= "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10:.3f} {11:.3f} {12} {13} {14} {15} {16} {17} {18} {19} {20:.3f} {21:.3f} {22} {23} {24} {25} {26} {27} {28} {29} {30} {31} {32} {33:.3f} {34:.3f} {35:.3f} {36:.3f} {37:.3f} {38:.3f} {39:.3f} {40:.3f} {41:.3f} {42:.3f} {43:.3f} {44:.3f} {45:.3f} {46:.3f} {47:.3f} {48:.3f} {49:.3f} {50:.3f} {51:.3f} {52:.3f} {53:.3f} {54:.3f} {55:.3f} {56:.3f} {57:.3f} {58:.3f} {59:.3f} {60:.3f} {61:.3f} {62:.3f} {63:.3f} {64:.3f} {65:.3f} {66:.3f} {67:.3f} {68:.3f} {69:.3f} {70:.3f} {71:.3f} {72:.3f} {73:.3f} {74:.3f} {75:.3f} {76:.3f} {77:.3f} {78:.3f} {79:.3f} {80:.3f} {81:.3f} {82} {83} {84} \n".format(obj.Num[indx], obj.RA[indx], obj.Dec[indx], obj.XPos[indx], obj.YPos[indx], obj.Mag[indx], obj.Kron[indx], obj.FluxRad[indx], obj.IsoArea[indx], obj.AIm[indx], obj.AR[indx], obj.Angle[indx], obj.Background[indx], obj.Class[indx], obj.Flag[indx], obj.XMin[indx], obj.XMax[indx], obj.YMin[indx], obj.YMax[indx], obj.Sersic[indx], obj.RSky[indx], obj.RKron[indx], obj.SkyFlag[indx], obj.OFFX[indx], obj.OFFY[indx], obj.InputImage[indx], obj.OutIm[indx], obj.Objx[indx], obj.PPPnum[indx], obj.Restart[indx], obj.Cluster[indx], obj.ra[indx], obj.dec[indx], obj.PosXSer[indx], obj.ErPosXSer[indx], obj.PosYSer[indx], obj.ErPosYSer[indx], obj.MagSer[indx], obj.ErMagSer[indx], obj.ReSer[indx], obj.ErReSer[indx], obj.NSer[indx], obj.ErNSer[indx], obj.AxisSer[indx], obj.ErAxisSer[indx], obj.PaSer[indx], obj.ErPaSer[indx], obj.KSer[indx], obj.MeanMeSer[indx], obj.ErMeanMeSer[indx], obj.MeSer[indx], obj.ErMeSer[indx], obj.PosXExp[indx], obj.ErPosXExp[indx], obj.PosYExp[indx], obj.ErPosYExp[indx], obj.MagExp[indx], obj.ErMagExp[indx], obj.RsExp[indx], obj.ErRsExp[indx], obj.AxisExp[indx], obj.ErAxisExp[indx], obj.PaExp[indx], obj.ErPaExp[indx], obj.MeanMeExp[indx], obj.ErMeanMeExp[indx], obj.MeExp[indx], obj.ErMeExp[indx], obj.MsExp[indx], obj.ErMsExp[indx], obj.Sky[indx], obj.ErSky[indx], obj.MagTotal[indx], obj.ErMagTotal[indx], obj.BulgeTotal[indx], obj.ErBt[indx], obj.ChiNu[indx], obj.Tidal[indx], obj.ObjChiNu[indx], obj.Bump[indx], obj.SNR[indx], obj.NDof[indx], obj.FitFlag[indx], obj.ReFit[indx])
                SEXOUT.write(lineout)

            else:
# 61 cols
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    7
                lineout= "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10:.3f} {11:.3f} {12} {13} {14} {15} {16} {17} {18} {19} {20:.3f} {21:.3f} {22} {23} {24} {25} {26} {27} {28} {29} {30} {31} {32} {33:.3f} {34:.3f} {35:.3f} {36:.3f} {37:.3f} {38:.3f} {39:.3f} {40:.3f} {41:.3f} {42:.3f} {43:.3f} {44:.3f} {45:.3f} {46:.3f} {47:.3f} {48:.3f} {49:.3f} {50:.3f} {51:.3f} {52:.3f} {53:.3f} {54:.3f} {55:.3f} {56:.3f} {57:.3f} {58:.3f} {59} {60} {61} \n".format(obj.Num[indx], obj.RA[indx], obj.Dec[indx], obj.XPos[indx], obj.YPos[indx], obj.Mag[indx], obj.Kron[indx], obj.FluxRad[indx], obj.IsoArea[indx], obj.AIm[indx], obj.AR[indx], obj.Angle[indx], obj.Background[indx], obj.Class[indx], obj.Flag[indx], obj.XMin[indx], obj.XMax[indx], obj.YMin[indx], obj.YMax[indx], obj.Sersic[indx], obj.RSky[indx], obj.RKron[indx], obj.SkyFlag[indx], obj.OFFX[indx], obj.OFFY[indx], obj.InputImage[indx], obj.OutIm[indx], obj.Objx[indx], obj.PPPnum[indx], obj.Restart[indx], obj.Cluster[indx], obj.ra[indx], obj.dec[indx], obj.PosXSer[indx], obj.ErPosXSer[indx], obj.PosYSer[indx], obj.ErPosYSer[indx], obj.MagSer[indx], obj.ErMagSer[indx], obj.ReSer[indx], obj.ErReSer[indx], obj.NSer[indx], obj.ErNSer[indx], obj.AxisSer[indx], obj.ErAxisSer[indx], obj.PaSer[indx], obj.ErPaSer[indx], obj.KSer[indx], obj.MeanMeSer[indx], obj.ErMeanMeSer[indx], obj.MeSer[indx], obj.ErMeSer[indx], obj.Sky[indx], obj.ErSky[indx], obj.ChiNu[indx], obj.Tidal[indx], obj.ObjChiNu[indx], obj.Bump[indx], obj.SNR[indx], obj.NDof[indx], obj.FitFlag[indx], obj.ReFit[indx] )



                SEXOUT.write(lineout)

        index+=1

    OBJIN.close()
#    BTIN.close()
    SEXOUT.close()

    return (flag)




def SelectColumns(infile, paramfile, outfile,headflag):
    "Selects columns from dgcg file to print in outfile"

    (name,ext,dgtrash)=infile.split(".")

    if ext == "ser":
        fitfunc = "Sersic"

    elif ext == "bd":
        fitfunc = "BD"
    else:
        print("SelectColumns: file extention not recognized\n")
        sys.exit()


    SexNumFlag = False
    SexRAFlag = False
    SexDecFlag = False
    SexXPosFlag = False
    SexYPosFlag = False
    SexMagFlag = False
    SexKronFlag = False
    SexFluxRadFlag = False
    SexIsoAreaFlag = False
    SexAImFlag = False
    SexARFlag = False
    SexAngleFlag = False
    SexBackgroundFlag = False
    SexClassFlag = False
    SexFlagFlag = False
    XMinFlag = False
    XMaxFlag = False
    YMinFlag = False
    YMaxFlag = False
    InSersicFlag = False
    RSkyFlag = False
    RKronFlag = False
    SkyFlagFlag = False
    OFFXFlag = False
    OFFYFlag = False
    InputImageFlag = False
    OutImageFlag = False
    InputFileFlag = False
    NumFlag = False
    RestartFlag = False
    ClusterFlag = False
    RAFlag = False
    DECFlag = False

    PosXSerFlag = False
    ErPosXSerFlag = False
    PosYSerFlag = False
    ErPosYSerFlag = False
    MagSerFlag = False
    ErMagSerFlag = False
    ReSerFlag = False
    ErReSerFlag = False
    NSerFlag = False
    ErNSerFlag = False
    AxisSerFlag = False
    ErAxisSerFlag = False
    PaSerFlag = False
    ErPaSerFlag = False
    KSerFlag = False
    MeanMeSerFlag = False
    MeSerFlag = False
    ErMeanMeSerFlag = False
    ErMeSerFlag = False



    PosXExpFlag = False
    ErPosXExpFlag = False
    PosYExpFlag = False
    ErPosYExpFlag = False
    MagExpFlag = False
    ErMagExpFlag = False
    RsExpFlag = False
    ErRsExpFlag = False
    AxisExpFlag = False
    ErAxisExpFlag = False
    PaExpFlag = False
    ErPaExpFlag = False
    MeanMeExpFlag = False
    MeExpFlag = False
    ErMeanMeExpFlag = False
    ErMeExpFlag = False
    MsExpFlag = False
    ErMsExpFlag = False



    SkyFlag = False
    ErSkyFlag = False
    MagTotalFlag = False
    ErMagTotalFlag = False
    BulgeTotalFlag = False
    ErBtFlag = False
    ChiNuFlag = False
    TidalFlag = False
    ObjChiNuFlag = False
    BumpFlag = False
    SNRFlag = False
    NDofFlag = False
    FitFlagFlag = False
    ReFitFlag = False



    IN  = open (infile)
    PARIN  = open (paramfile)
    OUT = open (outfile,"w")

    lines = (line.rstrip() for line in PARIN) # All lines including the blank ones
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    lines = (line.rstrip() for line in lines)     # remove lines containing only comments
    lines = (line for line in lines if line)  # Non-blank lines


    for line in lines:

        (params)=line.split()
        param=params[0]
        FoundFlag = False


        if param == "SexNum":
            SexNumFlag=True
            FoundFlag = True
        if param == "SexRA":
            SexRAFlag=True
            FoundFlag = True
        if param == "SexDec":
            SexDecFlag=True
            FoundFlag = True
        if param == "SexXPos":
            SexXPosFlag=True
            FoundFlag = True
        if param == "SexYPos":
            SexYPosFlag=True
            FoundFlag = True
        if param == "SexMag":
            SexMagFlag=True
            FoundFlag = True
        if param == "SexKron":
            SexKronFlag=True
            FoundFlag = True
        if param == "SexFluxRad":
            SexFluxRadFlag=True
            FoundFlag = True
        if param == "SexIsoArea":
            SexIsoAreaFlag=True
            FoundFlag = True
        if param == "SexAIm":
            SexAImFlag=True
            FoundFlag = True
        if param == "SexAR":
            SexARFlag=True
            FoundFlag = True
        if param == "SexAngle":
            SexAngleFlag=True
            FoundFlag = True
        if param == "SexBackground":
            SexBackgroundFlag=True
            FoundFlag = True
        if param == "SexClass":
            SexClassFlag=True
            FoundFlag = True
        if param == "SexFlag":
            SexFlagFlag=True
            FoundFlag = True
        if param == "XMin":
            XMinFlag=True
            FoundFlag = True
        if param == "XMax":
            XMaxFlag=True
            FoundFlag = True
        if param == "YMin":
            YMinFlag=True
            FoundFlag = True
        if param == "YMax":
            YMaxFlag=True
            FoundFlag = True
        if param == "InSersic":
            InSersicFlag=True
            FoundFlag = True
        if param == "RSky":
            RSkyFlag=True
            FoundFlag = True
        if param == "RKron":
            RKronFlag=True
            FoundFlag = True
        if param == "SkyFlag":
            SkyFlagFlag=True
            FoundFlag = True
        if param == "OFFX":
            OFFXFlag=True
            FoundFlag = True
        if param == "OFFY":
            OFFYFlag=True
            FoundFlag = True
        if param == "InputImage":
            InputImageFlag=True
            FoundFlag = True
        if param == "OutImage":
            OutImageFlag=True
            FoundFlag = True
        if param == "InputFile":
            InputFileFlag=True
            FoundFlag = True
        if param == "Num":
            NumFlag=True
            FoundFlag = True
        if param == "Restart":
            RestartFlag=True
            FoundFlag = True
        if param == "Cluster":
            ClusterFlag=True
            FoundFlag = True
        if param == "RA":
            RAFlag=True
            FoundFlag = True
        if param == "DEC":
            DECFlag=True
            FoundFlag = True
        if param == "PosXSer":
            PosXSerFlag=True
            FoundFlag = True
        if param == "ErPosXSer":
            ErPosXSerFlag=True
            FoundFlag = True
        if param == "PosYSer":
            PosYSerFlag=True
            FoundFlag = True
        if param == "ErPosYSer":
            ErPosYSerFlag=True
            FoundFlag = True
        if param == "MagSer":
            MagSerFlag=True
            FoundFlag = True
        if param == "ErMagSer":
            ErMagSerFlag=True
            FoundFlag = True
        if param == "ReSer":
            ReSerFlag=True
            FoundFlag = True
        if param == "ErReSer":
            ErReSerFlag=True
            FoundFlag = True
        if param == "NSer":
            NSerFlag=True
            FoundFlag = True
        if param == "ErNSer":
            ErNSerFlag=True
            FoundFlag = True
        if param == "AxisSer":
            AxisSerFlag=True
            FoundFlag = True
        if param == "ErAxisSer":
            ErAxisSerFlag=True
            FoundFlag = True
        if param == "PaSer":
            PaSerFlag=True
            FoundFlag = True
        if param == "ErPaSer":
            ErPaSerFlag=True
            FoundFlag = True
        if param == "KSer":
            KSerFlag=True
            FoundFlag = True
        if param == "MeanMeSer":
            MeanMeSerFlag=True
            FoundFlag = True
        if param == "MeSer":
            MeSerFlag=True
            FoundFlag = True
        if param == "ErMeanMeSer":
            ErMeanMeSerFlag=True
            FoundFlag = True
        if param == "ErMeSer":
            ErMeSerFlag=True
            FoundFlag = True
        if fitfunc == "BD":
            if param == "PosXExp":
                PosXExpFlag=True
                FoundFlag = True
            if param == "ErPosXExp":
                ErPosXExpFlag=True
                FoundFlag = True
            if param == "PosYExp":
                PosYExpFlag=True
                FoundFlag = True
            if param == "ErPosYExp":
                ErPosYExpFlag=True
                FoundFlag = True
            if param == "MagExp":
                MagExpFlag=True
                FoundFlag = True
            if param == "ErMagExp":
                ErMagExpFlag=True
                FoundFlag = True
            if param == "RsExp":
                RsExpFlag=True
                FoundFlag = True
            if param == "ErRsExp":
                ErRsExpFlag=True
                FoundFlag = True
            if param == "AxisExp":
                AxisExpFlag=True
                FoundFlag = True
            if param == "ErAxisExp":
                ErAxisExpFlag=True
                FoundFlag = True
            if param == "PaExp":
                PaExpFlag=True
                FoundFlag = True
            if param == "ErPaExp":
                ErPaExpFlag=True
                FoundFlag = True
            if param == "MeanMeExp":
                MeanMeExpFlag=True
                FoundFlag = True
            if param == "MeExp":
                MeExpFlag=True
                FoundFlag = True
            if param == "ErMeanMeExp":
                ErMeanMeExpFlag=True
                FoundFlag = True
            if param == "ErMeExp":
                ErMeExpFlag=True
                FoundFlag = True
            if param == "MsExp":
                MsExpFlag=True
                FoundFlag = True
            if param == "ErMsExp":
                ErMsExpFlag=True
                FoundFlag = True
            if param == "MagTotal":
                MagTotalFlag=True
                FoundFlag = True
            if param == "ErMagTotal":
                ErMagTotalFlag=True
                FoundFlag = True
            if param == "BulgeTotal":
                BulgeTotalFlag=True
                FoundFlag = True
            if param == "ErBt":
                ErBtFlag=True
                FoundFlag = True

        if param == "Sky":
            SkyGalFlag=True
            FoundFlag = True
        if param == "ErSky":
            ErSkyFlag=True
            FoundFlag = True
        if param == "ChiNu":
            ChiNuFlag=True
            FoundFlag = True
        if param == "Tidal":
            TidalFlag=True
            FoundFlag = True
        if param == "ObjChiNu":
            ObjChiNuFlag=True
            FoundFlag = True
        if param == "Bump":
            BumpFlag=True
            FoundFlag = True
        if param == "SNR":
            SNRFlag=True
            FoundFlag = True
        if param == "NDof":
            NDofFlag=True
            FoundFlag = True
        if param == "FitFlag":
            FitFlagFlag=True
            FoundFlag = True
        if param == "ReFit":
            ReFitFlag=True
            FoundFlag = True

        if FoundFlag == False:
            errmsg = "{} parameter not found".format(param)
            print (errmsg)


## Header

    head=""

    num=1


    if headflag == True:


        if SexNumFlag==True:
            head+="#   {:2d}  SexNum                Sextractor Running object number   \n".format(num)
            num+=1
        if SexRAFlag==True:
            head+="#   {:2d}  SexRA                 Sextractor Right Ascension \n".format(num)
            num+=1
        if SexDecFlag==True:
            head+="#   {:2d}  SexDec                Sextractor Declination \n".format(num)
            num+=1
        if SexXPosFlag==True:
            head+="#   {:2d}  SexXPos               Sextractor X position \n".format(num)
            num+=1
        if SexYPosFlag==True:
            head+="#   {:2d}  SexYPos               Sextractor Y position \n".format(num)
            num+=1
        if SexMagFlag==True:
            head+="#   {:2d}  SexMag                Sextractor Magnitude \n".format(num)
            num+=1
        if SexKronFlag==True:
            head+="#   {:2d}  SexKron               Sextractor Kron Radius\n".format(num)
            num+=1
        if SexFluxRadFlag==True:
            head+="#   {:2d}  SexFluxRad            Sextractor Flux Radius \n".format(num)
            num+=1
        if SexIsoAreaFlag==True:
            head+="#   {:2d}  SexIsoArea            Sextractor Isophotal Area  \n".format(num)
            num+=1
        if SexAImFlag==True:
            head+="#   {:2d}  SexAIm                Sextractor Major axis \n".format(num)
            num+=1
        if SexARFlag==True:
            head+="#   {:2d}  SexAR                 Sextractor Axis ratio \n".format(num)
            num+=1
        if SexAngleFlag==True:
            head+="#   {:2d}  SexAngle              Sextractor Angle respect x-axis \n".format(num)
            num+=1
        if SexBackgroundFlag==True:
            head+="#   {:2d}  SexBackground         Sextractor Background  \n".format(num)
            num+=1
        if SexClassFlag==True:
            head+="#   {:2d}  SexClass              Sextractor star/galaxy classifier \n".format(num)
            num+=1
        if SexFlagFlag==True:
            head+="#   {:2d}  SexFlag               Sextractor Extraction flags \n".format(num)
            num+=1
        if XMinFlag==True:
            head+="#   {:2d}  XMin                  Minimum x-coordinate among detected pixels\n".format(num)
            num+=1
        if XMaxFlag==True:
            head+="#   {:2d}  XMax                  Maximum x-coordinate among detected pixels\n".format(num)
            num+=1
        if YMinFlag==True:
            head+="#   {:2d}  YMin                  Minimum y-coordinate among detected pixels \n".format(num)
            num+=1
        if YMaxFlag==True:
            head+="#   {:2d}  YMax                  Maximum y-coordinate among detected pixels \n".format(num)
            num+=1
        if InSersicFlag==True:
            head+="#   {:2d}  InSersic              initial parameter for Sersic index \n".format(num)
            num+=1
        if RSkyFlag==True:
            head+="#   {:2d}  RSky                  Scale Radius object for sky computation \n".format(num)
            num+=1
        if RKronFlag==True:
            head+="#   {:2d}  RKron                 Kron Scale Radius object to determine simultaneous fitting \n".format(num)
            num+=1
        if SkyFlagFlag==True:
            head+="#   {:2d}  SkyFlag               Sky flag output \n".format(num)
            num+=1
        if OFFXFlag==True:
            head+="#   {:2d}  OFFX                  Offset x for tile image\n".format(num)
            num+=1
        if OFFYFlag==True:
            head+="#   {:2d}  OFFY                  Offset y for tile image\n".format(num)
            num+=1
        if InputImageFlag==True:
            head+="#   {:2d}  InputImage            Input image  \n".format(num)
            num+=1
        if OutImageFlag==True:
            head+="#   {:2d}  OutImage              output image with dgcg extention \n".format(num)
            num+=1
        if InputFileFlag==True:
            head+="#   {:2d}  InputFile             Input parameter file\n".format(num)
            num+=1
        if NumFlag==True:
            head+="#   {:2d}  Num                   Object Number \n".format(num)
            num+=1
        if RestartFlag==True:
            head+="#   {:2d}  Restart               Restart GALFIT File \n".format(num)
            num+=1
        if ClusterFlag==True:
            head+="#   {:2d}  Cluster               Name of the galaxy cluster \n".format(num)
            num+=1
        if RAFlag==True:
            head+="#   {:2d}  RA                    Right Ascension \n".format(num)
            num+=1
        if DECFlag==True:
            head+="#   {:2d}  DEC                   Declination \n".format(num)
            num+=1
        if PosXSerFlag==True:
            head+="#   {:2d}  PosXSer               Sersic X position \n".format(num)
            num+=1
        if ErPosXSerFlag==True:
            head+="#   {:2d}  ErPosXSer             Error Sersic X position\n".format(num)
            num+=1
        if PosYSerFlag==True:
            head+="#   {:2d}  PosYSer               Sersic Y position \n".format(num)
            num+=1
        if ErPosYSerFlag==True:
            head+="#   {:2d}  ErPosYser             Error Sersic Y position\n".format(num)
            num+=1
        if MagSerFlag==True:
            head+="#   {:2d}  MagSer                Sersic Magnitud\n".format(num)
            num+=1
        if ErMagSerFlag==True:
            head+="#   {:2d}  ErMagSer              Error Sersic Magnitud \n".format(num)
            num+=1
        if ReSerFlag==True:
            head+="#   {:2d}  ReSer                 Sersic effective radius \n".format(num)
            num+=1
        if ErReSerFlag==True:
            head+="#   {:2d}  ErReSer               Error Sersic effective radius\n".format(num)
            num+=1
        if NSerFlag==True:
            head+="#   {:2d}  NSer                  Sersic index\n".format(num)
            num+=1
        if ErNSerFlag==True:
            head+="#   {:2d}  ErNSer                Error Sersic index \n".format(num)
            num+=1
        if AxisSerFlag==True:
            head+="#   {:2d}  AxisSer               Sersic axis ratio\n".format(num)
            num+=1
        if ErAxisSerFlag==True:
            head+="#   {:2d}  ErAxisSer             Error Sersic axis ratio \n".format(num)
            num+=1
        if PaSerFlag==True:
            head+="#   {:2d}  PaSer                 Sersic position angular \n".format(num)
            num+=1
        if ErPaSerFlag==True:
            head+="#   {:2d}  ErPaSer               Error sersic position angular  \n".format(num)
            num+=1
        if KSerFlag==True:
            head+="#   {:2d}  KSer                  Sersic K parameter related to n \n".format(num)
            num+=1
        if MeanMeSerFlag==True:
            head+="#   {:2d}  MeanMeSer             Sersic Mean surface brightness at effective radius \n".format(num)
            num+=1
        if MeSerFlag==True:
            head+="#   {:2d}  MeSer                 Sersic surface brightness at effective radius \n".format(num)
            num+=1
        if ErMeanMeSerFlag==True:
            head+="#   {:2d}  ErMeanMeSer             Sersic Mean surface brightness at effective radius \n".format(num)
            num+=1
        if ErMeSerFlag==True:
            head+="#   {:2d}  ErMeSer                 Sersic surface brightness at effective radius \n".format(num)
            num+=1
        if PosXExpFlag==True:
            head+="#   {:2d}  PosXExp               Exponential X position\n".format(num)
            num+=1
        if ErPosXExpFlag==True:
            head+="#   {:2d}  ErPosXExp             Error Exponential X position\n".format(num)
            num+=1
        if PosYExpFlag==True:
            head+="#   {:2d}  PosYExp               Exponential Y position\n".format(num)
            num+=1
        if ErPosYExpFlag==True:
            head+="#   {:2d}  ErPosYExp             Error Exponential Y position \n".format(num)
            num+=1
        if MagExpFlag==True:
            head+="#   {:2d}  MagExp                 Exponential Magnitud \n".format(num)
            num+=1
        if ErMagExpFlag==True:
            head+="#   {:2d}  ErMagExp              Error Exponential  Magnitud \n".format(num)
            num+=1
        if RsExpFlag==True:
            head+="#   {:2d}  RsExp                 Exponential Scale Radius \n".format(num)
            num+=1
        if ErRsExpFlag==True:
            head+="#   {:2d}  ErRsExp               Error Exponential Scale Radius \n".format(num)
            num+=1
        if AxisExpFlag==True:
            head+="#   {:2d}  AxisExp               Exponential Axis Ratio \n".format(num)
            num+=1
        if ErAxisExpFlag==True:
            head+="#   {:2d}  ErAxisExp             Error Exponential Axis Ratio  \n".format(num)
            num+=1
        if PaExpFlag==True:
            head+="#   {:2d}  PaExp                 Exponential Postion angular \n".format(num)
            num+=1
        if ErPaExpFlag==True:
            head+="#   {:2d}  ErPaExp               Error Exponential angular \n".format(num)
            num+=1
        if MeanMeExpFlag==True:
            head+="#   {:2d}  MeanMeExp             Exponential mean surface brightness at effective radius \n".format(num)
            num+=1
        if MeExpFlag==True:
            head+="#   {:2d}  MeExp                 Exponential surface brightness at effective radius\n".format(num)
            num+=1
        if ErMeanMeExpFlag==True:
            head+="#   {:2d}  ErMeanMeExp             Exponential mean surface brightness at effective radius \n".format(num)
            num+=1
        if ErMeExpFlag==True:
            head+="#   {:2d}  ErMeExp                 Exponential surface brightness at effective radius\n".format(num)
            num+=1
        if MsExpFlag==True:
            head+="#   {:2d}  MsExp                 Exponential surface brightness at effective radius\n".format(num)
            num+=1
        if ErMsExpFlag==True:
            head+="#   {:2d}  ErMsExp                 Exponential surface brightness at effective radius\n".format(num)
            num+=1
        if SkyGalFlag==True:
            head+="#   {:2d}  Sky                   Background Sky value \n".format(num)
            num+=1
        if ErSkyFlag==True:
            head+="#   {:2d}  ErSky                 Error Background Sky value\n".format(num)
            num+=1
        if MagTotalFlag==True:
            head+="#   {:2d}  MagTotal              Total magnitude \n".format(num)
            num+=1
        if ErMagTotalFlag==True:
            head+="#   {:2d}  ErMagTotal            Error total magnitude \n".format(num)
            num+=1
        if BulgeTotalFlag==True:
            head+="#   {:2d}  BulgeTotal            Bulge to Total magnitude ratio \n".format(num)
            num+=1
        if ErBtFlag==True:
            head+="#   {:2d}  ErBt                  Error bulge to Total magnitude ratio \n".format(num)
            num+=1
        if ChiNuFlag==True:
            head+="#   {:2d}  ChiNu                 Chi square divided over degrees of freedom \n".format(num)
            num+=1
        if TidalFlag==True:
            head+="#   {:2d}  Tidal                 Tidal value (See manual) \n".format(num)
            num+=1
        if ObjChiNuFlag==True:
            head+="#   {:2d}  ObjChiNu              ChiNu computed inside Kron Radius \n".format(num)
            num+=1
        if BumpFlag==True:
            head+="#   {:2d}  Bump                  Bumpiness at Kron Radius \n".format(num)
            num+=1
        if SNRFlag==True:
            head+="#   {:2d}  SNR                   Signal to Noise Ratio computed inside Kron Radius \n".format(num)
            num+=1
        if NDofFlag==True:
            head+="#   {:2d}  NDof                  Number of deegres of freedom \n".format(num)
            num+=1
        if FitFlagFlag==True:
            head+="#   {:2d}  FitFlag               Output flags. See manual \n".format(num)
            num+=1
        if ReFitFlag==True:
            head+="#   {:2d}  ReFit                 Flag output. Object needs refit? \n".format(num)
            num+=1

        OUT.write(head)

###### select columns

    lines = (line.rstrip() for line in IN) # All lines including the blank ones
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    lines = (line.rstrip() for line in lines)     # remove lines containing only comments
    lines = (line for line in lines if line)  # Non-blank lines


    for line in lines:


        if (fitfunc == "BD"):
# 79 cols

            (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, PosXExp, ErPosXExp, PosYExp, ErPosYExp, MagExp, ErMagExp, RsExp, ErRsExp, AxisExp, ErAxisExp, PaExp, ErPaExp, MeanMeExp, ErMeanMeExp, MeExp, ErMeExp, MsExp, ErMsExp, Sky, ErSky, MagTotal, ErMagTotal, BulgeTotal, ErBt, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit)=line.split()
        else:
#64 cols

            (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, Sky, ErSky, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit)=line.split()


        outline=""

        if SexNumFlag==True:
            outline+="{} ".format(Num)

        if SexRAFlag==True:
            outline+="{} ".format(RA)

        if SexDecFlag==True:
            outline+="{} ".format(Dec)

        if SexXPosFlag==True:
            outline+="{} ".format(XPos)

        if SexYPosFlag==True:
            outline+="{} ".format(YPos)

        if SexMagFlag==True:
            outline+="{} ".format(Mag)

        if SexKronFlag==True:
            outline+="{} ".format(Kron)

        if SexFluxRadFlag==True:
            outline+="{} ".format(FluxRad)

        if SexIsoAreaFlag==True:
            outline+="{} ".format(IsoArea)

        if SexAImFlag==True:
            outline+="{} ".format(AIm)

        if SexARFlag==True:
            outline+="{} ".format(AR)

        if SexAngleFlag==True:
            outline+="{} ".format(Angle)

        if SexBackgroundFlag==True:
            outline+="{} ".format(Background)

        if SexClassFlag==True:
            outline+="{} ".format(Class)

        if SexFlagFlag==True:
            outline+="{} ".format(Flag)

        if XMinFlag==True:
            outline+="{} ".format(XMin)

        if XMaxFlag==True:
            outline+="{} ".format(XMax)

        if YMinFlag==True:
            outline+="{} ".format(YMin)

        if YMaxFlag==True:
            outline+="{} ".format(YMax)

        if InSersicFlag==True:
            outline+="{} ".format(Sersic)

        if RSkyFlag==True:
            outline+="{} ".format(RSky)

        if RKronFlag==True:
            outline+="{} ".format(RKron)

        if SkyFlagFlag==True:
            outline+="{} ".format(SkyFlag)

        if OFFXFlag==True:
            outline+="{} ".format(OFFX)

        if OFFYFlag==True:
            outline+="{} ".format(OFFY)

        if InputImageFlag==True:
            outline+="{} ".format(InputImage)

        if OutImageFlag==True:
            outline+="{} ".format(OutIm)

        if InputFileFlag==True:
            outline+="{} ".format(Objx)

        if NumFlag==True:
            outline+="{} ".format(PPPnum)

        if RestartFlag==True:
            outline+="{} ".format(Restart)

        if ClusterFlag==True:
            outline+="{} ".format(Cluster)

        if RAFlag==True:
            outline+="{} ".format(ra)

        if DECFlag==True:
            outline+="{} ".format(dec)

        if PosXSerFlag==True:
            outline+="{} ".format(PosXSer)

        if ErPosXSerFlag==True:
            outline+="{} ".format(ErPosXSer)

        if PosYSerFlag==True:
            outline+="{} ".format(PosYSer)

        if ErPosYSerFlag==True:
            outline+="{} ".format(ErPosYSer)

        if MagSerFlag==True:
            outline+="{} ".format(MagSer)

        if ErMagSerFlag==True:
            outline+="{} ".format(ErMagSer)

        if ReSerFlag==True:
            outline+="{} ".format(ReSer)

        if ErReSerFlag==True:
            outline+="{} ".format(ErReSer)

        if NSerFlag==True:
            outline+="{} ".format(NSer)

        if ErNSerFlag==True:
            outline+="{} ".format(ErNSer)

        if AxisSerFlag==True:
            outline+="{} ".format(AxisSer)

        if ErAxisSerFlag==True:
            outline+="{} ".format(ErAxisSer)

        if PaSerFlag==True:
            outline+="{} ".format(PaSer)

        if ErPaSerFlag==True:
            outline+="{} ".format(ErPaSer)

        if KSerFlag==True:
            outline+="{} ".format(KSer)

        if MeanMeSerFlag==True:
            outline+="{} ".format(MeanMeSer)

        if MeSerFlag==True:
            outline+="{} ".format(MeSer)

        if ErMeanMeSerFlag==True:
            outline+="{} ".format(ErMeanMeSer)

        if ErMeSerFlag==True:
            outline+="{} ".format(ErMeSer)

        if PosXExpFlag==True:
            outline+="{} ".format(PosXExp)

        if ErPosXExpFlag==True:
            outline+="{} ".format(ErPosXExp)

        if PosYExpFlag==True:
            outline+="{} ".format(PosYExp)

        if ErPosYExpFlag==True:
            outline+="{} ".format(ErPosYExp)

        if MagExpFlag==True:
            outline+="{} ".format(MagExp)

        if ErMagExpFlag==True:
            outline+="{} ".format(ErMagExp)

        if RsExpFlag==True:
            outline+="{} ".format(RsExp)

        if ErRsExpFlag==True:
            outline+="{} ".format(ErRsExp)

        if AxisExpFlag==True:
            outline+="{} ".format(AxisExp)

        if ErAxisExpFlag==True:
            outline+="{} ".format(ErAxisExp)

        if PaExpFlag==True:
            outline+="{} ".format(PaExp)

        if ErPaExpFlag==True:
            outline+="{} ".format(ErPaExp)

        if MeanMeExpFlag==True:
            outline+="{} ".format(MeanMeExp)

        if MeExpFlag==True:
            outline+="{} ".format(MeExp)


        if ErMeanMeExpFlag==True:
            outline+="{} ".format(ErMeanMeExp)

        if ErMeExpFlag==True:
            outline+="{} ".format(ErMeExp)


        if MsExpFlag==True:
            outline+="{} ".format(MsExp)

        if ErMsExpFlag==True:
            outline+="{} ".format(ErMsExp)


        if SkyGalFlag==True:
            outline+="{} ".format(Sky)

        if ErSkyFlag==True:
            outline+="{} ".format(ErSky)

        if MagTotalFlag==True:
            outline+="{} ".format(MagTotal)

        if ErMagTotalFlag==True:
            outline+="{} ".format(ErMagTotal)

        if BulgeTotalFlag==True:
            outline+="{} ".format(BulgeTotal)

        if ErBtFlag==True:
            outline+="{} ".format(ErBt)

        if ChiNuFlag==True:
            outline+="{} ".format(ChiNu)

        if TidalFlag==True:
            outline+="{} ".format(Tidal)

        if ObjChiNuFlag==True:
            outline+="{} ".format(ObjChiNu)

        if BumpFlag==True:
            outline+="{} ".format(Bump)

        if SNRFlag==True:
            outline+="{} ".format(SNR)

        if NDofFlag==True:
            outline+="{} ".format(NDof)

        if FitFlagFlag==True:
            outline+="{} ".format(FitFlag)

        if ReFitFlag==True:
            outline+="{} ".format(ReFit)

        outline+="\n"


        OUT.write(outline)



    IN.close()
    OUT.close()
    PARIN.close()

    return True




def Ascii2Fits(infile, filefits):
    "Converts an ascii table file to another fits table file"

    if os.path.exists(filefits):
        print("{}: filename exist! it will be overwritten ".format(filefits))
        os.remove(filefits)

    ban = 1


#>>> from astropy.table import Table
#>>> t = Table([[1, 2], [4, 5], [7, 8]], names=('a', 'b', 'c'))
#>>> t.write('table.fits', format='fits')

#>>> from astropy.io import fits
#>>> import numpy as np
#>>> c1 = fits.Column(name='a', array=np.array([1, 2]), format='K')
#>>> c2 = fits.Column(name='b', array=np.array([4, 5]), format='K')
#>>> c3 = fits.Column(name='c', array=np.array([7, 8]), format='K')
#>>> t = fits.BinTableHDU.from_columns([c1, c2, c3])
#>>> t.writeto('table.fits')

    (name,ext,dgtrash)=infile.split(".")

    if ext == "ser":
        fitfunc = "Sersic"

    elif ext == "bd":
        fitfunc = "BD"
    else:
        print("Ascii2Table: file extention not recognized\n")
        sys.exit()


    if (fitfunc == "BD"):
# 79 cols

        (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, PosXExp, ErPosXExp, PosYExp, ErPosYExp, MagExp, ErMagExp, RsExp, ErRsExp, AxisExp, ErAxisExp, PaExp, ErPaExp, MeanMeExp, ErMeanMeExp, MeExp, ErMeExp, MsExp, ErMsExp, Sky, ErSky, MagTotal, ErMagTotal, BulgeTotal, ErBt, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit) = np.genfromtxt(infile, delimiter="", unpack=True)


        col1 = fits.Column(name='SExtractor Num Object', format='J', array=Num)
        col2 = fits.Column(name='SExtractor RA', format='E', array=RA)
        col3 = fits.Column(name='SExtractor DEC', format='E', array=Dec)
        col4 = fits.Column(name='SExtractor Posx', format='E', array=XPos)
        col5 = fits.Column(name='SExtractor Posy', format='E', array=YPos)
        col6 = fits.Column(name='SExtractor Mag', format='E', array=Mag)
        col7 = fits.Column(name='SExtractor KronRadius', format='E', array=Kron)
        col8 = fits.Column(name='SExtractor FluxRadius', format='E', array=FluxRad)
        col9 = fits.Column(name='SExtractor IsoArea', format='E', array=IsoArea)
        col10 = fits.Column(name='SExtractor AImage', format='E', array=AIm)

        col11 = fits.Column(name='SExtractor AxisRatio', format='E', array=AR)
        col12 = fits.Column(name='SExtractor PositionAngle', format='E', array=Angle)
        col13 = fits.Column(name='SExtractor Background', format='E', array=Background)
        col14 = fits.Column(name='SExtractor ClassStar', format='E', array=Class)
        col15 = fits.Column(name='SExtractor Flag', format='J', array=Flag)
        col16 = fits.Column(name='Xmin', format='E', array=XMin)
        col17 = fits.Column(name='Xmax', format='E', array=XMax)
        col18 = fits.Column(name='Ymin', format='E', array=YMin)
        col19 = fits.Column(name='Ymax', format='E', array=YMax)
        col20 = fits.Column(name='SersicInput', format='E', array=Sersic)

        col21 = fits.Column(name='SkyKronRadius', format='E', array=RSky)
        col22 = fits.Column(name='KronRadius', format='E', array=RKron)
        col23 = fits.Column(name='SkyFlag', format='J', array=SkyFlag)
        col24 = fits.Column(name='OffX', format='J', array=OFFX)
        col25 = fits.Column(name='OffY', format='J', array=OFFY)
#        col26 = fits.Column(name='OutputFlag', format='J', array=lflag)
        col26 = fits.Column(name='Input SubImage', format='20A', array=InputImage)
        col27 = fits.Column(name='Output Image', format='40A', array=OutIm)
        col28 = fits.Column(name='Object File', format='10A', array=Objx)
        col29 = fits.Column(name='Num Object', format='J', array=PPPnum)
        col30 = fits.Column(name='Restart File', format='10A', array=Restart)


        col31 = fits.Column(name='Input Image', format='10A', array=Cluster)
        col32 = fits.Column(name='RA', format='30A', array=ra)
        col33 = fits.Column(name='DEC', format='30A', array=dec)
        col34 = fits.Column(name='PosxSer', format='E', array=PosXSer)
        col35 = fits.Column(name='PosxSerErr', format='E', array=ErPosXSer)
        col36 = fits.Column(name='PosySer', format='E', array=PosYSer)
        col37 = fits.Column(name='PosySerErr', format='E', array=ErPosYSer)
        col38 = fits.Column(name='MagSer', format='E', array=MagSer)
        col39 = fits.Column(name='MagSerErr', format='E', array=ErMagSer)
        col40 = fits.Column(name='ReSer', format='E', array=ReSer)
        col41 = fits.Column(name='ReSerErr', format='E', array=ErReSer)
        col42 = fits.Column(name='IndexSer', format='E', array=NSer)

        col43 = fits.Column(name='IndexSerErr', format='E', array=ErNSer)
        col44 = fits.Column(name='AxisRatioSer', format='E', array=AxisSer)
        col45 = fits.Column(name='AxisRatioSerErr', format='E', array=ErAxisSer)
        col46 = fits.Column(name='PositionAngleSer', format='E', array=PaSer)
        col47 = fits.Column(name='PositionAngleSerErr', format='E', array=ErPaSer)
        col48 = fits.Column(name='KSer', format='E', array=KSer)
        col49 = fits.Column(name='MeanMeSer', format='E', array=MeanMeSer)
        col50 = fits.Column(name='ErMeanMeSer', format='E', array=ErMeanMeSer)

        col51 = fits.Column(name='MeSer', format='E', array=MeSer)
        col52 = fits.Column(name='ErMeSer', format='E', array=ErMeSer)
        col53 = fits.Column(name='PosXExp', format='E', array=PosXExp)
        col54 = fits.Column(name='PosxExpErr', format='E', array=ErPosXExp)
        col55 = fits.Column(name='PosyExp', format='E', array=PosYExp)
        col56 = fits.Column(name='PosyExpErr', format='E', array=ErPosYExp)
        col57 = fits.Column(name='MagExp', format='E', array=MagExp)
        col58 = fits.Column(name='MagExpErr', format='E', array=ErMagExp)
        col59 = fits.Column(name='RsExp', format='E', array=RsExp)
        col60 = fits.Column(name='RsExpErr', format='E', array=ErRsExp)


        col61 = fits.Column(name='AxisRatioExp', format='E', array=AxisExp)
        col62 = fits.Column(name='AxisRatioExpErr', format='E', array=ErAxisExp)
        col63 = fits.Column(name='PositionAngleExp', format='E', array=PaExp)
        col64 = fits.Column(name='PositionAngleExpErr', format='E', array=ErPaExp)
        col65 = fits.Column(name='MeanMeExp', format='E', array=MeanMeExp)
        col66 = fits.Column(name='ErMeanMeExp', format='E', array=ErMeanMeExp)
        col67 = fits.Column(name='MeExp', format='E', array=MeExp)
        col68 = fits.Column(name='ErMeExp', format='E', array=ErMeExp)
        col69 = fits.Column(name='MsExp', format='E', array=MsExp)
        col70 = fits.Column(name='ErMsExp', format='E', array=ErMsExp)

        col71 = fits.Column(name='Sky', format='E', array=Sky)
        col72 = fits.Column(name='SkyErr', format='E', array=ErSky)
        col73 = fits.Column(name='MagTotal', format='E', array=MagTotal)
        col74 = fits.Column(name='MagTotalErr', format='E', array=ErMagTotal)
        col75 = fits.Column(name='BulgeTotal', format='E', array=BulgeTotal)
        col76 = fits.Column(name='BulgeTotalErr', format='E', array=ErBt)
        col77 = fits.Column(name='ChiNu', format='E', array=ChiNu)
        col78 = fits.Column(name='Tidal', format='E', array=Tidal)
        col79 = fits.Column(name='LocalChiNu', format='E', array=ObjChiNu)
        col80 = fits.Column(name='Bumpiness', format='E', array=Bump)

        col81 = fits.Column(name='Signal to Noise Ratio', format='E', array=SNR)
        col82 = fits.Column(name='Num Degrees Of Freedom', format='E', array=NDof)
        col83 = fits.Column(name='FitFlag', format='E', array=FitFlag)
        col84 = fits.Column(name='ReFit Flag', format='E', array=ReFit)


        t=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,
                         col41, col42, col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, col62, col63, col64, col65, col66, col67, col68, col69, col70, col71, col72, col73, col74, col75, col76, col77, col78, col79, col80, col81, col82, col83, col84])



    else:
#64 cols

        (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, Sky, ErSky, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit)= np.genfromtxt(infile, delimiter="", unpack=True)


        col1 = fits.Column(name='SExtractor Num Object', format='J', array=Num)
        col2 = fits.Column(name='SExtractor RA', format='E', array=RA)
        col3 = fits.Column(name='SExtractor DEC', format='E', array=Dec)
        col4 = fits.Column(name='SExtractor Posx', format='E', array=XPos)
        col5 = fits.Column(name='SExtractor Posy', format='E', array=YPos)
        col6 = fits.Column(name='SExtractor Mag', format='E', array=Mag)
        col7 = fits.Column(name='SExtractor KronRadius', format='E', array=Kron)
        col8 = fits.Column(name='SExtractor FluxRadius', format='E', array=FluxRad)
        col9 = fits.Column(name='SExtractor IsoArea', format='E', array=IsoArea)
        col10 = fits.Column(name='SExtractor AImage', format='E', array=AIm)

        col11 = fits.Column(name='SExtractor AxisRatio', format='E', array=AR)
        col12 = fits.Column(name='SExtractor PositionAngle', format='E', array=Angle)
        col13 = fits.Column(name='SExtractor Background', format='E', array=Background)
        col14 = fits.Column(name='SExtractor ClassStar', format='E', array=Class)
        col15 = fits.Column(name='SExtractor Flag', format='J', array=Flag)
        col16 = fits.Column(name='Xmin', format='E', array=XMin)
        col17 = fits.Column(name='Xmax', format='E', array=XMax)
        col18 = fits.Column(name='Ymin', format='E', array=YMin)
        col19 = fits.Column(name='Ymax', format='E', array=YMax)
        col20 = fits.Column(name='SersicInput', format='E', array=Sersic)


        col21 = fits.Column(name='SkyKronRadius', format='E', array=RSky)
        col22 = fits.Column(name='KronRadius', format='E', array=RKron)
        col23 = fits.Column(name='SkyFlag', format='J', array=SkyFlag)
        col24 = fits.Column(name='OffX', format='J', array=OFFX)
        col25 = fits.Column(name='OffY', format='J', array=OFFY)
#        col26 = fits.Column(name='OutputFlag', format='J', array=lflag)
        col26 = fits.Column(name='Input SubImage', format='20A', array=InputImage)
        col27 = fits.Column(name='Output Image', format='40A', array=OutIm)
        col28 = fits.Column(name='Object File', format='10A', array=Objx)
        col29 = fits.Column(name='Num Object', format='J', array=PPPnum)
        col30 = fits.Column(name='Restart File', format='10A', array=Restart)


        col31 = fits.Column(name='Input Image', format='10A', array=Cluster)
        col32 = fits.Column(name='RA', format='30A', array=ra)
        col33 = fits.Column(name='DEC', format='30A', array=dec)
        col34 = fits.Column(name='PosxSer', format='E', array=PosXSer)
        col35 = fits.Column(name='PosxSerErr', format='E', array=ErPosXSer)
        col36 = fits.Column(name='PosySer', format='E', array=PosYSer)
        col37 = fits.Column(name='PosySerErr', format='E', array=ErPosYSer)
        col38 = fits.Column(name='MagSer', format='E', array=MagSer)
        col39 = fits.Column(name='MagSerErr', format='E', array=ErMagSer)
        col40 = fits.Column(name='ReSer', format='E', array=ReSer)

        col41 = fits.Column(name='ReSerErr', format='E', array=ErReSer)
        col42 = fits.Column(name='IndexSer', format='E', array=NSer)
        col43 = fits.Column(name='IndexSerErr', format='E', array=ErNSer)
        col44 = fits.Column(name='AxisRatioSer', format='E', array=AxisSer)
        col45 = fits.Column(name='AxisRatioSerErr', format='E', array=ErAxisSer)
        col46 = fits.Column(name='PositionAngleSer', format='E', array=PaSer)
        col47 = fits.Column(name='PositionAngleSerErr', format='E', array=ErPaSer)
        col48 = fits.Column(name='KSer', format='E', array=KSer)
        col49 = fits.Column(name='MeanMeSer', format='E', array=MeanMeSer)
        col50 = fits.Column(name='ErMeanMeSer', format='E', array=ErMeanMeSer)


        col51 = fits.Column(name='MeSer', format='E', array=MeSer)
        col52 = fits.Column(name='ErMeSer', format='E', array=ErMeSer)
        col53 = fits.Column(name='Sky', format='E', array=Sky)
        col54 = fits.Column(name='SkyErr', format='E', array=ErSky)
        col55 = fits.Column(name='ChiNu', format='E', array=ChiNu)
        col56 = fits.Column(name='Tidal', format='E', array=Tidal)
        col57 = fits.Column(name='LocalChiNu', format='E', array=ObjChiNu)
        col58 = fits.Column(name='Bumpiness', format='E', array=Bump)
        col59 = fits.Column(name='Signal to Noise Ratio', format='E', array=SNR)
        col60 = fits.Column(name='Num Degrees Of Freedom', format='E', array=NDof)

        col61 = fits.Column(name='FitFlag', format='E', array=FitFlag)
        col62 = fits.Column(name='ReFit Flag', format='E', array=ReFit)


        t=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,
                         col41, col42, col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, col62])




# add comments to header

#tbhdu.header.add_comment('Flags value 0: good fit')
#tbhdu.header.add_comment('Flags value 1: ran into constraints during the fit')
#tbhdu.header.add_comment('Flags value 2: parameters values makes non sense ')
#tbhdu.header.add_comment('Flags value 4: need extra component (bar) ')
#tbhdu.header.add_comment('Flags value 8: bad fit')
#tbhdu.header.add_comment('Flags value 16: cD galaxy ')
#tbhdu.header.add_comment('Flags value 32: fit models are not compatible with galaxy. Example: galaxy distorsioned or need model to fit spiral arms.')
#tbhdu.header.add_comment('Flags value 64: Delete it. Probably a star, galaxy too faint or has saturated pixels')

    t.writeto(filefits)  # writing output

    return True







#############################################
#############################################
#############################################


def OldAscii2Table(file, filefits):
    "Converts ascii table to fits table. Deprecated"


    if os.path.exists(filefits):
        print("{}: filename exist! it will be overwritten ".format(filefits))
        os.remove(filefits)

    ban = 1


#>>> from astropy.table import Table
#>>> t = Table([[1, 2], [4, 5], [7, 8]], names=('a', 'b', 'c'))
#>>> t.write('table.fits', format='fits')

#>>> from astropy.io import fits
#>>> import numpy as np
#>>> c1 = fits.Column(name='a', array=np.array([1, 2]), format='K')
#>>> c2 = fits.Column(name='b', array=np.array([4, 5]), format='K')
#>>> c3 = fits.Column(name='c', array=np.array([7, 8]), format='K')
#>>> t = fits.BinTableHDU.from_columns([c1, c2, c3])
#>>> t.writeto('table.fits')



    f = open(file)
    try:
        for line in f:
            (num, ra, dec, xpos, ypos, mag, kron, fluxrad, isoarea, aim, ar, angle, background, star, sexflag, xmin, xmax, ymin, ymax, sersic, rsky, rkron, skyflag, offx, offy, flag, inputimage, outim, objx, restart, cluster, posxser, erposxser, posyser, erposyser, magser, ermagser, reser, erreser, nser, ernser, axiser, eraxiser, paser, erpaser,
         kser, meanmeser, meser, posxexp, erposxexp, posyexp, erposyexp, magexp, ermagexp, rsexp, errsexp, axisexp, eraxisexp, paexp, erpaexp, meanmeexp, meexp, sky, ersky, magtotal, ermagtotal, bulgetotal, erbt, chinu, tidal, objchinu, bump, meanmesky, snr, ndof, fitflag, ermagsersim, erresersim, ernsersim, ermagexpsim, errsexpsim) = line.split()

        if ban == 1:
            ban = 0
            lnum = [int(num)]
            lra = [ra]
            ldec = [dec]
            lxpos = [float(xpos)]
            lypos = [float(ypos)]
            lmag = [float(mag)]
            lkron = [float(kron)]
            lfluxrad = [float(fluxrad)]
            lisoarea = [float(isoarea)]
            laim = [float(aim)]
            lar = [float(ar)]
            langle = [float(angle)]
            lbackground = [float(background)]
            lstar = [float(star)]
            lsexflag = [int(sexflag)]
            lxmin = [float(xmin)]
            lxmax = [float(xmax)]
            lymin = [float(ymin)]
            lymax = [float(ymax)]
            lsersic = [float(sersic)]
            lrsky = [float(rsky)]
            lrkron = [float(rkron)]
            lskyflag = [int(skyflag)]
            loffx = [int(offx)]
            loffy = [int(offy)]
            lflag = [int(flag)]
            linputimage = [inputimage]
            loutim = [outim]
            lobjx = [objx]
            lrestart = [restart]  # 30
            lcluster = [cluster]
            lposxser = [float(posxser)]
            lerposxser = [float(erposxser)]
            lposyser = [float(posyser)]
            lerposyser = [float(erposyser)]
            lmagser = [float(magser)]
            lermagser = [float(ermagser)]
            lreser = [float(reser)]
            lerreser = [float(erreser)]
            lnser = [float(nser)]
            lernser = [float(ernser)]
            laxiser = [float(axiser)]
            leraxiser = [float(eraxiser)]
            lpaser = [float(paser)]
            lerpaser = [float(erpaser)]
            lkser = [float(kser)]
            lmeanmeser = [float(meanmeser)]
            lmeser = [float(meser)]
            lposxexp = [float(posxexp)]
            lerposxexp = [float(erposxexp)]  # 50
            lposyexp = [float(posyexp)]
            lerposyexp = [float(erposyexp)]
            lmagexp = [float(magexp)]
            lermagexp = [float(ermagexp)]
            lrsexp = [float(rsexp)]
            lerrsexp = [float(errsexp)]
            laxisexp = [float(axisexp)]
            leraxisexp = [float(eraxisexp)]
            lpaexp = [float(paexp)]
            lerpaexp = [float(erpaexp)]
            lmeanmeexp = [float(meanmeexp)]
            lmeexp = [float(meexp)]
            lsky = [float(sky)]
            lersky = [float(ersky)]
            lmagtotal = [float(magtotal)]
            lermagtotal = [float(ermagtotal)]
            lbulgetotal = [float(bulgetotal)]
            lerbt = [float(erbt)]
            lchinu = [float(chinu)]
            ltidal = [float(tidal)]    # 70
            lobjchinu = [float(objchinu)]
            lbump = [float(bump)]
            lmeanmesky = [float(meanmesky)]
            lsnr = [float(snr)]
            lndof = [float(ndof)]
            lfitflag = [float(fitflag)]
            lermagsersim = [float(ermagsersim)]
            lerresersim = [float(erresersim)]
            lernsersim = [float(ernsersim)]
            lermagexpsim = [float(ermagexpsim)]
            lerrsexpsim = [float(errsexpsim)]

        else:

            lnum.append(int(num))
            lra.append(ra)
            ldec.append(dec)
            lxpos.append(float(xpos))
            lypos.append(float(ypos))
            lmag.append(float(mag))
            lkron.append(float(kron))
            lfluxrad.append(float(fluxrad))
            lisoarea.append(float(isoarea))
            laim.append(float(aim))
            lar.append(float(ar))
            langle.append(float(angle))
            lbackground.append(float(background))
            lstar.append(float(star))
            lsexflag.append(int(sexflag))
            lxmin.append(float(xmin))
            lxmax.append(float(xmax))
            lymin.append(float(ymin))
            lymax.append(float(ymax))
            lsersic.append(float(sersic))
            lrsky.append(float(rsky))
            lrkron.append(float(rkron))
            lskyflag.append(int(skyflag))
            loffx.append(int(offx))
            loffy.append(int(offy))
            lflag.append(int(flag))
            linputimage.append(inputimage)
            loutim.append(outim)
            lobjx.append(objx)
            lrestart.append(restart)
            lcluster.append(cluster)
            lposxser.append(float(posxser))
            lerposxser.append(float(erposxser))
            lposyser.append(float(posyser))
            lerposyser.append(float(erposyser))
            lmagser.append(float(magser))
            lermagser.append(float(ermagser))
            lreser.append(float(reser))
            lerreser.append(float(erreser))
            lnser.append(float(nser))
            lernser.append(float(ernser))
            laxiser.append(float(axiser))
            leraxiser.append(float(eraxiser))
            lpaser.append(float(paser))
            lerpaser.append(float(erpaser))
            lkser.append(float(kser))
            lmeanmeser.append(float(meanmeser))
            lmeser.append(float(meser))
            lposxexp.append(float(posxexp))
            lerposxexp.append(float(erposxexp))
            lposyexp.append(float(posyexp))
            lerposyexp.append(float(erposyexp))
            lmagexp.append(float(magexp))
            lermagexp.append(float(ermagexp))
            lrsexp.append(float(rsexp))
            lerrsexp.append(float(errsexp))
            laxisexp.append(float(axisexp))
            leraxisexp.append(float(eraxisexp))
            lpaexp.append(float(paexp))
            lerpaexp.append(float(erpaexp))
            lmeanmeexp.append(float(meanmeexp))
            lmeexp.append(float(meexp))
            lsky.append(float(sky))
            lersky.append(float(ersky))
            lmagtotal.append(float(magtotal))
            lermagtotal.append(float(ermagtotal))
            lbulgetotal.append(float(bulgetotal))
            lerbt.append(float(erbt))
            lchinu.append(float(chinu))
            ltidal.append(float(tidal))
            lobjchinu.append(float(objchinu))
            lbump.append(float(bump))
            lmeanmesky.append(float(meanmesky))
            lsnr.append(float(snr))
            lndof.append(float(ndof))
            lfitflag.append(float(fitflag))
            lermagsersim.append(float(ermagsersim))
            lerresersim.append(float(erresersim))
            lernsersim.append(float(ernsersim))
            lermagexpsim.append(float(ermagexpsim))
            lerrsexpsim.append(float(errsexpsim))

    finally:
            f.close()

#$Num[$i] $RA[$i] $Dec[$i] $XPos[$i] $YPos[$i] $Mag[$i] $Kron[$i] $FluxRad[$i] $IsoArea[$i] $AIm[$i] $AR[$i] $Angle[$i] $Background[$i] $Star[$i] $Flag[$i] $XMin[$i] $XMax[$i] $YMin[$i] $YMax[$i] $Sersic[$i] $RSky[$i] $RKron[$i] $SkyFlag[$i] $OFFX[$i] $OFFY[$i] 1 $inputimage $outim $object $restart $cluster $posxser $erposxser $posyser $erposyser $magser $ermagser $reser $erreser $nser $ernser $axiser $eraxiser $paser $erpaser $kser $meanmeser $meser $posxexp $erposxexp $posyexp $erposyexp $magexp $ermagexp $rsexp $errsexp $axisexp $eraxisexp $paexp $erpaexp $meanmeexp $meexp $sky $ersky $magtotal $ermagtotal $bulgetotal $erbt $chinu $tidal $objchinu $bump $meanmesky $snr $ndof $fitflag $ermagsersim $erresersim $ernsersim $ermagexpsim $errsexpsim

    col1 = fits.Column(name='NumObject', format='J', array=lnum)
    col2 = fits.Column(name='RA', format='30A', array=lra)
    col3 = fits.Column(name='DEC', format='30A', array=ldec)
    col4 = fits.Column(name='SExtractorPosx', format='E', array=lxpos)
    col5 = fits.Column(name='SExtractorPosy', format='E', array=lypos)
    col6 = fits.Column(name='SExtractorMag', format='E', array=lmag)
    col7 = fits.Column(name='SExtractorKronRadius', format='E', array=lkron)
    col8 = fits.Column(name='SExtractorFluxRadius',
                       format='E', array=lfluxrad)
    col9 = fits.Column(name='SExtractorIsoArea', format='E', array=lisoarea)
    col10 = fits.Column(name='SExtractorAImage', format='E', array=laim)
    col11 = fits.Column(name='SExtractorAxisRatio', format='E', array=lar)
    col12 = fits.Column(name='SExtractorPositionAngle',
                        format='E', array=langle)
    col13 = fits.Column(name='SExtractorBackground',
                        format='E', array=lbackground)
    col14 = fits.Column(name='SExtractorClassStar', format='E', array=lstar)
    col15 = fits.Column(name='SExtractorFlag', format='J', array=lsexflag)
    col16 = fits.Column(name='Xmin', format='E', array=lxmin)
    col17 = fits.Column(name='Xmax', format='E', array=lxmax)
    col18 = fits.Column(name='Ymin', format='E', array=lymin)
    col19 = fits.Column(name='Ymax', format='E', array=lymax)
    col20 = fits.Column(name='SersicInput', format='E', array=lsersic)
    col21 = fits.Column(name='SkyKronRadius', format='E', array=lrsky)
    col22 = fits.Column(name='KronRadius', format='E', array=lrkron)
    col23 = fits.Column(name='SkyFlag', format='J', array=lskyflag)
    col24 = fits.Column(name='OffX', format='J', array=loffx)
    col25 = fits.Column(name='OffY', format='J', array=loffy)
    col26 = fits.Column(name='OutputFlag', format='J', array=lflag)
    col27 = fits.Column(name='InputSubImage',
                        format='20A', array=linputimage)
    col28 = fits.Column(name='OutputImage', format='40A', array=loutim)
    col29 = fits.Column(name='ObjectFile', format='10A', array=lobject)
    col30 = fits.Column(name='RestartFile', format='10A', array=lrestart)
    col31 = fits.Column(name='InputImage', format='10A', array=lcluster)
    col32 = fits.Column(name='PosxSer', format='E', array=lposxser)
    col33 = fits.Column(name='PosxSerErr', format='E', array=lerposxser)
    col34 = fits.Column(name='PosySer', format='E', array=lposyser)
    col35 = fits.Column(name='PosySerErr', format='E', array=lerposyser)
    col36 = fits.Column(name='MagSer', format='E', array=lmagser)
    col37 = fits.Column(name='MagSerErr', format='E', array=lermagser)
    col38 = fits.Column(name='ReSer', format='E', array=lreser)
    col39 = fits.Column(name='ReSerErr', format='E', array=lerreser)
    col40 = fits.Column(name='IndexSer', format='E', array=lnser)
    col41 = fits.Column(name='IndexSerErr', format='E', array=lernser)
    col42 = fits.Column(name='AxisRatioSer', format='E', array=laxiser)
    col43 = fits.Column(name='AxisRatioSerErr', format='E', array=leraxiser)
    col44 = fits.Column(name='PositionAngleSer', format='E', array=lpaser)
    col45 = fits.Column(name='PositionAngleSerErr',
                        format='E', array=lerpaser)
    col46 = fits.Column(name='KSer', format='E', array=lkser)
    col47 = fits.Column(name='MeanMeSer', format='E', array=lmeanmeser)
    col48 = fits.Column(name='MeSer', format='E', array=lmeser)
    col49 = fits.Column(name='PosxExp', format='E', array=lposxexp)
    col50 = fits.Column(name='PosxExpErr', format='E', array=lerposxexp)
    col51 = fits.Column(name='PosyExp', format='E', array=lposyexp)
    col52 = fits.Column(name='PosyExpErr', format='E', array=lerposyexp)
    col53 = fits.Column(name='MagExp', format='E', array=lmagexp)
    col54 = fits.Column(name='MagExpErr', format='E', array=lermagexp)
    col55 = fits.Column(name='RsExp', format='E', array=lrsexp)
    col56 = fits.Column(name='RsExpErr', format='E', array=lerrsexp)
    col57 = fits.Column(name='AxisRatioExp', format='E', array=laxisexp)
    col58 = fits.Column(name='AxisRatioExpErr', format='E', array=leraxisexp)
    col59 = fits.Column(name='PositionAngleExp', format='E', array=lpaexp)
    col60 = fits.Column(name='PositionAngleExpErr',
                        format='E', array=lerpaexp)
    col61 = fits.Column(name='MeanMeExp', format='E', array=lmeanmeexp)
    col62 = fits.Column(name='MeExp', format='E', array=lmeexp)
    col63 = fits.Column(name='Sky', format='E', array=lsky)
    col64 = fits.Column(name='SkyErr', format='E', array=lersky)
    col65 = fits.Column(name='MagTotal', format='E', array=lmagtotal)
    col66 = fits.Column(name='MagTotalErr', format='E', array=lermagtotal)
    col67 = fits.Column(name='BulgeTotal', format='E', array=lbulgetotal)
    col68 = fits.Column(name='BulgeTotalErr', format='E', array=lerbt)
    col69 = fits.Column(name='ChiNu', format='E', array=lchinu)
    col70 = fits.Column(name='Tidal', format='E', array=ltidal)
    col71 = fits.Column(name='LocalChiNu', format='E', array=lobjchinu)
    col72 = fits.Column(name='Bumpiness', format='E', array=lbump)
    col73 = fits.Column(name='MeanMeSky', format='E', array=lmeanmesky)
    col74 = fits.Column(name='SignalNoiseRatio', format='E', array=lsnr)
    col75 = fits.Column(name='NumDegreesOfFreedom', format='E', array=lndof)
    col76 = fits.Column(name='FitFlag', format='E', array=lfitflag)
    col77 = fits.Column(name='MagSerSimErr', format='E', array=lermagsersim)
    col78 = fits.Column(name='ReSerSimErr', format='E', array=lerresersim)
    col79 = fits.Column(name='IndexSerSimErr', format='E', array=lernsersim)
    col80 = fits.Column(name='MagExpSimErr', format='E', array=lermagexpsim)
    col81 = fits.Column(name='RsExpSimErr', format='E', array=lerrsexpsim)

    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,
                         col41, col42, col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, col62, col63, col64, col65, col66, col67, col68, col69, col70, col71, col72, col73, col74, col75, col76, col77, col78, col79, col80, col81])

    tbhdu = fits.BinTableHDU.from_columns(cols)


# add comments to header

#tbhdu.header.add_comment('Flags value 0: good fit')
#tbhdu.header.add_comment('Flags value 1: ran into constraints during the fit')
#tbhdu.header.add_comment('Flags value 2: parameters values makes non sense ')
#tbhdu.header.add_comment('Flags value 4: need extra component (bar) ')
#tbhdu.header.add_comment('Flags value 8: bad fit')
#tbhdu.header.add_comment('Flags value 16: cD galaxy ')
#tbhdu.header.add_comment('Flags value 32: fit models are not compatible with galaxy. Example: galaxy distorsioned or need model to fit spiral arms.')
#tbhdu.header.add_comment('Flags value 64: Delete it. Probably a star, galaxy too faint or has saturated pixels')

    tbhdu.writeto(filefits)  # writing output

    return True
