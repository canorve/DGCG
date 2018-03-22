
import os.path
import os
from astropy.io import fits
import numpy as np
import sys
import scipy
import scipy.special
import subprocess as sp

from lib import check




###############
# Deprecated. Do not use it anymore  #
###############
def Default():
    # k Check
    "initialize default values for all the variables (including configuration file)"
    "All these variables are global"


#   in variables
#    global Num, RA,Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle
#    global Background, Sky, Class, Flag, XMin, XMax, YMin, YMax, XSMin , XSMax, YSMin
#    global YSMax, Sersic, RSky, RKron, SkyFlag, Neighbors, IX, IY, XBuffer,YBuffer, OFFX,OFFY



# parameters file

    global Img
    global SexCat
    global SigImg
    global PsfDir
    global MagZpt
    global PlateScale
    global FitFunc
    global GalBot
    global GalTop
    global ConvBox
    global FitBox
    global MagDiff
    global MagCut
    global KronScale
    global SkyScale
    global Offset
    global SkyWidth
    global NSer
    global MaxFit
#    global MagMin
    global MagMax
    global FlagSex
    global ConsFile
    global Region
    global Bxmin
    global Bymin
    global Bxmax
    global Bymax
    global Split
#    global AutoSatRegion
    global SatRegionScale
    global Ds9SatReg
    global FileOut
    global SegFile
    global SkyFile
    global PixFile
    global PixPrefix
    global Ds9OutName
    global Ds9OutNum
    global Ds9FitReg
    global BoxOut
    global BoxSkyOut
    global SexSort
    global SexArSort
    global Erase
    global Nice
    global Overwrite
    global Execute
    global Contrast
    global Bias
    global SimulFit
    global ErrPer

# other variables

    global TempDir
    global PixDir
    global MaskPixDir
    global SkyDir
    global OutputDir
    global InputDir
    global OffsetPos  # offset position file
    global Crashes
    global Fitted
    global SkyCrashes
    global SkyFitted

    global NCol
    global NRow
    global Buffer
    global Failures
    global Success
    global Total
    global TotPix  # enough pixels to compute sky?
    global MinPix  # minimum pixels for objx
    global MaxTimes
    global ExpTime
    global Gain
    global Rdnoise
    global TimesRe
    global ListObjs
# outputs files

    global LogFile
    global BtFile
    global FitsFile

    Img = "image.fits"
    SexCat = "sex.cat"
    SigImg = "none"
    PsfDir = "psfdir"
    MagZpt = 21.471
    PlateScale = 1.0
    FitFunc = "BD"
    GalBot = 0
    GalTop = 0.5
    ConvBox = 60
    FitBox = 6
    MagDiff = 5
    MagCut = 3
    KronScale = 1.5
    SkyScale = 3
    Offset = 20
    SkyWidth = 20
    NSer = 1.5
    MaxFit = 10
#    MagMin = 0
    MagMax = 19
    FlagSex = 4
    ConsFile = "constraints"
    Region = 0
    Bxmin = 0
    Bymin = 0
    Bxmax = 2000
    Bymax = 2000
    Split = 5
#    AutoSatRegion = 0
    SatRegionScale = 2
    Ds9SatReg = "sat.reg"
    FileOut = "fits"
    SegFile = "seg.fits"
    SkyFile = "sky.fits"
    PixFile = "pix.fits"
    PixPrefix = "pixels"
    Ds9OutName = "images"
    Ds9OutNum = 40
    Ds9FitReg = "fit.reg"
    BoxOut = "box.reg"
    BoxSkyOut = "boxsky.reg"
    SexSort = "sexsort.cat"
    SexArSort = "sexarsort.cat"
    Erase = 0
    Nice = 0
    Overwrite = 1
    Execute = 1
    Contrast = 3
    Bias = 0.7
    SimulFit = True
    ErrPer = 0.3

# other variables

    TempDir = "tempfits"
    PixDir = "objpixels"
    MaskPixDir = "maskobjpixels"
    SkyDir = "objsky"
    OutputDir = "outputs"
    InputDir = "inputs"
    OffsetPos = "offsets"  # offset position file
    Crashes = "dgcg.crashes"
    Fitted = "dgcg.fitted"
    SkyCrashes = "sky.crashes"
    SkyFitted = "sky.fitted"

    NCol = 2000
    NRow = 2000
    Buffer = 200
    Failures = 0
    Success = 0
    Total = 0
    TotPix = 30  # enough pixels to compute sky?
    MinPix = 10  # minimum pixels for object
    MaxTimes = 3
    ExpTime = 1
    Gain = 1
    Rdnoise = 0
    TimesRe = 3

    ListObjs = "ListObjs"
# outputs files

    LogFile = FileOut + ".log"
    BtFile = FileOut + ".bt"
    FitsFile = FileOut + ".fits"



    return True


def ReadFile(parvar,File):
    "Read DGCG parameters from file"


    with open(File) as INCONF:

        # All lines including the blank ones
        lines = (line.rstrip() for line in INCONF)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            (params) = line.split()

            if params[0] == "Img":     # input image
                try:
                    (param, Img) = line.split()
                    parvar.Img = str(Img)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SexCat":     # Sextractor catalog
                try:
                    (param, SexCat) = line.split()
                    parvar.SexCat = str(SexCat)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SigImg":     # Sigma image
                try:
                    (param, SigImg) = line.split()
                    parvar.SigImg = str(SigImg)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "PsfDir":     # PSF directory
                try:
                    (param, PsfDir) = line.split()
                    parvar.PsfDir = str(PsfDir)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "MagZpt":     # Magnitude zeropoint (GALFIT)
                try:
                    (param, MagZpt) = line.split()
                    parvar.MagZpt = float(MagZpt)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "PlateScale":  # Plate Scale
                try:
                    (param, PlateScale) = line.split()
                    parvar.PlateScale = float(PlateScale)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "FitFunc":     # BD or Sersic decomposition?
                try:
                    (param, FitFunc) = line.split()
                    parvar.FitFunc = str(FitFunc)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "GalClas":     # Star (1.0) or galaxy (0.0)?
                try:
                    (param, GalClas) = line.split()
                    parvar.GalClas = float(GalClas)
#                    parvar.GalBot = float(GalBot)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "ConvBox":  # convolution box Size
                try:
                    (param, ConvBox) = line.split()
                    parvar.ConvBox = int(ConvBox)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "FitBox":  # times the galaxy Size
                try:
                    (param, FitBox) = line.split()
                    parvar.FitBox = float(FitBox)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            # Maximum magnitude difference among galaxies for simultaneous
            # fitting
            if params[0] == "MagDiff":
                try:
                    (param, MagDiff) = line.split()
                    parvar.MagDiff = float(MagDiff)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            # Maximum magnitude difference between galaxy and low mag limit for
            # simultaneous fitting
            if params[0] == "MagCut":
                try:
                    (param, MagCut) = line.split()
                    parvar.MagCut = float(MagCut)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "KronScale":  # Scale factor for Ellipses
                try:
                    (param, KronScale) = line.split()
                    parvar.KronScale = float(KronScale)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SkyScale":  # Scale factor for Ellipses for sky computing
                try:
                    (param, SkyScale) = line.split()
                    parvar.SkyScale = float(SkyScale)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Offset":  # Additional offset to scale factor
                try:
                    (param, Offset) = line.split()
                    parvar.Offset = int(Offset)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SkyWidth":  # Additional offset to Sky scale factor
                try:
                    (param, SkyWidth) = line.split()
                    parvar.SkyWidth = int(SkyWidth)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "NSer":  # Sersic index parameter
                try:
                    (param, NSer) = line.split()
                    parvar.NSer = float(NSer)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SimulFit":  # Do simultaneous fitting?
                try:
                    (param, SimulFit) = line.split()
                    SimulFit = int(SimulFit)
                    parvar.SimulFit = bool(SimulFit)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "MaxFit":  # Max number of simultaneous fitting
                try:
                    (param, MaxFit) = line.split()
                    parvar.MaxFit = int(MaxFit)
                except ValueError:
                    print >> (sys.stderr, "does not exist")
                    print >> (sys.stderr, "Exception: %s" % str(e))
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    sys.exit(1)

#            if params[0] == "MagRange":     # Acceptable mag range
#                try:
#                    (param, MagMin, MagMax) = line.split()
#                    parvar.MagMin = float(MagMin)
#                    parvar.MagMax = float(MagMax)
#                except:
#                    print("Unexpected error at reading param file:",
#                          sys.exc_info()[0])
#                    raise

            if params[0] == "MagMax":     # Acceptable mag range
                try:
                    (param, MagMax) = line.split()
                    parvar.MagMax = float(MagMax)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "FlagSex":  # Sextractor flags
                try:
                    (param, FlagSex) = line.split()
                    parvar.FlagSex = int(FlagSex)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "ConsFile":     # constraint file
                try:
                    (param, ConsFile) = line.split()
                    parvar.ConsFile = str(ConsFile)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Region":  # Whole image? or region of image
                try:
                    (param, Region) = line.split()
                    Region = int(Region)
                    parvar.Region = bool(Region)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Boundary":     # Boundary region if Region
                try:
                    (param, Bxmin, Bymin, Bxmax, Bymax) = line.split()
                    parvar.Bxmin = int(Bxmin)
                    parvar.Bxmax = int(Bxmax)
                    parvar.Bymin = int(Bymin)
                    parvar.Bymax = int(Bymax)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Split":  # Split image
                try:
                    (param, Split) = line.split()
                    parvar.Split = int(Split)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

#            if params[0] == "AutoSatRegion":  # Automatic size computation of saturaded regions
#                try:
#                    (param, AutoSatRegion) = line.split()
#                    AutoSatRegion = int(AutoSatRegion)
#                    parvar.AutoSatRegion = bool(AutoSatRegion)
#                except:
#                    print("Unexpected error at reading param file:",
#                          sys.exc_info()[0])
#                    raise

            if params[0] == "SatRegionScale":  # Scale factor for saturaded regions
                try:
                    (param, SatRegionScale) = line.split()
                    parvar.SatRegionScale = float(SatRegionScale)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Ds9SatReg":     # user input ds9 saturation file
                try:
                    (param, Ds9SatReg) = line.split()
                    parvar.Ds9SatReg = str(Ds9SatReg)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "FileOut":     # Preposition name for output files
                try:
                    (param, FileOut) = line.split()
                    parvar.FileOut = str(FileOut)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "ColPar":     # Sextractor catalog
                try:
                    (param, ColPar) = line.split()
                    parvar.ColPar = str(ColPar)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "HeadFlag":
                try:
                    (param, HeadFlag) = line.split()
                    HeadFlag = int(HeadFlag)
                    parvar.HeadFlag = bool(HeadFlag)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise





            if params[0] == "SegFile":     # Preposition name for output segmentation file
                try:
                    (param, SegFile) = line.split()
                    parvar.SegFile = str(SegFile)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SkyFile":     # Preposition name for output segmentation fits for sky
                try:
                    (param, SkyFile) = line.split()
                    parvar.SkyFile = str(SkyFile)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            # Preposition name for output  files containing pixels for every
            # object
            if params[0] == "PixPrefix":
                try:
                    (param, PixPrefix) = line.split()
                    parvar.PixPrefix = str(PixPrefix)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Ds9OutName":    # preposition name for ds9 output
                try:
                    (param, Ds9OutName) = line.split()
                    parvar.Ds9OutName = str(Ds9OutName)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Ds9OutNum":    # Maximum number of objects for ds9
                try:
                    (param, Ds9OutNum) = line.split()
                    parvar.Ds9OutNum = int(Ds9OutNum)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Ds9FitReg":    # Ds9 region output file
                try:
                    (param, Ds9FitReg) = line.split()
                    parvar.Ds9FitReg = str(Ds9FitReg)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Contrast":    # Constrast value ds9
                try:
                    (param, Contrast) = line.split()
                    parvar.Contrast = float(Contrast)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Bias":    # Bias value
                try:
                    (param, Bias) = line.split()
                    parvar.Bias = float(Bias)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "BoxOut":    # Ds9 box region file
                try:
                    (param, BoxOut) = line.split()
                    parvar.BoxOut = str(BoxOut)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "BoxSkyOut":    # Ds9 box sky region file

                try:
                    (param, BoxSkyOut) = line.split()
                    parvar.BoxSkyOut = str(BoxSkyOut)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SexSort":    # Ds9 box sky region file

                try:
                    (param, SexSort) = line.split()
                    parvar.SexSort = str(SexSort)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "SexArSort":    # Ds9 box sky region file
                try:
                    (param, SexArSort) = line.split()
                    parvar.SexArSort = str(SexArSort)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "ErrPer":    # Percentage error limit
                try:
                    (param, ErrPer) = line.split()
                    parvar.ErrPer = float(ErrPer)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Erase":    # erase unnecesary files
                try:
                    (param, Erase) = line.split()
                    Erase = int(Erase)
                    parvar.Erase = bool(Erase)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Nice":    # nice command
                try:
                    (param, Nice) = line.split()
                    Nice = int(Nice)
                    parvar.Nice = bool(Nice)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Overwrite":    # Overwrite files?
                try:
                    (param, Overwrite) = line.split()
                    Overwrite = int(Overwrite)
                    parvar.Overwrite = bool(Overwrite)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise

            if params[0] == "Execute":    # How to Execute (check manual)
                try:
                    (param, Execute) = line.split()
                    parvar.Execute = int(Execute)
                except:
                    print("Unexpected error at reading param file:",
                          sys.exc_info()[0])
                    raise



    return True


def joinsexcat(maincat, secondcat, output, KronScale):
    "merges two Sextractor catalogs"

    f_out = open(output, "w")

    N, Alpha, Delta, X, Y, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg = np.genfromtxt(
        maincat, delimiter="", skip_header=15, unpack=True)

    N = N.astype(int)
    Flg = Flg.astype(int)


    AR = 1 - E
    RKron = KronScale * Ai * Kr

    maskron = RKron <= 0
    if maskron.any():
        RKron[maskron] = 1

    maskar = AR <= 0.005
    if maskar.any():
        AR[maskar] = 0.005

    for idx, item in enumerate(N):

        line = "{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[
                                                                                                        idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

        f_out.write(line)

    total = len(N)

    NewN = total + 1


# second cat
    N2, Alpha2, Delta2, X2, Y2, Mg2, Kr2, Fluxr2, Isoa2, Ai2, E2, Theta2, Bkgd2, Idx2, Flg2 = np.genfromtxt(
        secondcat, delimiter="", skip_header=15, unpack=True)

    N2 = N2.astype(int)
    Flg2 = Flg2.astype(int)


    for idx2, item2 in enumerate(N2):

        flag = False
        for idx, item in enumerate(N):

            flag = CheckKron(X2[idx2], Y2[idx2], X[idx], Y[
                             idx], RKron[idx], Theta[idx], AR[idx])

            if flag:   # boolean value
                break

        if not flag:

            line = "{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(NewN, Alpha2[idx2], Delta2[idx], X2[
                                                                                                            idx2], Y2[idx2], Mg2[idx2], Kr2[idx2], Fluxr2[idx2], Isoa2[idx2], Ai2[idx2], E2[idx2], Theta2[idx2], Bkgd2[idx2], Idx2[idx2], Flg2[idx2])

            f_out.write(line)

            NewN += 1

    f_out.close()





def ds9satbox(satfileout, output, satscale, satoffset):
    "Creates a file for ds9 which selects bad saturated regions"

    scaleflag = 1
    offsetflag = 1
    regfileflag = 1
    magflag = 1
    clasflag = 1

    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag
    checkflag = False
    regflag = 0  # flag for saturaded regions

    f_out = open(satfileout, "w")

    N, Alpha, Delta, X, Y, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg = np.genfromtxt(
        output, delimiter="", unpack=True, dtype="i8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,i8")

    N = N.astype(int)
    Flg = Flg.astype(int)



    for idx, item in enumerate(N):

        bi = Ai[idx] * (1 - E[idx])

        Theta[idx] = Theta[idx] * np.pi / 180  # rads!!!

        Rkronx = satscale * 2 * Ai[idx] * Kr[idx] + satoffset
        Rkrony = satscale * 2 * bi * Kr[idx] + satoffset

        if Rkronx == 0:
            Rkronx = 1

        if Rkrony == 0:
            Rkrony = 1

        # check if object has saturated regions
        checkflag = CheckFlag(Flg[idx], flagsat, maxflag)

        if (checkflag):

            line = "box({0},{1},{2},{3},0) # color=red move=0 \n".format(
                X[idx], Y[idx], Rkronx, Rkrony)
            f_out.write(line)

            line2 = "point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(
                X[idx], Y[idx], N[idx])
            f_out.write(line2)

    f_out.close()


def ds9kron(sexfile, sexfile2, regfile):
    "Put flags on objects which are inside saturated regions"

    f_out = open(sexfile2, "w")

    scale = 1
    offset = 0

    checkflag=False

    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

    N, Alpha, Delta, X, Y, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg = np.genfromtxt(
        sexfile, delimiter="", unpack=True)

    N = N.astype(int)
    Flg = Flg.astype(int)


    for idx, item in enumerate(N):

        Rkron = scale * Ai[idx] * Kr[idx] + offset

        if Rkron == 0:

            Rkron = 1

        q = (1 - E)
        bim = q * Rkron


        checkflag = check.CheckFlag(Flg[idx], flagsat)

        regflag = check.CheckSatReg(X[idx], Y[idx], regfile,
                              Rkron, Theta[idx], E[idx])

        if (checkflag == False) and (regflag == True):

            Flg[idx] = Flg[idx] + 4

            line = "{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[
                                                                                                            idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

            f_out.write(line)

        else:

            line = "{0:.0f} {1} {2} {3} {4} {5} {6} {7} {8:.0f} {9} {10} {11} {12} {13} {14:.0f} \n".format(N[idx], Alpha[idx], Delta[
                                                                                                            idx], X[idx], Y[idx], Mg[idx], Kr[idx], Fluxr[idx], Isoa[idx], Ai[idx], E[idx], Theta[idx], Bkgd[idx], Idx[idx], Flg[idx])

            f_out.write(line)

    f_out.close()



def CatSort(parvar):
    # k Check
    "sort the sextractor catalog by magnitude,"
    "get sizes for objects and write it in a new file"

# The sextractor catalog must contain the following parameters:
#   1 NUMBER                 Running object number
#   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
#   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
#   4 X_IMAGE                Object position along x                                    [pixel]
#   5 Y_IMAGE                Object position along y                                    [pixel]
#   6 MAG_APER               Fixed aperture magnitude vector                            [mag]
#   7 KRON_RADIUS            Kron apertures in units of A or B
#   8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
#   9 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
#  10 A_IMAGE                Profile RMS along major axis                               [pixel]
#  11 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
#  12 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
#  13 BACKGROUND             Background at centroid position                            [count]
#  14 CLASS_STAR             S/G classifier output
#  15 FLAGS                  Extraction flags

    print("Sorting and getting sizes for objects \n")

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        parvar.SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    Rkron = parvar.KronScale * ai * kr

    Rwsky = parvar.SkyScale * ai * kr + parvar.Offset + parvar.SkyWidth

#   considering to use only  KronScale instead of SkyScale:
#    Rwsky = parvar.KronScale * ai * kr + parvar.Offset + parvar.SkyWidth

    Bim = (1 - e) * Rkron

    Area = np.pi * Rkron * Bim*(-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, parvar.NCol, parvar.NRow)

    (sxsmin, sxsmax, sysmin, sysmax) = GetSize(
        xx, yy, Rwsky, theta, e, parvar.NCol, parvar.NRow)

    f_out = open(parvar.SexSort, "w")

    index = mg.argsort()

    for i in index:

        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n".format(n[i], alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i], theta[i], bkgd[i], idx[i], flg[i], np.int(
            np.round(sxmin[i])), np.int(np.round(sxmax[i])), np.int(np.round(symin[i])), np.int(np.round(symax[i])), np.int(np.round(sxsmin[i])), np.int(np.round(sxsmax[i])), np.int(np.round(sysmin[i])), np.int(np.round(sysmax[i])))
        f_out.write(line)

    f_out.close()

    return len(n)


def CatArSort(parvar):
    # k Check

    # sort the sextractor
    # catalog by magnitude,
    # get sizes for objects
    # and write it in a new file

    # The sextractor catalog must contain the following parameters:
    #   1 NUMBER                 Running object number
    #   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
    #   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
    #   4 X_IMAGE                Object position along x                                    [pixel]
    #   5 Y_IMAGE                Object position along y                                    [pixel]
    #   6 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   7 KRON_RADIUS            Kron apertures in units of A or B
    #   8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
    #   9 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
    #  10 A_IMAGE                Profile RMS along major axis                               [pixel]
    #  11 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
    #  12 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
    #  13 BACKGROUND             Background at centroid position                            [count]
    #  14 CLASS_STAR             S/G classifier output
    #  15 FLAGS                  Extraction flags


    print("Sorting and getting sizes for objects \n")

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        parvar.SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    Rkron = parvar.KronScale * ai * kr

    Rwsky = parvar.SkyScale * ai * kr + parvar.Offset + parvar.SkyWidth


#   considering to use only  KronScale instead of SkyScale
#    Rwsky = parvar.KronScale * ai * kr + parvar.Offset + parvar.SkyWidth

    Bim = (1 - e) * Rkron

    Area = np.pi * Rkron * Bim *(-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, parvar.NCol, parvar.NRow)

    (sxsmin, sxsmax, sysmin, sysmax) = GetSize(
        xx, yy, Rwsky, theta, e, parvar.NCol, parvar.NRow)


    f_out = open(parvar.SexArSort, "w")

    index = Area.argsort()
    for i in index:

        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(n[i], alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i], theta[i], bkgd[i], idx[i], flg[i], np.int(
            np.round(sxmin[i])), np.int(np.round(sxmax[i])), np.int(np.round(symin[i])), np.int(np.round(symax[i])), np.int(np.round(sxsmin[i])), np.int(np.round(sxsmax[i])), np.int(np.round(sysmin[i])), np.int(np.round(sysmax[i])))

        f_out.write(line)

    f_out.close()

    return len(n)



def GetSize(x, y, R, theta, ell, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    # k Check
    q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

# getting size

    xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        if isinstance(xmin,np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol

    if mask.any():
        if isinstance(xmax,np.ndarray):
            xmax[mask] = ncol
        else:
            xmax = ncol

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin,np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax,np.ndarray):
            ymax[mask] = nrow
        else:
            ymax = nrow


    return (xmin, xmax, ymin, ymax)

###############
## not used:
###############
def MakeSatDs9(parvar):
    "Make a box DS9 region which contains saturated regions"
# k Check

    checkflag = False
    flagsat = 4
    maxflag = 128

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        parvar.SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    f_out = open(parvar.Ds9SatReg, "w")

    for idx, val in enumerate(n):

        # check if object doesn't has saturated regions
        checkflag = CheckFlag(flg[idx], flagsat, maxflag)

        if (checkflag == True):

            Rsat = parvar.SatRegionScale * ai * kr + 1
            (xmin, xmax, ymin, ymax) = GetSize(
                xx, yy, Rsat, theta, e, parvar.NCol, parvar.NRow)

            sizex = xmax - xmin
            sizey = ymax - ymin

            line = "box({},{},{},{},0) \n".format(
                xx[idx], yy[idx], sizex[idx], sizey[idx])

            f_out.write(line)

    f_out.close()

    return True


def BoxDs9(parvar):
    "Creates a file containing the size for every object in region ds9 file"
# k Check

    checkflag = False
    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, xlo, xhi, ylo, yhi, xslo, xshi, yslo, yshi = np.genfromtxt(
        parvar.SexSort, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    f1_out = open(parvar.BoxOut, "w")
    f2_out = open(parvar.BoxSkyOut, "w")

    xsize = xhi - xlo
    ysize = yhi - ylo

    xssize = xshi - xslo
    yssize = yshi - yslo

    for idx, val in enumerate(n):

        # check if object doesn't has saturated regions
        checkflag = check.CheckFlag(flg[idx], flagsat)

        if (checkflag == False):

            line1 = "box({},{},{},{},0) \n".format(
                xx[idx] + 1, yy[idx] + 1, xsize[idx] + 1, ysize[idx] + 1)
            line2 = "point({},{}) \n".format(xx[idx] + 1, yy[idx] + 1)

            line3 = "box({},{},{},{},0) \n".format(
                xx[idx] + 1, yy[idx] + 1, xssize[idx] + 1, yssize[idx] + 1)
            line4 = "point({},{}) \n".format(xx[idx] + 1, yy[idx] + 1)

            f1_out.write(line1)
            f1_out.write(line2)
            f2_out.write(line3)
            f2_out.write(line4)

    f1_out.close()
    f2_out.close()

    return True


def FindNeighbors():  # modificar para adaptarla al codigo DGCG
    # k Check

    #  This subroutine find neighbors for every galaxy
    #  note for myself: make a tree code of this in the future

    # abrir el archivo en variable global
        #borraer total=Num+1  # num era #num le quite el  #

    KronScale = 1  # debe quitarse

    overflag = 0

    checkflag = False
    checkflag2 = False
    flagsat = 4
    maxflag = 128

    catfile = "skysort"

    N, Alpha, Delta, XPos, YPos, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg, Xlo, Xhi, Ylo, Yhi, Xslo, Xshi, Yslo, Yshi = np.genfromtxt(
        catfile, delimiter="", unpack=True)  # sorteado

    N = N.astype(int)
    Flg = Flg.astype(int)


    Angle = Theta - 90
    AR = 1 - E  # agregado
    RKron = KronScale * Ai * Kr

    for idx, vali in enumerate(N):

        for jdx, valj in enumerate(N):

            checkflag = CheckFlag(Flg[idx], flagsat, maxflag)
            checkflag2 = CheckFlag(Flg[jdx], flagsat, maxflag)

            if (checkflag == False and checkflag2 == False):

                if (idx != jdx):

                    overflag = CheckOverlap(XPos[idx], YPos[idx], RKron[idx], Angle[idx], AR[
                                            idx], XPos[jdx], YPos[jdx], RKron[jdx], Angle[jdx], AR[jdx])

                    if (overflag == True):

                        print ("Elipse {} se traslapa con {}".format(N[idx], N[jdx]))

    return True


def FindNeighbors2(parvar,obj,n):
    # k Check

    #  This subroutine find neighbors for every galaxy
    #  note for myself: make a tree code of this in the future

    # Hace lo mismo que FindNeighbors pero para un solo objeto
    # abrir el archivo en variable global

    parvar.KronScale = 1  # debe quitarse

    overflag = 0

    flagcheck = 0
    flagcheck2 = 0
    flagsat = 4
#    maxflag = 128

    #catfile = "skysort"

    #N, Alpha, Delta, XPos, YPos, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg, Xlo, Xhi, Ylo, Yhi, Xslo, Xshi, Yslo, Yshi = np.genfromtxt(
    #    catfile, delimiter="", unpack=True)  # sorteado

#    N = N.astype(int)
#    Flg = Flg.astype(int)

#    mask = obj.Num == n

#    idx = np.where(mask)[0][0]
    idx=n


#    Angle = Theta - 90
#    AR = 1 - E  # agregado
#    RKron = KronScale * Ai * Kr

    flagcheck = check.CheckFlag(obj.Flag[idx], flagsat)
    A = np.empty((0))

    for jdx, vali in enumerate(obj.Num):

        flagcheck2 = check.CheckFlag(obj.Flag[jdx], flagsat)

        if (flagcheck == False and flagcheck2 == False):

            if (idx != jdx):

                overflag = check.CheckOverlap(obj.XPos[idx], obj.YPos[idx], obj.RKron[idx], obj.Angle[idx], obj.AR[
                                        idx], obj.XPos[jdx], obj.YPos[jdx], obj.RKron[jdx], obj.Angle[jdx], obj.AR[jdx])


                if (overflag == True):
                    #A = np.append(A, obj.Num[jdx])
                    A = np.append(A, jdx)



    A = A.astype(int)
    return A


def UpdateSatFlags(parvar):
    "update flags of saturated regions"

    regflag = 0
    checkflag = False
    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

#    fileflag = 1

    changeflag = False

    runcmd = "cp {} cattemp".format(parvar.SexSort)
    errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errcp)

    print("updating flags for saturaded ds9 box regions \n")

    Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, E, Theta, Background, Class, Flag, XMin, XMax, YMin, YMax, XSMin, XSMax, YSMin, YSMax = np.genfromtxt(
        "cattemp", delimiter="", unpack=True)  # sorted

    Num = Num.astype(int)
    Flag = Flag.astype(int)

    Angle = Theta - 90
    AR = 1 - E
#    RSky = parvar.SkyScale * AIm * Kron + parvar.Offset + parvar.SkyWidth
#    RKron = parvar.KronScale * AIm * Kron
    RSky =  AIm * Kron + parvar.Offset + parvar.SkyWidth
    RKron =  AIm * Kron


    Sky = Background

    masky = RSky <= 0
    if masky.any():
        RSky[masky] = 1

    maskron = RKron <= 0
    if maskron.any():
        RKron[maskron] = 1

    for idx, item in enumerate(Num):

        regflag = check.CheckSatReg(XPos[idx], YPos[idx], parvar.Ds9SatReg, RKron[idx], Theta[idx], E[idx])
        # check if object doesn't has saturated regions
        checkflag = check.CheckFlag(Flag[idx], flagsat)

        if(checkflag == False and regflag == True):

            Flag[idx] += flagsat
            changeflag = True

    os.remove("cattemp")

# writing new catalog if some flag had changed

    if (changeflag == True):

        f_out = open(parvar.SexSort, "w")

        for idx, item in enumerate(Num):

            line = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22} \n".format(Num[idx], RA[idx], Dec[idx], XPos[idx], YPos[idx], Mag[idx], Kron[
                idx], FluxRad[idx], IsoArea[idx], AIm[idx], E[idx], Theta[idx], Background[idx], Class[idx], Flag[idx], XMin[idx], XMax[idx], YMin[idx], YMax[idx], XSMin[idx], XSMax[idx], YSMin[idx], YSMax[idx])
            f_out.write(line)


        f_out.close()

    return True




######

#not checked
# moved to output.py
#def ReadConstraints(file, numobj):
#    "Read constraints file and get constraints for parameters"

# Read Constraints File
# and get contraints for
# parameters


#    flag = 1

# Constraints file

#    magmin = -99
#    magmax = 99

#    remin = -10000
#    remax = 10000

#    nmin = -100
#    nmax = 100

#    with open(file) as INCONS:
#        try:

            # All lines including the blank ones
#            lines = (line.rstrip() for line in INCONS)
#            lines = (line.split('#', 1)[0]
#                     for line in lines)  # remove comments
            # remove lines containing only comments
#            lines = (line.rstrip() for line in lines)
#            lines = (line for line in lines if line)  # Non-blank lines

#            lines = list(lines)
#            index = 0
#            while index < len(lines):

#                line = lines[index]
#                (params) = line.split()

#                if len(params) > 3:

#                    if params[3] == "to":     # input image
#                        if params[0] == numobj:

#                            if(params[1] == "mag" or params[1] == 3):

#                                magmin = params[2]
#                                magmax = params[4]

#                            if(params[1] == "re" or params[1] == 4 or params[1] == "rs"):

#                                remin = params[2]
#                                remax = params[4]

 #                           if(params[1] == "n" or params[1] == 5):

#                                nmin = params[2]
#                                nmax = params[4]

#                index += 1

#        except:
#            print("Can't open constraints file or it was not provided \n")

#    return (magmin, magmax, remin, remax, nmin, nmax)
