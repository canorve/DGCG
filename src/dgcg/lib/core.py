
import os.path
import os
from astropy.io import fits
import numpy as np
import sys
import scipy
import scipy.special
import subprocess as sp


from dgcg.lib import image
from dgcg.lib import check
from dgcg.lib import catfil



######################### Classes ####################

class ParamFile:
    """ Object class stores variables for every of the image"""

    def __init__(self):
        """ Create a new object """

        self.Img = "image.fits"
        self.SexCat = "sex.cat"
        self.SigImg = "none"
        self.PsfDir = "psfdir"
        self.MagZpt = 21.471
        self.PlateScale = 1.0
        self.FitFunc = "BD"

        self.GalClas = 0.5
        self.ConvBox = 60
        self.FitBox = 6
        self.MagDiff = 5
        self.MagCut = 3
        self.KronScale = 1.5
        self.SkyScale = 3
        self.Offset = 20
        self.SkyWidth = 20
        self.NSer = 1.5
        self.MaxFit = 10
#        self.MagMin = 0
        self.MagMax = 19
        self.FlagSex = 4
        self.ConsFile = "constraints"
        self.Region = False
        self.Bxmin = 0
        self.Bymin = 0
        self.Bxmax = 2000
        self.Bymax = 2000
        self.Split = 5
#        self.AutoSatRegion = 0
        self.SatRegionScale = 2
        self.Ds9SatReg = "sat.reg"
        self.SegFile = "seg.fits"
        self.SkyFile = "sky.fits"
#        self.PixFile = "pix.fits"  # Old
        self.PixPrefix = "pixels"
        self.Ds9OutName = "images"
        self.Ds9OutNum = 40
        self.Ds9FitReg = "fit.reg"
        self.BoxOut = "box.reg"
        self.BoxSkyOut = "boxsky.reg"
        self.SexSort = "sexsort.cat"
        self.SexArSort = "sexarsort.cat"
        self.Erase = False
        self.Nice = False
        self.Overwrite = True
        self.Execute = 1
        self.Contrast = 3
        self.Bias = 0.7
        self.SimulFit = True
        self.ErrPer = 0.3

# other variables

        self.Rdnoise = 5.2
        self.ExpTime = 1
        self.Gain    = 7
        self.NCol    = 2000
        self.NRow    = 2000


        self.TempDir = "tempfits"
        self.PixDir = "objpixels"
        self.MaskPixDir = "maskobjpixels"
        self.MaskDir = "Masks"
        self.SkyDir = "Sky"
        self.OutputDir = "Outputs"
        self.InputDir = "Inputs"
        self.OffsetPos = "offsets"  # offset position file
        self.Crashes = "dgcg.crashes"
        self.Fitted = "dgcg.fitted"
        self.SkyCrashes = "sky.crashes"
        self.SkyFitted = "sky.fitted"
        self.RunDir = "RunGALFIT"
        self.ReRunDir = "ReRunGALFIT"



        self.NCol = 2000
        self.NRow = 2000
        self.Buffer = 200
        self.Failures = 0
        self.Success = 0
        self.Total = 0
        self.TotPix = 30  # enough pixels to compute sky?
        self.MinPix = 10  # minimum pixels for object
        self.MaxTimes = 3
        self.ExpTime = 1
        self.Gain = 1
        self.Rdnoise = 0
        self.TimesRe = 3

        self.ListObjs = "ListObjs"

# outputs files

        self.LogFile =  "fits.log"
        self.BtFile =  "fits.bt"
        self.FitsFile = "fits.fits"
        self.ColPar = "selcol.param"
        self.ColFile = "columns.fits"
        self.SexOut = "sexout.txt"
        self.TableFits = "table.fits"
        self.FileOut = "fits.txt"
        self.HeadFlag=True

#  EllipSectGalfit       

        self.ranx=0.7  ## is this Ok?
        self.rany=1  ## is this Ok?
        self.flagsub=False
        self.dpi=100
        self.flagout=False
        self.flagpix=True






class Object:
    """ Object class stores information for every object """

    def __init__(self):
        """ Create a new object """

        self.Num = 0
        self.RA = 0
        self.Dec = 0
        self.XPos = 0
        self.YPos = 0
        self.Mag  = 99
        self.Kron = 0
        self.FluxRad =0
        self.IsoArea =0
        self.AIm  =0
        self.E =  0
        self.Theta =0
        self.Background = 0.0
        self.Class =  0
        self.Flag = 99
        self.XMin  =0
        self.XMax =0
        self.YMin =0
        self.YMax =0
        self.XSMin = 0
        self.XSMax =0
        self.YSMin =0
        self.YSMax =0

        self.Angle = 0

        self.AR = 1

        self.RKron = 0
        self.Sky = 0.0

        self.Sersic = 0

        self.RSky =0
        self.RKron =0

        self.SkyFlag = False

        self.Neighbors = 0

        self.IX = 0
        self.IY = 0


        self.XBuffer =0
        self.YBuffer =0


## OFFX and OFFY help to translate from the big image to the tile image:
        self.OFFX = 0
        self.OFFY = 0

## same as above, but it translate from big image to GALFIT output image:
        self.gOFFX = 0
        self.gOFFY = 0

#####     positions of the galfit output box:
        self.gXMIN = 0
        self.gXMAX = 0
        self.gYMIN = 0
        self.gYMAX = 0
#####

##     These variables store the GALFIT outputs:

        self.ra = "none"
        self.dec = "none"

        self.InputImage = "none"
        self.OutIm = "none"
        self.Objx = "none"
        self.PPPnum = 0
        self.Restart = "none"
        self.Cluster = "none"
        self.PosXSer = 0.0
        self.ErPosXSer = 0.0
        self.PosYSer = 0.0
        self.ErPosYSer = 0.0
        self.MagSer = 99
        self.ErMagSer = 99
        self.ReSer = 0.0
        self.ErReSer = 0.0
        self.NSer = 0.0
        self.ErNSer = 0.0
        self.AxisSer = 0.0
        self.ErAxisSer = 0.0
        self.PaSer = 0.0
        self.ErPaSer = 0.0
        self.KSer = 0.0
        self.MeanMeSer = 99
        self.ErMeanMeSer = 99
        self.MeSer = 99
        self.ErMeSer = 99
        self.PosXExp = 0.0
        self.ErPosXExp = 0.0
        self.PosYExp = 0.0
        self.ErPosYExp = 0.0
        self.MagExp = 99
        self.ErMagExp = 99
        self.RsExp = 0.0
        self.ErRsExp = 0.0
        self.AxisExp = 0.0
        self.ErAxisExp = 0.0
        self.PaExp = 0.0
        self.ErPaExp = 0.0
        self.MeanMeExp = 99
        self.ErMeanMeExp = 99
        self.MeExp = 99
        self.ErMeExp = 99
        self.MsExp = 99
        self.ErMsExp = 99
        self.Sky = 0
        self.ErSky = 0.0
        self.MagTotal = 0.0
        self.ErMagTotal = 0.0
        self.BulgeTotal = 0.0
        self.ErBt = 0.0
        self.ChiNu = 0.0
        self.Tidal = 0.0
        self.ObjChiNu = 0.0
        self.Bump = 0.0
#        self.MeanMesky = 0.0
        self.SNR = 0.0
        self.NDof = 0.0
        self.FitFlag = 0

        self.ReFit = False


#### End of Classes

def RunSky(parvar,obj):
    "Run GALFIT to compute sky for every object"

    fskycrash = open(parvar.SkyCrashes, "w")

    checkflag=False
    flagsat=4
#    maxflag=128

    rangeflag=False

#    Tot = len(Obj.Num)

    for idx, item in enumerate(obj.Num):

        print ("Computing sky value for object {}  \n".format(obj.Num[idx]))

        flagpix=False
        count=0
        countsize=0


        xlo = obj.XSMin[idx]
        xhi = obj.XSMax[idx]

        ylo = obj.YSMin[idx]
        yhi = obj.YSMax[idx]

        xsize = xhi - xlo
        ysize = yhi - ylo


        if parvar.SkyFile.find(".") != -1:
            (TSKY, trash) = parvar.SkyFile.split(".")
        else:
            TSKY=parvar.SkyFile

        if parvar.Img.find(".") != -1:
            (TNAM, trash) = parvar.Img.split(".")
        else:
            TNAM=parvar.Img

        if parvar.SigImg.find(".") != -1:
            (TSIG, trash) = parvar.SigImg.split(".")
        else:
            TSIG=parvar.SigImg

#        (TSKY, trash) = SkyFile.split(".")
#        (TNAM, trash) = Img.split(".")
#        (TSIG, trash) = SigImg.split(".")


        skyname  = TSKY  + "-" + obj.IX[idx].astype(str) + "-" + obj.IY[idx].astype(str) + ".fits"
        imgname  = TNAM  + "-" + obj.IX[idx].astype(str) + "-" + obj.IY[idx].astype(str) + ".fits"
        rmsname  = TSIG  + "-" + obj.IX[idx].astype(str) + "-" + obj.IY[idx].astype(str) + ".fits"

        lindir="{}/{}".format(parvar.TempDir,skyname)

        (ncol,nrow)=image.GetAxis(lindir)


        xfit = obj.XPos[idx] - obj.OFFX[idx]
        yfit = obj.YPos[idx] - obj.OFFY[idx]


        xlo = np.int (xfit - xsize/2)
        if (xlo <= 0):
            xlo = 1

        xhi = np.int (xfit + xsize/2)
        if( xhi > ncol):
            xhi=ncol

        ylo = np.int (yfit - ysize/2)
        if (ylo <= 0):
            ylo = 1

        yhi = np.int (yfit + ysize/2)
        if(yhi > nrow):
            yhi = nrow

## checking for enough pixels

        while (countsize <= parvar.MaxTimes):

            fil="{}/{}".format(parvar.TempDir,skyname)
            count= image.CountPix(fil,xlo,xhi,ylo,yhi,0)

            if (count >= parvar.TotPix):   # change this in the future for a an appropiate number of pixels according to object size

                flagpix = True
                countsize=parvar.MaxTimes + 1

            if (flagpix == True):
                checkflag   = check.CheckFlag(obj.Flag[idx],flagsat)

                if (checkflag == False):
#                    if (obj.Mag[idx] >= parvar.MagMin and obj.Mag[idx] <= parvar.MagMax):
                    if (obj.Mag[idx] <= parvar.MagMax):
                        if (obj.Class[idx] <= parvar.GalClas):
                            if (parvar.Region == 0 or (parvar.Region == 1 and obj.XPos[idx] >= parvar.Bxmin and obj.XPos[idx] <= parvar.Bxmax and obj.YPos[idx] >= parvar.Bymin and obj.YPos[idx] <= parvar.Bymax)):

                                obj.Sky[idx] = GetSky(parvar,obj,idx,skyname,imgname,rmsname,xlo,xhi,ylo,yhi)

                            else:
                                obj.Sky[idx] = obj.Background[idx]
                                obj.SkyFlag[idx] = False
                                print ("Skipping object {}, it is outside region ({},{},{},{}) \n".format(obj.Num[idx],parvar.Bxmin,parvar.Bxmax,parvar.Bymin,parvar.Bymax))
                        else:
                            obj.Sky[idx]= obj.Background[idx]
                            obj.SkyFlag[idx] = False
                            print ("Skipping object {}, Maybe it is not a galaxy ({}) \n".format(obj.Num[idx],parvar.GalClas))
                    else:
                        obj.Sky[idx]= obj.Background[idx]
                        obj.SkyFlag[idx] = False
                        print ("Skipping object {}, {} is greater than magnitude limit MagMagx ({})  \n".format(obj.Num[idx],obj.Mag[idx],parvar.MagMax))
                        if (obj.Mag[idx] > parvar.MagMax):
                            rangeflag =True
                else:
                    obj.Sky[idx] = obj.Background[idx]
                    obj.SkyFlag[idx] = False
                    print ("Skipping object {}, it has saturated regions  \n".format(obj.Num[idx]))
                    #print ("Flag = ",obj.Flag[idx])

            elif(countsize < parvar.MaxTimes):
                # double the size every time
                xsize = 2 * (xhi - xlo)
                ysize = 2 * (yhi - ylo)

                # new sizes:

                xlo = np.int (xfit - xsize/2)
                xhi = np.int (xfit + xsize/2)

                ylo = np.int (yfit -  ysize/2)
                yhi = np.int (yfit +  ysize/2)


                if (xlo <= 0):
                    xlo = 1

                if (xhi > ncol):
                    xhi = ncol

                if (ylo <= 0):
                    ylo = 1

                if (yhi > nrow):
                    yhi = nrow

                countsize+=1

            else:
                objsky = "Sky-{}".format(obj.Num[idx])
                line = "{} \n".format(objsky)
                fskycrash.write(line)
                print ("Not enough pixels to compute sky for object {}. Sextractor value will be used. \n".format(obj.Num[idx]))
                obj.SkyFlag[idx] = False
                break

        if(rangeflag == True):
            break

    fskycrash.close()

    return True



def GetSky(parvar,obji,i, skyname, imgname, rmsname, xlo, xhi, ylo, yhi):
    "Run GALFIT to compute sky value"

    # k Check

    # Run galfit to compute sky value
    # if the fit crash, it takes the SExtractor sky value
    # This function is called by RunSky

#    catfile = "sexsort.cat"

    #N, Alpha, Delta, obj.XPos, obj.YPos, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg, Xlo, Xhi, Ylo, Yhi, Xslo, Xshi, Yslo, Yshi = np.genfromtxt(
    #    catfile, delimiter="", unpack=True)  # sorteado

    #N = N.astype(int)
    #Flg = Flg.astype(int)
    Z = 0
    fit = 1
    iden = obji.Num[i]
    back = obji.Background[i]  # le puse Bkgd por el catalogo
    skyval = 0
#    fileflag = 1

    findflag=False

# open file to run on GALFIT

#    objsky = "Sky-{}".format(iden) #oldie
    objsky = "{}/Sky-{}".format(parvar.SkyDir,iden)

    objskyy = "Sky-{}-out.fits".format(iden)  # complete con -out.fits
#    objskyy = "{}/Sky-{}-out.fits".format(parvar.SkyDir,iden)  # complete con -out.fits



    T1 = "{}/{}".format(parvar.TempDir,imgname)
    T2 = "{}/{}".format(parvar.TempDir,rmsname)
    T3 = "{}/{}".format(parvar.TempDir,skyname)

    SKYOUT = open(objsky, "w")

    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    PrintHeader(SKYOUT, T1, objskyy, T2, "NONE", 1, T3, "NONE", xlo, xhi, ylo,
                yhi, parvar.ConvBox, parvar.ConvBox, parvar.MagZpt, parvar.PlateScale, parvar.PlateScale, "regular", 0, 0)

# creates sky component

    PrintSky(SKYOUT, iden, back, Z, fit)
    SKYOUT.close()

#    f = open("skytemp.txt", "w")
#    err = sp.call(["galfit", objsky], stdout=f)  # Run GALFIT
#    f.close()

#####
    runcmd = "galfit {} > skytemp.txt".format(objsky) #run GALFIT
    errno = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)
#####


    findflag=SearchObj(iden)

    if (findflag == False):

        # print SkyCrashed "{} \n".format(objsky)
        skyval = back
        # SkyFlag[i]=False

        # else:

        # print SkyFitted "{} \n".format(objsky)
        # open fit.log and  get sky value
    else:
        skyval=FitLogSky(objsky)

# erasing unnecessary files


    os.remove(objskyy)
    os.remove("skytemp.txt")

    return float(skyval)


def FitLogSky(Objsky):
    "Reads fit.log file and gets sky value"
    # K

    fitlog = "fit.log"

    val = 0

    with open(fitlog) as INCONF:

        # All lines including the blank ones
        lines = (line.rstrip() for line in INCONF)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines
        #it= iter(lines)
        lines = list(lines)
        index = 0
        # print lines
        while index < len(lines):
            # for lines in lines:
            line = lines[index]
            (params) = line.split()

            if params[0] == "Init.":     # input image
                if params[4] == Objsky:

                    #print (line)

                    index += 3  # go forward 3 lines

                    line = lines[index]

                    #print (line)

                    # print  line
                    (left,right) = line.split(",")


                    (trash1, param, trash2, trash3) = right.split()


                    param = float(param)
                    val = param

                    break

            index += 1

    return val


def PrintHeader(hdl, A, B, C, D, E, F, G, xlo, xhi, ylo, yhi, convx, convy, J, platedx, platedy, O, P, S):
    "print GALFIT header in a file"

    # k Check
    # print to filehandle
    # the header for GALFIT

    lineZ = "==================================================================================================\n"
    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}                                   # Input Data image (FITS file)                            \n".format(A)
    lineB = "B) {}                                   # Output data image block                                 \n".format(B)
    lineC = "C) {}                                   # Sigma image name (made from data if blank or \"none\")  \n".format(C)
    lineD = "D) {}                                   # Input PSF image and (optional) diffusion kernel         \n".format(D)
    lineE = "E) {}                                   # PSF fine sampling factor relative to data               \n".format(E)
    lineF = "F) {}                                   # Bad pixel mask (FITS image or ASCII coord list)         \n".format(F)
    lineG = "G) {}                                   # File with parameter constraints (ASCII file)            \n".format(G)
    lineH = "H) {} {} {} {}                          # Image region to fit (xmin xmax ymin ymax)               \n".format(xlo, xhi, ylo, yhi)
    lineI = "I) {} {}                                # Size of the convolution box (x y)                       \n".format(convx, convy)
    lineJ = "J) {}                                   # Magnitude photometric zeropoint                         \n".format(J)
    lineK = "K) {} {}                                # Plate scale (dx dy). [arcsec per pixel]               \n".format(platedx, platedy)
    lineO = "O) {}                                   # Display type (regular, curses, both)                    \n".format(O)
    lineP = "P) {}                                   # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps      \n".format(P)
    lineS = "S) {}                                   # Modify/create objects interactively?                    \n".format(S)
    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS                                                     \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:                                      \n"
    line3 = "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,             \n"
    line4 = "#       ferrer, powsersic, sky, and isophote.                                    \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified:                    \n"
    line7 = "#      C0 (diskyness/boxyness),                                                  \n"
    line8 = "#      Fn (n=integer, Azimuthal Fourier Modes),                                  \n"
    line9 = "#      R0-R10 (PA rotation, for creating spiral structures).                     \n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number                                                               \n"
    line12 = "# column 2:                                                                                 \n"
    line13 = "#          -- Parameter 0:    the allowed functions are: sersic, nuker, expdisk             \n"
    line14 = "#                             edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky    \n"
    line15 = "#          -- Parameter 1-10: value of the initial parameters                               \n"
    line16 = "#          -- Parameter C0:   For diskiness/boxiness                                        \n"
    line17 = "#                             <0 = disky                                                    \n"
    line18 = "#                             >0 = boxy                                                     \n"
    line19 = "#          -- Parameter Z:    Outputting image options, the options are:                    \n"
    line20 = "#                             0 = normal, i.e. subtract final model from the data to create \n"
    line21 = "#                             the residual image                                            \n"
    line22 = "#                             1 = Leave in the model -- do not subtract from the data       \n"
    line23 = "#                                                                                           \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0)                                       \n"
    line25 = "# column 4: comment                                                                         \n"
    line26 = " \n"

    line27 = "==================================================================================================\n"

    hdl.write(lineZ)
    hdl.write(lineX)
    hdl.write(lineA)
    hdl.write(lineB)
    hdl.write(lineC)
    hdl.write(lineD)
    hdl.write(lineE)
    hdl.write(lineF)
    hdl.write(lineG)
    hdl.write(lineH)
    hdl.write(lineI)
    hdl.write(lineJ)
    hdl.write(lineK)
    hdl.write(lineO)
    hdl.write(lineP)
    hdl.write(lineS)
    hdl.write(lineY)

    hdl.write(line0)
    hdl.write(line1)
    hdl.write(line2)
    hdl.write(line3)
    hdl.write(line4)
    hdl.write(line5)
    hdl.write(line6)
    hdl.write(line7)
    hdl.write(line8)
    hdl.write(line9)
    hdl.write(line10)

    hdl.write(line11)
    hdl.write(line12)
    hdl.write(line13)
    hdl.write(line14)
    hdl.write(line15)
    hdl.write(line16)
    hdl.write(line17)
    hdl.write(line18)
    hdl.write(line19)
    hdl.write(line20)
    hdl.write(line21)
    hdl.write(line22)
    hdl.write(line23)
    hdl.write(line24)
    hdl.write(line25)
    hdl.write(line26)
    hdl.write(line27)

    return True


def PrintSky(hdl, ncomp, sky, Z, fit):
    "Print GALFIT sky function to filehandle"

    # k Check

    line00 = "# Object number: {}                                                             \n".format(ncomp)
    line01 = " 0)      sky            #    Object type                                        \n"
    line02 = " 1) {}         {}       # sky background        [ADU counts]                    \n".format(sky, fit)
    line03 = " 2) 0.000      0        # dsky/dx (sky gradient in x)                           \n"
    line04 = " 3) 0.000      0        # dsky/dy (sky gradient in y)                           \n"
    line05 = " Z) {}                  # Skip this model in output image?  (yes=1, no=0)       \n".format(Z)
    line06 = "\n"
    line07 = "================================================================================\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)

    return True


def PrintSersic(hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit):
    "print GALFIT Sersic function to filehandle"
    # k Check

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}                                                             \n".format(
            ncomp)
    line01 = " 0)     sersic               #  Object type                                     \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}            #  position x, y     [pixel]                       \n".format(
        xpos, ypos, fit, fit)
    line03 = " 3) {:.2f}       {}              #  total magnitude                                 \n".format(
        magser, fit)
    line04 = " 4) {:.2f}       {}              #  R_e         [Pixels]                            \n".format(
            reser, fit)
    line05 = " 5) {}       {}              #  Sersic exponent (deVauc=4, expdisk=1)           \n".format(
        nser, fit)
    line06 = " 9) {:.2f}       {}              #  axis ratio (b/a)                                \n".format(
        axratser, fit)
    line07 = "10) {:.2f}       {}              #  position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        angleser, fit)
    line08 = " Z) {}                       #  Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line09 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)
    hdl.write(line09)

    return True

def PrintExp(hdl, ncomp, xpos, ypos, magexp, rsexp, axratexp, angleexp, Z, fit):
    "print GALFIT exponential function to filehandle"

    # k Check

    # print to filehandle
    # a exponential function given
    # by the parameters

    line00 = "# Object number: $ncomp                                                        \n"
    line01 = " 0)     expdisk              # Object type                                     \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}           # position x, y     [pixel]                       \n".format(
        xpos, ypos, fit, fit)
    line03 = " 3) {:.2f}        {}             # total magnitude                                 \n".format(
        magexp, fit)
    line04 = " 4) {:.2f}        {}             #      Rs  [Pixels]                               \n".format(
        rsexp, fit)
    line05 = " 9) {:.2f}        {}             # axis ratio (b/a)                                \n".format(
        axratexp, fit)
    line06 = "10) {:.2f}        {}             # position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        angleexp, fit)
    line07 = " Z) {}                       # Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line08 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)

    return True




def DGCG(parvar,obj,flog,fobjs,fout2,fout3,fout4):
    "Formats Sersic and exponential models in files and run GALFIT"

    checkflag = False
    flagsat = 4
#    maxflag = 128

    ii=0
    for idx, item in enumerate(obj.Num):

        checkflag = check.CheckFlag(obj.Flag[idx], flagsat)

        if (checkflag == False):
#            if (obj.Mag[idx] >= parvar.MagMin and obj.Mag[idx] <= parvar.MagMax):
            if (obj.Mag[idx] <= parvar.MagMax):
                if (obj.Class[idx] <= parvar.GalClas):
                    if (parvar.Region == 0 or (parvar.Region == 1 and obj.XPos[idx] >= parvar.Bxmin and obj.XPos[idx] <= parvar.Bxmax and obj.YPos[idx] >= parvar.Bymin and obj.YPos[idx] <= parvar.Bymax)):

                        print("Making input file for object number: {} \n".format(obj.Num[idx]))
                        line = "{} {} \n".format(obj.Num[idx], idx)
                        fobjs.write(line)

###################################### Creates obj file for galfit
                        objfile = MakeObjFile(parvar,obj,idx,flog,fout2)
######################################

                        objid = obj.Num[idx]
                        outfile = parvar.OutputDir + "/" + "out" + "-" + str(objid)

                        print("Galaxy number {} \n".format(ii))
                        print("Running GALFIT on object number: {} \n".format(obj.Num[idx]))

######################################## RUN GALFIT
                        GALFIT(objfile, outfile, parvar, objid, idx, fout3, fout4)
#########################################

                    else:
                        print("Object {} is rejected because its position ({},{}) is outside of Region ({},{},{},{})\n".format(
                            obj.Num[idx], obj.XPos[idx], obj.YPos[idx], parvar.Bxmin, parvar.Bxmax, parvar.Bymin, parvar.Bymax))
                else:
                    print("Object {} is rejected because maybe it's not a Galaxy {} is greater than range ({}) \n".format(
                        obj.Num[idx], obj.Class[idx],  parvar.GalClas))
            else:
                print("Object {} is rejected due to Magnitude {} is greater than magnitude limit ({}) \n".format(
                    obj.Num[idx], obj.Mag[idx], parvar.MagMax))
                if (obj.Mag[idx] > parvar.MagMax):
                    break     # finishing loop.  catalog is magnitude sorted.
        else:
            print("Object {} is rejected because it has saturated regions \n".format(obj.Num[idx]))

        ii+=1
    return True


def MakeObjFile(parvar,obj,nobj,flog,fout2):
    "makes input param files for GALFIT"

    objid = obj.Num[nobj]
    bkg = obj.Sky[nobj]
    skyflag = obj.SkyFlag[nobj]

    ncomp = 0
    parfile =  "obj" + "-" + str(objid)

# run on directory:
    parfile =  parvar.RunDir + "/" + parfile

#    print("ncomp: ",ncomp)


    fout1 = open(parfile, "w")

#    xsize = obj.XMax[nobj] - obj.XMin[nobj]
#    ysize = obj.YMax[nobj] - obj.YMin[nobj]

#   enlarge fit area

#    xsize = parvar.FitBox * xsize
#    ysize = parvar.FitBox * ysize

##  30 pixels is the minimum area to fit (arbitrary number):

#    if (xsize <= 30):
#        xsize = 30

#    if (ysize <= 30):
#        ysize = 30

    if parvar.SegFile.find(".") != -1:
        (TSEG, trash) = parvar.SegFile.split(".")
    else:
        TSEG=parvar.SegFile

    if parvar.Img.find(".") != -1:
        (TNAM, trash) = parvar.Img.split(".")
    else:
        TNAM=parvar.Img

    if parvar.SigImg.find(".") != -1:
        (TSIG, trash) = parvar.SigImg.split(".")
    else:
        TSIG=parvar.SigImg


    segname = TSEG + "-" + str(obj.IX[nobj]) + "-" + str(obj.IY[nobj]) + ".fits"
    imgname = TNAM + "-" + str(obj.IX[nobj]) + "-" + str(obj.IY[nobj]) + ".fits"
    rmsname = TSIG + "-" + str(obj.IX[nobj]) + "-" + str(obj.IY[nobj]) + ".fits"

    line = "{}/{}".format(parvar.TempDir, imgname)
    (ncol, nrow) = image.GetAxis(line)

#   Figure out which PSF is closest to object center

    psfname = PSFDist(obj.XPos[nobj], obj.YPos[nobj])


#   Construct postage stamp name

    (RAcel, DECcel) = Deg2Celestial(obj.RA[nobj], obj.Dec[nobj])


    outname = Uniqname(RAcel, DECcel, TNAM)


    outname = outname + "-" + str(objid)

    outname =  parvar.OutputDir + "/" + outname

    maskfile = "mask-" + str(objid) + ".fits"

# maskfile lives at:
    maskfile =  parvar.MaskDir + "/" + maskfile


    runcmd = "cp {}/{} {}".format(parvar.TempDir, segname, maskfile)
    errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)

    #  Calculate the (x,y) position of the current object relative to
    #  the tile in which it lives.

    xfit = obj.XPos[nobj] - obj.OFFX[nobj]
    yfit = obj.YPos[nobj] - obj.OFFY[nobj]

    #  Calculate fitting box needed to plug into galfit header:

#    xlo = int(xfit - xsize / 2)
#    if (xlo <= 0):
#        xlo = 1

#    xhi = int(xfit + xsize / 2)
#    if (xhi > ncol):
#        xhi = ncol

#    ylo = int(yfit - ysize / 2)
#    if (ylo <= 0):
#        ylo = 1

#    yhi = int(yfit + ysize / 2)
#    if (yhi > nrow):
#        yhi = nrow


    # Calculate information needed to plug back into the image header
    # at the end of the fit.

    xlo = obj.gXMIN[nobj]  - obj.OFFX[nobj]
    xhi = obj.gXMAX[nobj]  - obj.OFFX[nobj]
    ylo = obj.gYMIN[nobj]  - obj.OFFY[nobj]
    yhi = obj.gYMAX[nobj]  - obj.OFFY[nobj]

#    obj.gXMIN[nobj] = xlo + obj.OFFX[nobj]  # The [xmin:xmax,ymin:ymax] of the box
#    obj.gXMAX[nobj] = xhi + obj.OFFX[nobj]  # relative to the big image from which
#    obj.gYMIN[nobj] = ylo + obj.OFFY[nobj]  # the current tile was extracted.
#    obj.gYMAX[nobj] = yhi + obj.OFFY[nobj]


# xlo, ylo, xhi, yhi correspond to the small image
# gXMIN, gYMIN, gXMAX, gYMAX correspond to the big image

    ####################################
    #  Create GALFIT input file header #
    ####################################


    T1 = "{}/{}".format(parvar.TempDir, imgname)
    T2 = outname + "-out.fits"
    T3 = "{}/{}".format(parvar.TempDir, rmsname)

    PrintHeader(fout1, T1, T2, T3, psfname, 1, maskfile, parvar.ConsFile, xlo, xhi, ylo,
                yhi, parvar.ConvBox, parvar.ConvBox, parvar.MagZpt, parvar.PlateScale, parvar.PlateScale, "regular", 0, 0)

    # create components for the main galaxy:

    ncomp = Components(fout1, obj, nobj, parvar.FitFunc, xfit, yfit, ncomp)
#    print("ncomp: ",ncomp)


    # create components or mask for the neighbours galaxies:

    ncomp = Patch(fout1,flog, maskfile, obj, ncomp, nobj, parvar)

#    print("ncomp: ",ncomp)

    if (skyflag == True):
        PrintSky(fout1, objid, bkg, 0, 0)
#        ncomp += 1  # exclude sky as a component to fit
    elif(skyflag == False):
        PrintSky(fout1, objid, bkg, 0, 1)
#        ncomp += 1  # exclude sky as a component to fit

    print("A total number of {} components for object {} will be fitted \n".format(ncomp, objid))

    line = "A total number of {} components for object {} will be fitted \n".format(ncomp, objid)

    flog.write(line)


    line = "{} {}-out   {} {}  [{}:{},{}:{}] {} {} \n".format(
        objid, outname, obj.OFFX[nobj], obj.OFFY[nobj], obj.gXMIN[nobj], obj.gXMAX[nobj], obj.gYMIN[nobj], obj.gYMAX[nobj], obj.XPos[nobj], obj.YPos[nobj])


    fout2.write(line)
    flog.write(line)  # output was changed to log file

    fout1.close()


    return parfile


def PSFDist(xg, yg):
    "computes the nearest PSF to the object"

# computes the nearest PSF to the object
# Choose the appropriate convolution PSF.  Given a list of PSF names with
# names PSF-x-y.fits, located at position (x, y), figure out for a galaxy
# at position xgal, ygal, which PSF is the closest to this object.


#    psfs = np.genfromtxt("psf.temp", delimiter="", dtype= np.dtype((str, 23)) , unpack=True)
    psfs = np.genfromtxt("psf.temp", delimiter="", dtype= np.dtype(str) , unpack=True)



    neardist = 1.e6

    for idx, item in enumerate(psfs):


        (dum, xp, yp) = psfs[idx].split("-")
        (yp,trash) = yp.split(".fits")

        xp=float(xp)
        yp=float(yp)

        dist = np.sqrt((xp - xg)**2 + (yp - yg)**2)

        if (dist < neardist):

            neardist = dist
            usepsf = psfs[idx]


    return (usepsf)



def Deg2Celestial(ra, dec):
    "returns celestial coordinates from input of ra and dec in decimal degree format"

    # Call w/ ($ra_celestial,$dec_celestial) = deg2celestial($ra,$dec);
    # returns celestial coordinates from input of ra and dec in decimal degree
    # format.  This subroutine returns all significant digits in seconds of
    # RA and Dec; thus simpler version of deg2celestial.sub w/out roundoff.
    # #NOTE: Run this subroutine in a loop, only does one transf. at a time.
    # 6/9/02: fixed bug for -1<Dec<0.

    # subroutine deg2cel.precise  (version 06/09/02)

####################
    # RA portion
####################

    hr = ra / 15.
    hr = str(hr)
    (rah, junk) = hr.split(".")

    stnum = "." + junk
    mins = 60. * float(stnum)
    mins = str(mins)
    (ram, junk) = mins.split(".")

    stnum = "." + junk
    ras = 60. * (float(stnum))
    ras = "{0:.2f}".format(ras)
#    ras = str(ras)


    # always want double digits
    if float(rah) < 10:
        rah = "0" + rah
    if float(ram) < 10:
        ram = "0" + ram
    if float(ras) < 10:
        ras = "0" + ras

# rounding to two digits


    ra_celestial = rah + ":" + ram + ":" + ras

#################
    # Dec portion
#################

    dec=str(dec)
    (degree, junk) = dec.split(".")
    if (float(degree) >= 0.0):
        sign = "p"
        if (degree.find("+") != -1):
            (junk1, deg) = degree.split("+")
        else:
            deg = degree

    if (float(degree) <= 0.0 and float(dec) < 0.0):
        sign = "m"
        if (degree.find("-") != -1):
            (junk1, deg) = degree.split("-")
        else:
            deg = degree


    stnum = "." + junk
    dmin = 60. * (float(stnum))
    dmin = str(dmin)

    (amin, junk) = dmin.split(".")

    stnum = "." + junk
    asec = 60. * (float(stnum))
    asec = "{0:.2f}".format(asec)
#    asec = str(asec)

    # always want double digits
    if (float(deg) < 10):
        deg = "0" + deg
    if (float(amin) < 10):
        amin = "0" + amin
    if (float(asec) < 10):
        asec = "0" + asec

    # put sign back on
    if (sign == "m"):
        dec_celestial = "-" + deg + ":" + amin + ":" + asec
    if (sign == "p"):
        dec_celestial = "+" + deg + ":" + amin + ":" + asec

    return (ra_celestial, dec_celestial)


def Uniqname(ra_celest, dec_celest, tname):
    "returns one name for ra and dec"

    (ra1, ra2, ra3) = ra_celest.split(":")
#    $ra3 = sprintf("%.2f", $ra3);      #rounds sec to 0.01
    if (float(ra3) < 10.):
        ra3 = "0" + ra3

    ra_name = ra1 + ra2 + ra3
    (d1, d2, d3) = dec_celest.split(":")
#    $d3 = sprintf("%.1f", $d3);        #rounds asec to 0.1
    if (float(d3) < 10.):
        d3 = "0" + d3

    if (float(d1) < 0.):
        sign = "m"

    if (float(d1) >= 0.):
        sign = "p"

#    d1 = substr(d1,1,2)         #deletes sign
    d1 = d1[1:3]        # equivalent of substr
    dec_name = sign + d1 + d2 + d3


##
    out_name = tname + '-' + ra_name + '-' + dec_name
###


    return out_name


def Components (filehandle,obj,i,function,xpos,ypos,ncomp):
    "It creates components for bulge and disk from SEXtractror parameters"

## It creates components for bulge and disk
## or sersic functions it redistribuites the
## SEXtractror parameters in each component.

#    my ($nser,$magser,$reser,$axratser,$angleser);
#    my ($magexp,$rsexp,$axratexp,$angleexp,$Z,$fit,$objid);

    Z=0
    fit=1

    objid = obj.Num[i]

    if (function == "sersic"):
        magser   = obj.Mag[i]
        reser    = obj.FluxRad[i]
        nser     = obj.Sersic[i]
        axratser = obj.AR[i]
        angleser = obj.Angle[i]

    elif (function == "BD"):

        #Bulge  initial parameters
        magser   = obj.Mag[i]     + 0.5
        reser    = obj.FluxRad[i] * 0.7
        nser     = obj.Sersic[i]
        axratser = 0.8
        angleser = obj.Angle[i]

        # Disk initial parameters
        magexp    = obj.Mag[i]     + 0.5
        rsexp    =  obj.FluxRad[i]
        axratexp =  obj.AR[i]
        angleexp  = obj.Angle[i]

    ncomp+=1
    PrintSersic(filehandle,objid,xpos,ypos,magser,reser,nser,axratser,angleser,Z,fit)

    if (function == "BD"):
        ncomp+=1
        PrintExp(filehandle,objid,xpos,ypos,magexp,rsexp,axratexp,angleexp,Z,fit)


    return ncomp


def Patch(filehdl,flog, maskfile, obj, ncomp, nobj, parvar):

    "This subroutine make mask or create components for neighbours objects of the main one"

    #  this subroutine make mask or
    #  create components for the neighbours
    #  objects of the main galaxy.

    checkflag = False
    flagsat = 4
#    maxflag = 128

    mainid = obj.Num[nobj]


##### OLD:
#    pixfile = parvar.PixPrefix + "-" + str(mainid)


#################################

#    neigh = obj.Neighbors[nobj]  # change this for FindNeighbors2  ## Old
    neigh = catfil.FindNeighbors2(parvar,obj,nobj)

##################################3


    count = 1

    Z = 0


    xlo = obj.gXMIN[nobj] - obj.OFFX[nobj]
    xhi = obj.gXMAX[nobj] - obj.OFFX[nobj]
    ylo = obj.gYMIN[nobj] - obj.OFFY[nobj]
    yhi = obj.gYMAX[nobj] - obj.OFFY[nobj]



    line="removing main galaxy mask {} ".format(mainid)
    print(line)
    print(maskfile, mainid, xlo, xhi , ylo , yhi)
    print(maskfile, mainid, xlo, obj.gXMIN[nobj] ,obj.OFFX[nobj])

    image.SetZero(maskfile, mainid, xlo, xhi, ylo, yhi)  # remove main object from mask


    if (parvar.SimulFit == True):
        for idx, item in enumerate(neigh):  # loop on neighbors
            flagser = False
            flagexp = False

##################  object
            # this could be substituted by item
            neighn = item #neigh[idx]  # neighbors

            neighid = obj.Num[neighn]

            neighflag = obj.Flag[neighn]

            diff = obj.Mag[neighn] - obj.Mag[nobj]

            magcutdiff = obj.Mag[neighn] - parvar.MagMax

#            pixfile = parvar.PixPrefix + "-" + str(neighid)

            checkflag = check.CheckFlag(neighflag, flagsat)


            if ((diff <= parvar.MagDiff) and (magcutdiff <= parvar.MagCut) and (count <= parvar.MaxFit) and (checkflag == False) and (obj.Class[neighn] <= parvar.GalClas) ):
                # do simultaneous fitting
                count += 1

                line = "Galaxy {} is being simultaneously fitted with main galaxy: {} \n".format(obj.Num[neighn], mainid)
                flog.write(line)


                xlon = obj.gXMIN[neighn] - obj.OFFX[neighn]
                xhin = obj.gXMAX[neighn] - obj.OFFX[neighn]
                ylon = obj.gYMIN[neighn] - obj.OFFY[neighn]
                yhin = obj.gYMAX[neighn] - obj.OFFY[neighn]

                line="removing galaxy mask {} for main galaxy mask {}".format(obj.Num[neighn],mainid)
                print(line)
                print(maskfile, obj.Num[neighn], xlon, xhin , ylon , yhin)
                print(maskfile, obj.Num[neighn], xlon, obj.gXMIN[neighn] ,obj.OFFX[neighn] )
                image.SetZero(maskfile, obj.Num[neighn], xlon, xhin , ylon , yhin)  # remove galaxy which will be simultaneously fitted from mask

                if os.path.isfile("fit.log"):
                    (flagser, xpos, ypos, magser, reser, nser, axratser, angleser) = SearchSersic(neighid)  #checar
                else:
                    flagser=False


                if(flagser == True):   # initial parameters from fit.log

                    if (parvar.FitFunc == "BD"):
                        (flagexp, xposexp, yposexp, magexp, rsexp,axratexp, angleexp) = SearchExp(neighid)  #checar

                    else:
                        flagexp=False

                    if(flagexp == True):   ## BULGE DISK COMPONENT OR ...

#                        PrintSersic(filehdl, neighid, xpos, ypos, magser, reser, nser, axratser, angleser, Z, 1)
                        ncomp += 1
                        PrintSersic(filehdl, neighid, xpos, ypos, magser-0.5, 2*reser, nser, axratser, angleser, Z, 1) #test

#                        PrintExp(filehdl, neighid, xposexp, yposexp,magexp, rsexp, axratexp, angleexp, Z, 1)
#                        ncomp += 1
                        line = "Object {} will be simultaneously fitted with their parameters free \n".format(neighid)
                        flog.write(line)

                    elif(flagexp == False):  ### ... SERSIC COMPONENT

                        ncomp += 1
                        PrintSersic(filehdl, neighid, xpos, ypos, magser,reser, nser, axratser, angleser, Z, 1)
                        line = "Object {} will be simultaneously fitted with their parameters free \n".format(neighid)
                        flog.write(line)

#############################  THIS PART OF PARAMETERS FIXED WAS REMOVED

#                    elif(flagexp == False):
#                        PrintSersic(filehdl, neighid, xpos, ypos, magser,
#                                    reser, nser, axratser, angleser, Z, 0)
#                        ncomp += 1
#                        line = "Object {} will be simultaneouly fitted with their parameters fixed \n".format(
#                            neighid)
#                        flog.write(line)

                else:   ## initial parameters from sextractor
                    xpos = obj.XPos[neighn] - obj.OFFX[nobj]
                    ypos = obj.YPos[neighn] - obj.OFFY[nobj]

                    magser = obj.Mag[neighn]
                    reser = obj.FluxRad[neighn]
                    nser = obj.Sersic[neighn]
                    axratser = obj.AR[neighn]
                    angleser = obj.Angle[neighn]

                    line = "Object {} will be simultaneouly fitted with their parameters free \n".format(
                        neighid)
                    flog.write(line)
                    ncomp += 1
                    PrintSersic(filehdl, neighid, xpos, ypos, magser,reser, nser, axratser, angleser, Z, 1)
            else:
                line = "Galaxy {} is too faint ({} > {}) compared to main galaxy {} or it's not a galaxy \n".format(
                    obj.Num[neighn], obj.Mag[neighn], obj.Mag[nobj], mainid)
                flog.write(line)
#                print("(diff <= parvar.MagDiff) and (magcutdiff <= parvar.MagCut) and (count <= parvar.MaxFit) and (checkflag == False) and (obj.Class[neighn] <= parvar.GalClas)" )
#                line="({} > {}) and ({} > {}) and ({} > {}) and ({} != False) and ({} > {})".format(diff, parvar.MagDiff, magcutdiff, parvar.MagCut, count, parvar.MaxFit, checkflag,obj.Class[neighn],parvar.GalClas)
                print (line)

    return ncomp


def SearchSersic(neighid):
    "Search Sersic component in fit.log file and returns parameters"
    # Search sersic component
    # in fit.log file
    # and returns parameters

    # outputs = sxpos,sypos,smagser,sreser,snser,saxratser,sangleser

#    flag = 0
    errflag = False
#    fitlogflag = False
    findflag = False


# putting some limit to errors

    neighobj = "obj" + "-" + str(neighid)

    # eliminar:

    #consflag = ReadConstraints(constraints,1,consmagmin1,consmagmax1,consremin1,consremax1,consnmin1,consnmax1)
    #consflag = ReadConstraints(constraints,2,consmagmin2,consmagmax2,consremin2,consremax2,consnmin2,consnmax2)

    with open("fit.log") as FITLOG:


        # All lines including the blank ones
        lines = (line.rstrip() for line in FITLOG)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0
        while index < len(lines):

            line = lines[index]
            (params) = line.split()

            if params[0] == "Init.":     # input image

                if params[4] == neighobj:
                    # print line
                    findflag = True
                    index += 3  # go forward 3 lines

                    line = lines[index]
                    # print  line
                    words = line.split()
                    if words[0] == "sersic":

                        xpos = words[3]
                        ypos = words[4]
                        magser = words[5]
                        reser = words[6]
                        nser = words[7]
                        axratser = words[8]
                        angleser = words[9]

            # removing comma  & brackets

                        xpos = xpos.strip(",")  # equivalent in python
                        ypos = ypos.strip(")")

            # errors

                        index += 1

                        line = lines[index]

                        words = line.split()

                        erxpos = words[1]
                        erypos = words[2]
                        ermagser = words[3]
                        erreser = words[4]
                        ernser = words[5]
                        eraxiser = words[6]
                        erangleser = words[7]

                        erxpos = erxpos.strip(",")  # equivalent in python
                        erypos = erypos.strip(")")

                            # cleaning variables from [] {}

                        xpos = xpos.strip("[")
                        xpos = xpos.strip("]")
                        xpos = xpos.strip("{")
                        xpos = xpos.strip("}")
                        ypos = ypos.strip("[")
                        ypos = ypos.strip("]")
                        ypos = ypos.strip("{")
                        ypos = ypos.strip("}")

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

                        axratser = axratser.strip("[")
                        axratser = axratser.strip("]")
                        axratser = axratser.strip("{")
                        axratser = axratser.strip("}")

                        angleser = angleser.strip("[")
                        angleser = angleser.strip("]")
                        angleser = angleser.strip("{")
                        angleser = angleser.strip("}")

                            # cleaning err variables

                        erxpos = erxpos.strip("[")
                        erxpos = erxpos.strip("]")
                        erxpos = erxpos.strip("{")
                        erxpos = erxpos.strip("}")
                        erypos = erypos.strip("[")
                        erypos = erypos.strip("]")
                        erypos = erypos.strip("{")
                        erypos = erypos.strip("}")

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

                        erangleser = erangleser.strip("[")
                        erangleser = erangleser.strip("]")
                        erangleser = erangleser.strip("{")
                        erangleser = erangleser.strip("}")

                            # Checking if variables have * * values

                        if "*" in magser:
                            magser = magser.strip("*")
                            errflag = True
                            findflag = False

                        if "*" in reser:
                            reser = reser.strip("*")
                            errflag = True
                            findflag = False

                        if "*" in nser:
                            nser = nser.strip("*")
                            errflag = True
                            findflag = False

                        if "*" in axratser:
                            axratser = axratser.strip("*")
                            errflag = True
                            findflag = False

                        if "*" in angleser:
                            angleser = angleser.strip("*")
                            errflag = True
                            findflag = False


                        if ((ermagser == "nan") or (erreser == "nan") or (ernser == "nan") or (erangleser == "nan") or (eraxiser == "nan") or (ermagser == "-nan") or (erreser == "-nan") or (ernser == "-nan") or (erangleser == "-nan") or (eraxiser == "-nan")):
                            errflag = True
                            findflag = False

                        break

            index += 1
#        except:
#            print("Can't open fit.log Sersic \n")
#            fitlogflag = True

#        if errflag == 0 and fitlogflag == 0 and findflag == 1:
    if findflag == True:
        return findflag, float(xpos), float(ypos), float(magser), float(reser), float(nser), float(axratser), float(angleser)
    else:
        return findflag, 0, 0, 99, 0.1, 1, 1, 0


# def FitLogSersic(fitlog, Objsky):
#    "Reads fit.log file and returns Sersic parameters "
    # K
#    pass


def SearchExp(neighid):
    "Search exponential component in fit.log file and returns parameters"
    # Search exponential component
    # in fit.log file
    # and returns parameters

    # outputs = sxpos,sypos,smagser,sreser,snser,saxratser,sangleser

#    flag = False
    errflag = False
#    fitlogflag = False
    findflag = False

# putting some limit to errors

    neighobj = "obj" + "-" + str(neighid)


    # eliminar:

    #consflag = ReadConstraints(constraints,1,consmagmin1,consmagmax1,consremin1,consremax1,consnmin1,consnmax1)
    #consflag = ReadConstraints(constraints,2,consmagmin2,consmagmax2,consremin2,consremax2,consnmin2,consnmax2)

    with open("fit.log") as FITLOG:


            # All lines including the blank ones
        lines = (line.rstrip() for line in FITLOG)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
            # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0
        while index < len(lines):

            line = lines[index]
            (params) = line.split()

            if params[0] == "Init.":     # input image
                if params[4] == neighobj:
                    findflag = True

                    # print line

                    index += 5  # go forward 5 lines

                    line = lines[index]
                        # print  line
                    words = line.split()
                    if words[0] == "expdisk":

                        xpos = words[3]
                        ypos = words[4]
                        magexp = words[5]
                        rsexp = words[6]
                        axratexp = words[7]
                        angleexp = words[8]

            # removing comma  & brackets

                        xpos = xpos.strip(",")  # equivalent in python
                        ypos = ypos.strip(")")

            # errors

                        index += 1

                        line = lines[index]

                        words = line.split()

                        erxposexp = words[1]
                        eryposexp = words[2]
                        ermagexp = words[3]
                        errsexp = words[4]
                        eraxisexp = words[5]
                        erangleexp = words[6]

                        erxposexp = erxposexp.strip(",")  # equivalent in python
                        eryposexp = eryposexp.strip(")")

                            # cleaning variables from [] {}

                        xpos = xpos.strip("[")
                        xpos = xpos.strip("]")
                        xpos = xpos.strip("{")
                        xpos = xpos.strip("}")
                        ypos = ypos.strip("[")
                        ypos = ypos.strip("]")
                        ypos = ypos.strip("{")
                        ypos = ypos.strip("}")

                        magexp = magexp.strip("[")
                        magexp = magexp.strip("]")
                        magexp = magexp.strip("{")
                        magexp = magexp.strip("}")

                        rsexp = rsexp.strip("[")
                        rsexp = rsexp.strip("]")
                        rsexp = rsexp.strip("{")
                        rsexp = rsexp.strip("}")

                        axratexp = axratexp.strip("[")
                        axratexp = axratexp.strip("]")
                        axratexp = axratexp.strip("{")
                        axratexp = axratexp.strip("}")

                        angleexp = angleexp.strip("[")
                        angleexp = angleexp.strip("]")
                        angleexp = angleexp.strip("{")
                        angleexp = angleexp.strip("}")

                            # cleaning err variables

                        erxposexp = erxposexp.strip("[")
                        erxposexp = erxposexp.strip("]")
                        erxposexp = erxposexp.strip("{")
                        erxposexp = erxposexp.strip("}")
                        eryposexp = eryposexp.strip("[")
                        eryposexp = eryposexp.strip("]")
                        eryposexp = eryposexp.strip("{")
                        eryposexp = eryposexp.strip("}")

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

                        erangleexp = erangleexp.strip("[")
                        erangleexp = erangleexp.strip("]")
                        erangleexp = erangleexp.strip("{")
                        erangleexp = erangleexp.strip("}")

                            # Checking if variables have * * values

                        if "*" in magexp:
                            magexp = magexp.strip("*")
                            errflag = True
                            findflag = False


                        if "*" in rsexp:
                            rsexp = rsexp.strip("*")
                            errflag = True
                            findflag = False

                        if "*" in axratexp:
                            axratexp = axratexp.strip("*")
                            errflag = True
                            findflag = False


                        if "*" in angleexp:
                            angleexp = angleexp.strip("*")
                            errflag = True
                            findflag = False


                        if ((ermagexp == "nan") or (errsexp == "nan") or (erangleexp == "nan") or (eraxisexp == "nan") or (ermagexp == "-nan") or (errsexp == "-nan") or (erangleexp == "-nan") or (eraxisexp == "-nan")):
                            errflag = True
                            findflag = False

                        break

            index += 1
    #    except:
    #        print("Can't open fit.log expdisk \n")
    #        fitlogflag = True

#        if errflag == 0 and fitlogflag == 0 and findflag == 1:
    if findflag == True:
        return findflag, float(xpos), float(ypos), float(magexp), float(rsexp), float(axratexp), float(angleexp)
    else:
        return findflag, 0, 0, 99, 0.1, 1, 1, 0


def GALFIT(gfile,out,parvar,objid,n,fout3,fout4):
    "Just run GALFIT and update neighbors"


    sigfile="sigma-"+str(objid) + ".fits"

# sigma file lives at:
    sigfile =  parvar.MaskDir + "/" + sigfile


    if (parvar.Nice == 0):

        runcmd = "galfit -outsig {} > {}".format(gfile, out) #run GALFIT
#        runcmd = "galfit -outsig {}/{} > {}".format(parvar.RunDir,gfile, out) #run GALFIT

        errno = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)
#       CheckError(errcmp)


        runcmd = "mv sigma.fits {}".format(sigfile)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)
#       CheckError(errmv)


        if os.path.isfile("fit.log"):
            findflag= SearchObj(objid)
        else:
            findflag=False


        if (findflag == False):  ## change this number
            line = "{} \n".format(gfile)
            fout3.write(line)
            parvar.Failures+=1

        else:
            line = "{} \n".format(gfile)
            fout4.write(line)
            parvar.Success+=1
##            UpdateNeighbors(objid,n)  # removed
############################################################


    elif(parvar.Nice == 1):

        runcmd = "nice galfit -outsig {} > {}".format(gfile, out) #run GALFIT
#        runcmd = "nice galfit -outsig {}/{} > {}".format(parvar.RunDir,gfile, out) #run GALFIT

        errno = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)
#       CheckError(errcmp)

        runcmd = "mv sigma.fits {}".format(sigfile) #run GALFIT
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)
#       CheckError(errmv)

#        findflag= SearchObj(objid)
        if os.path.isfile("fit.log"):
            findflag= SearchObj(objid)
        else:
            findflag=False



        if (findflag == False):
            line = "{} \n".format(gfile)
            fout3.write(line)
            parvar.Failures+=1
        else:
            line = "{} \n".format(gfile)
            fout4.write(line)
            parvar.Success+=1
#            UpdateNeighbors(objid,n)  # Removed
#######################################

    return True

### Deprecated, not needed
def UpdateNeighbors(index,n):
    "update flags for objects"

#    global Num

    checkflag=False

    flagsat=4
    maxflag=128

    total=len(Num)+1
    overflag=0

    flagser=SearchSersic(index)  # Search component

    (flagexp)=SearchExp(index)

    if (flagser == True):

        for i, item in enumerate(Num):
            upflag=0
            if (i != n):



                ids = Neighbors[i].split("-")  # change for FindNeighbors2

                for j, item2 in enumerate(ids):

                    neighn=ids[j] # neighbors

                    if (neighn == n):
                        upflag=1


                checkflag   = CheckFlag(obj.Flag[i],flagsat,maxflag)

                if ((checkflag == False) and (upflag ==0)):

                    xx=obj.XPos[i] - obj.OFFX[n]  # correct positions for the obj subpanel
                    yy=obj.YPos[i] - obj.OFFY[n]

                    overflag=CheckOverlap(xx,yy,RKron[i],Angle[i],AR[i],xpos,ypos,TimesRe*reser,angleser,axratser)

                    if(flagexp == 0):

                        overflagexp=CheckOverlap(xx,yy,RKron[i],Angle[i],AR[i],xposexp,yposexp,TimesRe*(1.678347)*rsexp,angleexp,axratexp)


                    if (overflag == 1 or overflagexp == 1):

                        Neighbors[i] = Neighbors[i]+"-"+n


    else:
        print ("can't found obj {} for updating Neighbors \n".format(index))



    return True


def SearchObj(idx):
    "Search obj component in fit.log and returns True if it is found"
    # Search sersic component
    # in fit.log file
    # and returns parameters

    # outputs = sxpos,sypos,smagser,sreser,snser,saxratser,sangleser

    flag = False

# putting some limit to errors

    objidx = "obj" + "-" + str(idx)

    objnum=int(idx)

    # eliminar:

    #consflag = ReadConstraints(constraints,1,consmagmin1,consmagmax1,consremin1,consremax1,consnmin1,consnmax1)
    #consflag = ReadConstraints(constraints,2,consmagmin2,consmagmax2,consremin2,consremax2,consnmin2,consnmax2)

    FITLOG = open("fit.log")

    # All lines including the blank ones
    lines = (line.rstrip() for line in FITLOG)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0
    while index < len(lines):

        line = lines[index]
        (params) = line.split()

        if params[0] == "Init.":     # input image
            (name,num)=params[4].split("-")
            num= int(num)
            if num == objnum:
                flag=True
                break

        index += 1

#        except:
#            print("Can't open fit.log \n")

    FITLOG.close()
    return flag


#############################################################################
######################### End of program  ###################################
#############################################################################

#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/


#############################################################################
#############################################################################
#############################################################################
