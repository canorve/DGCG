#! /usr/bin/env python3


import numpy as np
import sys
import os
import stat
import subprocess as sp
import os.path
from astropy.io import fits
import scipy
import scipy.special
from timeit import default_timer as timer
import argparse 

from dgcg.lib import core 
from dgcg.lib import check
from dgcg.lib import image
from dgcg.lib import output
from dgcg.lib import catfil
from dgcg.lib import config 


############################################################################################
#                                      DGCG                                                #
#                                                                                          #
#                      Driver for GALFIT on Cluster Galaxies                               #
#                                                                                          #
#                       DGCG is a wrapper script for GALFIT                                #
#                                                                                          #
#                                                                                          #
#                          written by C. Añorve                                            #
#                                                                                          #
############################################################################################

# Last Version: 5/Jul/2011
# MaxFit and MagDiff modified

# Last Version: 10/Mar/2011
# number of degrees of freedom is added to output

# Last Version: 1/Mar/2011
# flag = 4 is activated when a nan error is occurred

# Last Version: 24/Feb/2011
# a new pixel directory was added. Now the final parameters are computed within 1 kron radius

# Last Version: 19/Feb/2011
# JoinSexOut was modified. Now provide 99 mag to galaxies which were not fitted by GALFIT

# Last Version: 28/Jan/2011
# MakeOutput and Tidal functions were corrected. principally the computation of the SNR

# Last Version: 9/Dic/2010
# CheckSatRegion function was corrected. Also run time is print before to create final files

# Last Version: 4/Nov/2010
# Tidal, bumpiness, local chinu, snr are computed inside the Scale * Kron radius

# Last Version: 17/Oct/2016
# Program was translated to Python
# Now it is called pyDGCG

# Last Version: April/2017
# files and folders have a tree structure
# variables are stored in a Class




################
#  main code:  #
################

def mainDGCG():

    global StartRun,EndRun

    # starting to count time
    StartRun = timer()


    parser = argparse.ArgumentParser(description="DGCG: Driver for GALFIT on Cluster Galaxies")

    # required arguments
    parser.add_argument("ConfigFile", help="DGCG configuration file ")


    args = parser.parse_args()

    InFile = args.ConfigFile 

    ParVar = config.read_config(InFile)

    #   initialize default variables
    #################################################
    #   creating a Object for input file parameters
    #ParVar = core.ParamFile()
    #################################################

    # read parameter file
    #catfil.ReadFile(ParVar,InFile)

    # verify parameters have sane values
    # temporalmente desabilitado
    #check.CheckSaneValues(ParVar)


    #########################################
    ### deleting any previous files run by DGCG
    #mover a una funcion

    runcmd = "rm {}".format(ParVar.Crashes)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm {}".format(ParVar.Fitted)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm inforegion"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm fits.log"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm fitlog.dgcg"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm objflags"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm psf.temp"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm {}".format(ParVar.OffsetPos)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm {}".format(ParVar.SexSort)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm {}".format(ParVar.SexArSort)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm {}".format(ParVar.SkyCrashes)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)

    #    runcmd = "rm {}".format(ParVar.SkyFitted)
    #    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
    #                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm {}".format(ParVar.ListObjs)
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm {}*".format(ParVar.FileOut)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm galfit.*"
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm {}*".format(ParVar.Ds9OutName)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)


    #    runcmd = "rm -r {}".format(ParVar.InputDir)
    #    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
    #                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm -r {}".format(ParVar.SkyDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm -r {}".format(ParVar.OutputDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm -r {}".format(ParVar.TempDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm -r {}".format(ParVar.RunDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)

    runcmd = "rm -r {}".format(ParVar.ReRunDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm -r {}".format(ParVar.MaskDir)
    errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)


    runcmd = "rm fit.log"
    errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                   stderr=sp.PIPE, universal_newlines=True)



#######################################

    import pdb;pdb.set_trace()

    (ParVar.NCol, ParVar.NRow) = image.GetAxis(ParVar.Img)
    (ParVar.ExpTime)           = image.GetExpTime(ParVar.Img)
    (ParVar.Gain)              = image.GetGain(ParVar.Img)
#    (ParVar.Rdnoise)           = image.GetRdnoise(ParVar.Img)


    ParVar.Total = catfil.CatArSort(ParVar)

    ParVar.Total = catfil.CatSort(ParVar)


##### segmentation mask

    if ((not(os.path.isfile(ParVar.SegFile))) or (ParVar.Overwrite == 1)):

        image.MakeImage(ParVar.SegFile, ParVar.NCol, ParVar.NRow)

#        if (ParVar.AutoSatRegion == 1):
#            catfil.MakeSatDs9(ParVar)

#       OldMakeMask(SegFile, SexArSort, KronScale,Ds9SatReg)
        image.MakeMask(ParVar.SegFile, ParVar.SexArSort, ParVar.KronScale,0,ParVar.Ds9SatReg)  # offset set to 0
        image.MakeSatBox(ParVar.SegFile, ParVar.Ds9SatReg, ParVar.Total + 1, ParVar.NCol, ParVar.NRow)
    else:
        print("Using old mask image {} \n".format(ParVar.SegFile))


########### Sky segmentation annuli Mask  ##########

    if ((not(os.path.isfile(ParVar.SkyFile))) or (ParVar.Overwrite == 1)):

        image.MakeImage(ParVar.SkyFile, ParVar.NCol, ParVar.NRow)

#        if (ParVar.AutoSatRegion == 1):
#            catfil.MakeSatDs9(ParVar)

#        MakeSkyMask(SkyFile, SexArSort, SkyScale, Offset, Ds9SatReg) # Old

        image.MakeMask(ParVar.SkyFile, ParVar.SexArSort, ParVar.SkyScale, ParVar.Offset, ParVar.Ds9SatReg)
        image.MakeSatBox(ParVar.SkyFile, ParVar.Ds9SatReg, ParVar.Total + 1, ParVar.NCol, ParVar.NRow)


    else:
        print("Using old mask image {} \n".format(ParVar.SkyFile))


    catfil.UpdateSatFlags(ParVar)

################ THIS SECTION WAS REMOVED BECAUSE PixFile IS NOT NEEDED ANYMORE ##################
# This section is needed to compute final parameters
# check how this will be affected at the end

#    if ((not(os.path.isfile(PixFile))) or (Overwrite == 1)):

#        MakeImage(PixFile, NCol, NRow)

#        if (AutoSatRegion == 1):
#            MakeSatDs9(SexCat, SatRegionScale, NCol, NRow, Ds9SatReg)

#       OldMakeMask(PixFile, SexArSort, 1,Ds9SatReg)  # 1 Kron radius
#        MakeMask(PixFile, SexArSort, 1,0,Ds9SatReg)  # 1 Kron radius
#        MakeSatBox(PixFile, Ds9SatReg, Total + 1, NCol, NRow)
#    else:
#        print("Using old mask image {} \n".format(PixFile))

###################


### creating directories...

#    if not os.path.exists(ParVar.InputDir):
#        os.makedirs(ParVar.InputDir)

    if not os.path.exists(ParVar.MaskDir):
        os.makedirs(ParVar.MaskDir)

    if not os.path.exists(ParVar.OutputDir):
        os.makedirs(ParVar.OutputDir)

    if not os.path.exists(ParVar.SkyDir):
        os.makedirs(ParVar.SkyDir)

    if not os.path.exists(ParVar.RunDir):
        os.makedirs(ParVar.RunDir)


#    if not os.path.exists(ParVar.PixDir):
#        os.makedirs(ParVar.PixDir)

#    if not os.path.exists(ParVar.MaskPixDir):
#        os.makedirs(ParVar.MaskPixDir)


    XChunk = int(ParVar.NCol / ParVar.Split)
    YChunk = int(ParVar.NRow / ParVar.Split)


# REMOVED NOT NEEDED ANYMORE
# mask pixels

#    print ("Getting pixels from every object to remove masks \n")
#    GetPixels(SegFile,SexSort,KronScale,PixPrefix,XChunk,YChunk,Buffer)

#    errpix = system("mv PixPrefix* MaskPixDir/\. ")
#    CheckError(errpix)

# object pixels

#    print ("Getting pixels from every object to compute final parameters \n")
#    GetPixels(PixFile,SexSort,1,PixPrefix,XChunk,YChunk,Buffer)

#    errpix = system("mv PixPrefix* PixDir/\. ")
#    CheckError(errpix)
####################################################

    print("Creating Ds9 Box region of all objects \n")
    catfil.BoxDs9(ParVar)


# splitting image files
    print("Splitting images \n")
    image.SplitImage(ParVar)

#    f = open("psf.temp", "w")
#    optarg = "{}/PSF*.fits".format(PsfDir)
#    errpsf = sp.call(["ls", optarg], stdout=f)

    runcmd = "ls {}/PSF*.fits > psf.temp".format(ParVar.PsfDir)
    errpsf = sp.run([runcmd], shell=True, stdout=sp.PIPE,
               stderr=sp.PIPE, universal_newlines=True)


#    f.close()
# CheckError(errpsf)



    print("====================== DGCG In a NutShell =========================== \n")
    print(" This script takes the pySEx output, and formats it for GALFIT.  \n")
    print(" DGCG: Driver for GALFIT on Cluster Galaxies                          \n")
    print(" Created by Christopher Añorve et al.                         	             \n")
    print("===================================================================== \n")


    flog = open(ParVar.LogFile, "w")

# print DGCG options in Log file

    output.PrintVar(ParVar,flog)

#    OffsetFile = "OffsetPos"

    fobjs = open(ParVar.ListObjs, "w")

    fout2 = open(ParVar.OffsetPos, "w")

    fout3 = open(ParVar.Crashes, "w")
    fout4 = open(ParVar.Fitted, "w")


#    fskycrash = open(dgcg.SkyCrashes, "w")
#    fskyfit = open(ParVar.SkyFitted, "w")


############################
# defining a object Class to store all the variables
    Obj = core.Object()
############################

###############  Read in sextractor sorted data ################

    Obj.Num, Obj.RA, Obj.Dec, Obj.XPos, Obj.YPos, Obj.Mag, Obj.Kron, Obj.FluxRad, Obj.IsoArea, Obj.AIm, Obj.E, Obj.Theta, Obj.Background, Obj.Class, Obj.Flag, Obj.XMin, Obj.XMax, Obj.YMin, Obj.YMax, Obj.XSMin, Obj.XSMax, Obj.YSMin, Obj.YSMax = np.genfromtxt(
	    ParVar.SexSort, delimiter="", unpack=True)  # sorted

##
    Obj.Angle = Obj.Theta - 90
    Obj.AR = 1 - Obj.E
    Obj.RKron = ParVar.KronScale * Obj.AIm * Obj.Kron
    Obj.Sky = Obj.Background

    Obj.Num = Obj.Num.astype(int)
    Obj.Flag = Obj.Flag.astype(int)

# other stuff:

#    Tot = len(Obj.Num)
    Tot = ParVar.Total

    Obj.Sersic = [ParVar.NSer] * Tot
    Obj.Sersic = np.array(Obj.Sersic)

    Obj.RSky = ParVar.SkyScale * Obj.AIm * Obj.Kron + ParVar.Offset + ParVar.SkyWidth
    Obj.RKron = ParVar.KronScale * Obj.AIm * Obj.Kron


    masky = Obj.RSky <= 0
    if masky.any():
        Obj.RSky[masky] = 1

    maskron = Obj.RKron <= 0
    if maskron.any():
        Obj.RKron[maskron] = 1

    Obj.SkyFlag = [True] * Tot
    Obj.SkyFlag = np.array(Obj.SkyFlag)

    Obj.Neighbors = Obj.Num


# subpanel stuff:

    Obj.IX = (Obj.XPos / XChunk) + 1
    Obj.IY = (Obj.YPos / YChunk) + 1

    Obj.IX = Obj.IX.astype(int)
    Obj.IY = Obj.IY.astype(int)


    Obj.XMin = Obj.XMin.astype(int)
    Obj.XMax = Obj.XMax.astype(int)
    Obj.YMin = Obj.YMin.astype(int)
    Obj.YMax = Obj.YMax.astype(int)

    Obj.XSMin = Obj.XSMin.astype(int)
    Obj.XSMax = Obj.XSMax.astype(int)
    Obj.YSMin = Obj.YSMin.astype(int)
    Obj.YSMax = Obj.YSMax.astype(int)



    maskblkx = Obj.IX > ParVar.Split
    if maskblkx.any():
        Obj.IX[maskblkx]= Obj.IX[maskblkx] -1


    maskblky = Obj.IY > ParVar.Split
    if maskblky.any():
        Obj.IY[maskblky]= Obj.IY[maskblky] -1


#   Make sure the object coordinate in the subpanel is correct

# create arrays

    Obj.XBuffer=np.array([0]*Tot)
    Obj.YBuffer=np.array([0]*Tot)


    maskix = Obj.IX == 1
    if maskix.any():
        Obj.XBuffer[maskix] = 0

    maskix = Obj.IX != 1
    if maskix.any():
        Obj.XBuffer[maskix] = ParVar.Buffer

    maskiy = Obj.IY == 1
    if maskiy.any():
        Obj.YBuffer[maskiy] = 0

    maskiy = Obj.IY != 1
    if maskiy.any():
        Obj.YBuffer[maskiy] = ParVar.Buffer

# Obj.OFFX and Obj.OFFY transform coordinates
#  from big image to tile image --> IM-X-Y.fits
    Obj.OFFX = (Obj.IX - 1) * XChunk - Obj.XBuffer
    Obj.OFFY = (Obj.IY - 1) * YChunk - Obj.YBuffer


##############################################################
##############################################################
# creating empty arrays

#    Obj.gXMIN = np.array([0]*Tot)
#    Obj.gXMAX = np.array([0]*Tot)
#    Obj.gYMIN = np.array([0]*Tot)
#    Obj.gYMAX = np.array([0]*Tot)


    XSize = Obj.XMax - Obj.XMin
    YSize = Obj.YMax - Obj.YMin

#   enlarge fit area

    XSize = ParVar.FitBox * XSize
    YSize = ParVar.FitBox * YSize

##  30 pixels is the minimum area to fit (arbitrary number):

    masksize = XSize < 30

    if masksize.any():
        XSize[masksize] = 30

    masksize = YSize < 30

    if masksize.any():
        YSize[masksize] = 30

    #  Calculate the (x,y) position of the current object relative to
    #  the tile in which it lives.

    XFit = Obj.XPos - Obj.OFFX
    YFit = Obj.YPos - Obj.OFFY

    #  Calculate fitting box needed to plug into galfit header:

    XLo = XFit - XSize / 2
    XLo = XLo.astype(int)

    maskxy = XLo <= 0
    if maskxy.any():
        XLo[maskxy] = 1


    XHi = XFit + XSize / 2
    XHi = XHi.astype(int)

    maskxy = XHi > ParVar.NCol  #This does not affect the code at all
    if maskxy.any():
        XHi[maskxy] = ParVar.NCol


    YLo = YFit - YSize / 2
    YLo = YLo.astype(int)

    maskxy = YLo <= 0
    if maskxy.any():
        YLo[maskxy] = 1


    YHi = YFit + YSize / 2
    YHi = YHi.astype(int)

    maskxy = YHi > ParVar.NRow  # same as above but for y axis
    if maskxy.any():
        YHi[maskxy] = ParVar.NRow


    # Calculate information needed to plug back into the image header
    # at the end of the fit.

    Obj.gXMIN = XLo + Obj.OFFX  # The [gXMIN:gXMAX,gYMIN:gYMAX] of the box
    Obj.gXMAX = XHi + Obj.OFFX  # relative to the big image from which
    Obj.gYMIN = YLo + Obj.OFFY  # the current tile was extracted.
    Obj.gYMAX = YHi + Obj.OFFY

    Obj.gOFFX = Obj.gXMIN - 1
    Obj.gOFFY = Obj.gYMIN - 1


############################################################
###########################################################
# XLo, YLo, XHi, YHi correspond to the small image
# gXMIN, gYMIN, gXMAX, gYMAX correspond to the big image (original)
##############################################################
##############################################################


# Creating empty arrays for output variables
    Obj.ra = np.array(["none"]*Tot, dtype=object)
    Obj.dec = np.array(["none"]*Tot, dtype=object)

    Obj.InputImage = np.array(["none"]*Tot,dtype=object)
    Obj.OutIm = np.array(["none"]*Tot,dtype=object)
    Obj.Objx = np.array(["none"]*Tot,dtype=object)
    Obj.PPPnum = np.array([0]*Tot)
    Obj.Restart = np.array(["none"]*Tot,dtype=object)
    Obj.Cluster = np.array(["none"]*Tot,dtype=object)


    Obj.PosXSer = np.array([0.0]*Tot)
    Obj.ErPosXSer = np.array([0.0]*Tot)
    Obj.PosYSer = np.array([0.0]*Tot)
    Obj.ErPosYSer = np.array([0.0]*Tot)
    Obj.MagSer =  np.array([99.9]*Tot)
    Obj.ErMagSer = np.array([0.0]*Tot)
    Obj.ReSer = np.array([0.0]*Tot)
    Obj.ErReSer = np.array([0.0]*Tot)
    Obj.NSer = np.array([0.0]*Tot)
    Obj.ErNSer = np.array([0.0]*Tot)
    Obj.AxisSer = np.array([0.0]*Tot)
    Obj.ErAxisSer = np.array([0.0]*Tot)
    Obj.PaSer = np.array([0.0]*Tot)
    Obj.ErPaSer = np.array([0.0]*Tot)
    Obj.KSer = np.array([0.0]*Tot)

    Obj.MeanMeSer = np.array([99.9]*Tot)
    Obj.ErMeanMeSer = np.array([99.9]*Tot)
    Obj.MeSer = np.array([99.9]*Tot)
    Obj.ErMeSer = np.array([99.9]*Tot)


    Obj.PosXExp = np.array([0.0]*Tot)
    Obj.ErPosXExp = np.array([0.0]*Tot)
    Obj.PosYExp = np.array([0.0]*Tot)
    Obj.ErPosYExp = np.array([0.0]*Tot)

    Obj.MagExp = np.array([99.9]*Tot)
    Obj.ErMagExp = np.array([99.9]*Tot)

    Obj.RsExp = np.array([0.0]*Tot)
    Obj.ErRsExp = np.array([0.0]*Tot)
    Obj.AxisExp = np.array([0.0]*Tot)
    Obj.ErAxisExp = np.array([0.0]*Tot)
    Obj.PaExp = np.array([0.0]*Tot)
    Obj.ErPaExp = np.array([0.0]*Tot)
    Obj.MeanMeExp =np.array([99.9]*Tot)
    Obj.MeExp = np.array([99.9]*Tot)
    Obj.ErMeanMeExp =np.array([99.9]*Tot)
    Obj.ErMeExp = np.array([99.9]*Tot)
    Obj.MsExp = np.array([99.9]*Tot)
    Obj.ErMsExp = np.array([99.9]*Tot)

    Obj.Sky = np.array([0.0]*Tot)
    Obj.ErSky = np.array([0.0]*Tot)

    Obj.MagTotal = np.array([99.9]*Tot)
    Obj.ErMagTotal = np.array([99.9]*Tot)

    Obj.BulgeTotal = np.array([0.0]*Tot)
    Obj.ErBt = np.array([0.0]*Tot)
    Obj.ChiNu = np.array([0.0]*Tot)
    Obj.Tidal = np.array([0.0]*Tot)
    Obj.ObjChiNu = np.array([0.0]*Tot)
    Obj.Bump = np.array([0.0]*Tot)
#    Obj.MeanMesky = np.array([99.9]*Tot) # eliminado
    Obj.SNR = np.array([0.0]*Tot)
    Obj.NDof = np.array([0]*Tot)
    Obj.FitFlag = np.array([0]*Tot)
    Obj.ReFit =  np.array([False]*Tot)   #rerun if didn't fit
###########################


#################
# posible to remove
#    print("Finding Neighbors for every object \n")

#    FindNeighbors()  # maybe remove this function
###########

# 0: Sort catalog and create segmentation images
# 1: Execute everything compute sky, create input files, run galfit and create output file;
# 2: only computes sky
# 3: Execute everything except creation of output file


    if ParVar.Execute != 0:

        #   Sky fitting
        #   DGCG computes the sky first and leaves it fixed for galaxy fitting.

        print("DGCG is going to compute sky \n")


# remove fit.log before to compute sky
        runcmd = "rm fit.log"
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)

# remove fit.log before to compute sky
        runcmd = "rm galfit.*"
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)



###############
        core.RunSky(ParVar,Obj)
###############


        runcmd = "mv fit.log {}/fit.log.sky".format(ParVar.SkyDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

        runcmd = "mv Sky-* {}/.".format(ParVar.SkyDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)


        runcmd = "mv galfit.* {}/.".format(ParVar.SkyDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

        print("Done sky fitting \n")


        if (ParVar.Execute != 2):

#########################################################
            # here comes the serious stuff:
            core.DGCG(ParVar,Obj,flog,fobjs,fout2,fout3,fout4)
##########################################################


            print("Done GALFITting :) \n")



#################################################################

# Closing files :

    flog.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fobjs.close()

    #fskycrash.close()
    #fskyfit.close()

    outflag=False
    finalflag=False
    posflag=False  #change to True
    joinflag=False

    if (ParVar.Execute != 3 and ParVar.Execute != 2 and ParVar.Execute != 0):

        print("Creating outputs files \n")


        outflag = output.ReadFitlog2("fitlog.dgcg")

        if (outflag == True):

            if (ParVar.FitFunc == "BD"):
                ParVar.BtFile = ParVar.FileOut + ".bd"

            else:
                ParVar.BtFile = ParVar.FileOut + ".ser"

            finalflag = output.MakeOutput("fitlog.dgcg", ParVar, Obj)

            if (finalflag == True):

#                posflag = output.PosCor(OffsetFile, ParVar.BtFile)  # correct positions ## Deprecated

                posflag=True

                if (posflag == True):

#                    ParVar.SexOut = ParVar.FileOut + ".dgcg"
                    ParVar.SexOut = ParVar.BtFile + ".dgcg"

                    # Join output file with sextractor catalog only for the fitted
                    # objects


                    joinflag = output.JoinSexOut(ParVar, Obj)


                    FilColOut = ParVar.FileOut + ".dgcg"

                    outcolflag=output.SelectColumns(ParVar.SexOut, ParVar.ColPar, FilColOut, ParVar.HeadFlag)


                    if (joinflag == True):

                        print("Ascii output file was succesfully created\n")

                        ParVar.TableFits = ParVar.FileOut + ".dgcg.fits"

                        print("Creating table fits file.. \n")

                        outfitsflag = output.Ascii2Fits(ParVar.SexOut, ParVar.TableFits)

                        print("Table fits file created\n")



                    else:
                        print("Can't create Ascii final output file \n")

                else:
                    print("Can't create dgcg output file \n")

            else:
                print("Can't create bt output file \n")

        else:
            print("Can't create output files \n")


# erasing files

    if (ParVar.Erase == 1):

        print("Erasing unnecesary files \n")

#        runcmd = "rm -r {}".format(ParVar.TempDir)
#        errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                      stderr=sp.PIPE, universal_newlines=True)
        # CheckError(errm)

        os.remove("psf.temp")

        os.remove("inforegion")

        os.remove("ListObjs")

        os.remove("objflags")

        os.remove(ParVar.OffsetPos)

        os.remove(ParVar.SkyCrashes)

        os.remove(ParVar.Crashes)
        os.remove(ParVar.Fitted)




#   $errno = system("rm $InputDir/mask-*");
#   CheckError($errno);

        runcmd = "rm {}/Sky-*".format(ParVar.SkyDir)
        errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)
        # CheckError(errm)

        runcmd = "rm {}/galfit.*".format(ParVar.SkyDir)
        errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)
        # CheckError(errm)


        runcmd = "rm -r {}".format(ParVar.SkyDir)
        errm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True)


## removed after the introduction of RunDir
#            runcmd = "mv obj-* {}/.".format(ParVar.InputDir)
#            errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                           stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)
#####################


        runcmd = "mv sigma-* {}/.".format(ParVar.MaskDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

#            runcmd = "mv out-* {}/.".format(ParVar.OutputDir)
#            errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                           stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

        runcmd = "mv galfit.* {}/.".format(ParVar.OutputDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

#            runcmd = "mv mask-* {}/.".format(ParVar.InputDir)
#            errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                           stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

        runcmd = "mv *-out.fits {}/.".format(ParVar.OutputDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errmv)

        runcmd = "mv *.png {}/.".format(ParVar.OutputDir)
        errmv = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                       stderr=sp.PIPE, universal_newlines=True)

#            mkdir = "{}/{}".format(ParVar.InputDir, ParVar.PsfDir)
#            if not os.path.exists(mkdir):
#                os.makedirs(mkdir)

#            runcmd = "cp {}/* {}/.".format(ParVar.PsfDir, mkdir)
#            errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                           stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errcp)

#            runcmd = "cp {} {}/.".format(ParVar.ConsFile, ParVar.InputDir)
#            errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,
#                           stderr=sp.PIPE, universal_newlines=True)
#        CheckError(errcp)



    if (ParVar.Execute != 0 and ParVar.Execute != 2):

        GalTot = ParVar.Failures + ParVar.Success

        print("DGCG had {} success out of a total of {} \n".format(ParVar.Success, GalTot))



    if not os.path.exists(ParVar.ReRunDir):
        os.makedirs(ParVar.ReRunDir)


    maskrun= Obj.ReFit == True


    ind=np.where(maskrun == True)

    indx=ind[0]

    for idx in enumerate(indx):

        objid = Obj.Num[idx[1]]

        parfile =  "obj" + "-" + str(objid)

        pmsg="copying file {} for refit. ReFit = {}; FitFlag = {} ".format(parfile, Obj.ReFit[idx[1]], Obj.FitFlag[idx[1]])
        print(pmsg)

        dirparfile =  ParVar.RunDir + "/" + parfile

        redirparfile =  ParVar.ReRunDir + "/" + parfile

        runcmd = "cp {} {}".format(dirparfile,redirparfile)
        errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

    print("copying files to refit done..")


################################
################################
################################

# computing running time

    EndRun = timer()

    RunTime = EndRun - StartRun

    RunTimeMin = RunTime/60

    RunTimeHour = RunTimeMin/60

    RunTimeDays = RunTimeHour/24


    if (RunTimeDays >= 1 ):

        print ("Job took {0:.2f} days \n".format(RunTimeDays))


    elif (RunTimeHour >= 1):

        print ("Job took {0:.2f} hours \n".format(RunTimeHour))


    elif (RunTimeMin >= 1):

        print ("Job took {0:.2f} minutes \n".format(RunTimeMin))

    else:
        print ("Job took {0:.2f} seconds \n".format(RunTime))



    print ("Done everything! \n")



#############################################################################
#              End of program  ###################################


#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/


##############################################################################


