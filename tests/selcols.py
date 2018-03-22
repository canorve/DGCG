#!/usr/bin/python

import numpy as np
import sys
import os
import stat
import subprocess as sp
import os.path
#from astropy.io import fits
import scipy
import scipy.special
from timeit import default_timer as timer



def main():

    if len(sys.argv[1:]) < 3 or len(sys.argv[1:]) > 4 :
        print ('Missing arguments')
        print ("Usage:\n %s [InFile] [ConfigFile] [OutFile] [Header]" % (sys.argv[0]))

        sys.exit()

#    else:
#        print ("DGCG Version: {}           \n".format(Version))


    InFile= sys.argv[1]
    ParamFile= sys.argv[2]
    OutFile= sys.argv[3]

    if len(sys.argv[1:]) == 4:
        Flag=int(sys.argv[4])
    else:
        Flag=0

    if Flag == 1:
        HeadFlag = True
    else:
        HeadFlag = False



    SelectColumns(InFile, ParamFile, OutFile, HeadFlag)





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






#end of program
if __name__ == '__main__':
    main()
