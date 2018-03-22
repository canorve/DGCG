#!/usr/bin/python3

import numpy as np
import sys
import os
import stat
import subprocess as sp
import os.path
from astropy.io import fits
import scipy
import scipy.special
#from timeit import default_timer as timer



def main():

#    if len(sys.argv[1:]) < 3 or len(sys.argv[1:]) > 4 :
#        print ('Missing arguments')
#        print ("Usage:\n %s [InFile] [ConfigFile] [OutFile] [Header]" % (sys.argv[0]))

#        sys.exit()

#    else:
#        print ("DGCG Version: {}           \n".format(Version))


#    InFile= sys.argv[1]
#    ParamFile= sys.argv[2]
#    OutFile= sys.argv[3]

#    if len(sys.argv[1:]) == 4:
#        Flag=int(sys.argv[4])
#    else:
#        Flag=0

#    if Flag == 1:
#        HeadFlag = True
#    else:
#        HeadFlag = False


#    posxser = posxser  +  obj.OFFX[indx]
#    posyser = posyser  +  obj.OFFY[indx]

#    xfit = obj.XPos[nobj] - obj.OFFX[nobj]
#    yfit = obj.YPos[nobj] - obj.OFFY[nobj]

    pppnum=222

    posxser = 1053.950
    posyser = 1137.98

    posx=1053.576
    posy=1138.105

    xmin = 1028
    xmax = 1079

    ymin = 1111
    ymax = 1165

    offx = 618
    offy = 618

#    xser  = posxser - xmin
#    yser  = posyser - ymin

    xser  = posx
    yser  = posy


    xser  = xser - offx
    yser  = yser - offy


    xmin = xmin - offx
    xmax = xmax - offx


    ymin = ymin - offy
    ymax = ymax - offy


    ReSer =  1.986

    pixreser=ReSer*0.68

    sky = 1065.840

    imsig = "sigma-{}".format(pppnum)
    immask = "seg-3-3.fits"

    pixreexp=0
    bulgetotal =1
    disktotal = 1 - bulgetotal
    bdradius = pixreser * bulgetotal + pixreexp * disktotal
    btflag=False

    imres="A85-004148.23-m0917003.50-222-out.fits"

    rmin = 2
    rmax = 2*bdradius

    FitBox=6



    (tidal,objchinu,bump,snr,ndof) = Tidal(imres, imsig, immask, pppnum,xser, yser, xmin, xmax, ymin, ymax, rmin, rmax,sky,btflag,FitBox)


    outs="Tidal = {}, objchinu = {}, bump = {}, snr = {}, ndof = {}".format(tidal,objchinu,bump,snr,ndof)
    print(outs)




def Tidal(imgout,imsig,immask,num,xser,yser,xlo,xhi,ylo,yhi,rmin,rmax,sky,btflag,fitbox):
    "compute Tidal, local Chinu and Bumpiness values as defined"
    "in Tal et al. 2009 AJ and Blakeslee 2006 ApJ value are defined between"
    "rmin and rmax defined by the user"

# imgout = galfit output image
# imsig = sigma image (The one created by GALFIT)
# immask = Mask image: The original
# xlo, xhi, ylo, yhi = coordinate of object

#    pixflag=True
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

    tidal=0
    objchinu=0
    bump=0
    snr=0

# ATENTION NEED SOME COORDINATE CORRECTION FOR XLO,YLO,XHI,YHI

    xsize = xhi - xlo
    ysize = yhi - ylo

#   enlarge fit area

    xsize = fitbox * xsize
    ysize = fitbox * ysize


    XLo = xser - xsize / 2
    XLo = int(XLo)

    XHi = xser + xsize / 2
    XHi = int(XHi)

    YLo = yser - ysize / 2
    YLo = int(YLo)

    YHi = yser + ysize / 2
    YHi = int(YHi)



## Check eliminate file argument from Tidal

    hdug = fits.open(imgout) # galaxy
    datg = hdug[1].data

#    hdum = fits.open(model) #model
    datm = hdug[2].data

    hdus = fits.open(imsig) #sigma
    dats = hdus[0].data


#### correct here for mask image

    hdu = fits.open(immask) #mask
    dat = hdu[0].data

    maskm = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # tile coordinate

# mask for tidal
    ypos, xpos = np.mgrid[ylo - 1:yhi, xlo - 1:xhi]

    dx = xpos - xser
    dy = ypos - yser

    dist = np.sqrt(dx**2 + dy**2)


# mask including rmin

# tile coordinate:
    masktid = dat[ylo - 1:yhi, xlo - 1:xhi] == num #and ( dist > rmin )  #--ATENTION center do not work

###############

    xlo = xlo - XLo  # transforming coordinates to small output image
    ylo = ylo - YLo

    xhi = xhi - XLo
    yhi = yhi - YLo

    xser = xser - XLo
    yser = yser - YLo

    if masktid.any():

        sumflux  = np.sum(datg[ylo - 1:yhi, xlo - 1:xhi][masktid] - sky) #small coordinate
        sumsig  = np.sum(dats[ylo - 1:yhi, xlo - 1:xhi][masktid]) #small coordinate


        sigflux  = dats[ylo - 1:yhi, xlo - 1:xhi][masktid] # small coordinate
        galflux  = datg[ylo - 1:yhi, xlo - 1:xhi][masktid] # small coordinate
        modflux  = datm[ylo - 1:yhi, xlo - 1:xhi][masktid] # small coordinate

# for tidal, correction for rmin
        galfluxtid  = datg[ylo - 1:yhi, xlo - 1:xhi][masktid] # small coordinate
        modfluxtid  = datm[ylo - 1:yhi, xlo - 1:xhi][masktid] # small coordinate


        resflux=(galflux - modflux)**2

#        res= (galflux - modflux)**2 #DANGER
        varchi=sigflux**2
        chinu= np.sum(resflux/varchi)
        #$sumchinu = $sumchinu + $chinu;

        pixcountchi=np.size(datm[ylo - 1:yhi, xlo - 1:xhi][masktid]) # small coordinate

        if(pixcountchi > 11):
#  local chinu
#  free parameters: sersic = 7  Sersic + disk = 11

            if (btflag==False):
                ndof=pixcountchi-7

            elif(btflag==True):
                ndof=pixcountchi-11


            objchinu= chinu / ndof
        else:
            objchinu=-1
            ndof=-1


        if(np.size(dats[ylo - 1:yhi, xlo - 1:xhi][masktid]) > 0): # small coordinate
# snr
            meanflux = sumflux  / np.size(datg[ylo - 1:yhi, xlo - 1:xhi][masktid]) # small coordinate
            sigma = sumsig / np.size(dats[ylo - 1:yhi, xlo - 1:xhi][masktid]) # small coordinate
            snr = meanflux / sigma

        else:
            snr=-1

# Tidal parameter

        tgal = np.abs((galfluxtid)/(modfluxtid) - 1)
        sumtidal=np.sum(tgal)
        pixcount=np.size(datm[ylo - 1:yhi, xlo - 1:xhi][masktid])  #small coordinate

        if pixcount > 0:
#            tidal = 10 * (tgal / pixcount)  #check the factor of 10
            tidal = (sumtidal / pixcount)
        else:
            tidal=-1

#  Tidal: Check the new formule for Tidal  in Tal. 2009 for correction for large radius


# bumpiness
        resbump  = (galflux - modflux)**2
        sflux    = np.sum(np.abs(modflux - sky))
        #sumsflux = sumsflux + sflux
        varres   = sigflux * sigflux
        numbump  = np.sum(resbump - varres)

#        sumres   = sumres + numbump
#		    $sumvarres= $sumvarres + $varres;

#        $pixcount++;

# Bumpiness
        if pixcount > 0:
            meansflux = sflux/pixcount
            meanres   = numbump/pixcount
#	    $varb      = $sumvarres/$pixcount;
            if (meanres < 0):
                meanres=0

            bump      = (np.sqrt(meanres))/meansflux
        #    bump      = 10*bump  # check this factor of 10

        else:
            bump=-1



###    return True

#    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

#    dx = xpos - x
#    dy = ypos - y

#####

#    hdum.close()
    hdug.close()
    hdus.close()
    hdu.close()


    return (tidal,objchinu,bump,snr,ndof)






#end of program
if __name__ == '__main__':
    main()
