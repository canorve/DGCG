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

    fitbox=6

#    xser  = posxser - xmin
#    yser  = posyser - ymin

    xser  = posx
    yser  = posy

#    xser  = xser - offx
#    yser  = yser - offy

#    xmin = xmin - offx
#    xmax = xmax - offx

#    ymin = ymin - offy
#    ymax = ymax - offy

    xsize = xmax - xmin
    ysize = ymax - ymin

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

#   coordinate transformation big image --> galfit image

    goffx = XLo - 1
    goffy = YLo - 1

    ReSer =  1.986

    pixreser=ReSer*0.68

    sky = 1065.840

    imsig = "sigma-{}".format(pppnum)
    immask = "seg-3-3.fits"

    immask = "seg.fits"  # large image

    pixreexp=0
    bulgetotal =1
    disktotal = 1 - bulgetotal
    bdradius = pixreser * bulgetotal + pixreexp * disktotal
    btflag=False

    imres="A85-004148.23-m0917003.50-222-out.fits"

    rmin = 2
    rmax = 2*bdradius


    (tidal,objchinu,bump,snr,ndof) = Tidal(imres, imsig, immask, pppnum, xser, yser, xmin, xmax, ymin, ymax, goffx, goffy, rmin, sky, btflag)

    outs="Tidal = {}, objchinu = {}, bump = {}, snr = {}, ndof = {}".format(tidal,objchinu,bump,snr,ndof)
    print(outs)



def Tidal(imgout,imsig,immask,num,xser,yser,xlo,xhi,ylo,yhi,goffx,goffy,rmin,sky,btflag):
    "compute Tidal, local Chinu and Bumpiness values as defined"
    "in Tal et al. 2009 AJ and Blakeslee 2006 ApJ value are defined between"
    "rmin and rmax defined by the user"

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

#        print(maskbum)
#        print(xglo,xghi,yglo,yghi)

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




#end of program
if __name__ == '__main__':
    main()
