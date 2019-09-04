
import os.path
import os
from astropy.io import fits
import numpy as np
import sys
import scipy
import scipy.special
import subprocess as sp

from lib import check



def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
#    ncol = hdu[0].header["NAXIS1"]
#    nrow = hdu[0].header["NAXIS2"]
    ncol = hdu[0].header.get("NAXIS1",2000) # return 2000 if not found
    nrow = hdu[0].header.get("NAXIS2",2000) # return 2000 if not found

    hdu.close()
    return ncol, nrow


def GetGain(Image):
    # k Check
    "Get gain from the image"

    hdu = fits.open(Image)
#    gain = hdu[0].header["GAIN"]
    gain = hdu[0].header.get("GAIN",7) # return 7 if not found

    hdu.close()
    return gain


def GetRdnoise(Image):
    # k Check
    "Get Rdnoise from the image"

    hdu = fits.open(Image)
#    rdnoise = hdu[0].header["RDNOISE"]
    rdnoise = hdu[0].header.get("RDNOISE",5.2) # return 5.2 if not found

    hdu.close()
    return rdnoise


def GetExpTime(Image):
    # k Check
    "Get exposition time from the image"

    hdu = fits.open(Image)
#    exptime = hdu[0].header["EXPTIME"]
    exptime = hdu[0].header.get("EXPTIME",1) # return 1 if not found

    hdu.close()
    return exptime


############################################# splitimage.pl ##############


def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"
# k Check


    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.open(Image)
    dat = hdu[0].data[ylo - 1:yhi, xlo - 1:xhi]
    hdu[0].data = dat
    hdu.writeto(Imageout, overwrite=True)
    hdu.close()



def SplitImage(parvar):
    "Cut image in NSPLIT pieces "

# Removing ".fits" part

    if parvar.Img.find(".") != -1:
        (TIMG, trash) = parvar.Img.split(".")
    else:
        TIMG=parvar.Img

    if parvar.SegFile.find(".") != -1:
        (TMASK, trash) = parvar.SegFile.split(".")
    else:
        TMASK=parvar.SegFile

    if parvar.SkyFile.find(".") != -1:
        (TSKY, trash) = parvar.SkyFile.split(".")
    else:
        TSKY=parvar.SkyFile

    if parvar.SigImg.find(".") != -1:
        (TSIG, trash) = parvar.SigImg.split(".")
    else:
        TSIG=parvar.SigImg

    if not os.path.isdir(parvar.TempDir):
        os.mkdir(parvar.TempDir)

    print("\n")
    print("Putting all sub-panel images into the directory: {}\n".format(parvar.TempDir))
    print("\n")


################  Split image into nsplit panels along a side  ################

    nsubx = int(parvar.NCol / parvar.Split)
    nsuby = int(parvar.NRow / parvar.Split)

# Get rid of the last / in the directory name:

    for j in range(1, parvar.Split + 1):

        for i in range(1, parvar.Split + 1):

            if (nsubx * (i - 1) - parvar.Buffer) <= 1:
                xlo = 1
            else:
                xlo = (nsubx * (i - 1) - parvar.Buffer + 1)

            if (nsuby * (j - 1) - parvar.Buffer) <= 1:
                ylo = 1
            else:
                ylo = (nsuby * (j - 1) - parvar.Buffer + 1)
            if (nsubx * i + parvar.Buffer >= parvar.NCol):
                xhi = parvar.NCol
            else:
                xhi = (nsubx * i + parvar.Buffer)
            if (nsuby * j + parvar.Buffer >= parvar.NRow):
                yhi = parvar.NRow
            else:
                yhi = (nsuby * j + parvar.Buffer)

# make temporarily fits files.

            imgout = "{}/{}-{}-{}.fits".format(parvar.TempDir, TIMG, i, j)
            maskout = "{}/{}-{}-{}.fits".format(parvar.TempDir, TMASK, i, j)
            skyout = "{}/{}-{}-{}.fits".format(parvar.TempDir, TSKY, i, j)
            sigout = "{}/{}-{}-{}.fits".format(parvar.TempDir, TSIG, i, j)


#            imgout = "{}-{}-{}.fits".format(TIMG, i, j)
#            maskout = "{}-{}-{}.fits".format(TMASK, i, j)
#            skyout = "{}-{}-{}.fits".format(TSKY, i, j)
#            sigout = "{}-{}-{}.fits".format(TSIG, i, j)



            if os.path.isfile(parvar.Img):
                if not os.path.isfile(imgout):
                    GetFits(parvar.Img, imgout, xlo, xhi, ylo, yhi)
                    # falta CheckError

            if os.path.isfile(parvar.SegFile):
                if not os.path.isfile(maskout):
                    GetFits(parvar.SegFile, maskout, xlo, xhi, ylo, yhi)
                    # falta CheckError

            if os.path.isfile(parvar.SkyFile):
                if not os.path.isfile(skyout):
                    GetFits(parvar.SkyFile, skyout, xlo, xhi, ylo, yhi)
                    # falta CheckError

            if os.path.isfile(parvar.SigImg):
                if not os.path.isfile(sigout):
                    GetFits(parvar.SigImg, sigout, xlo, xhi, ylo, yhi)
                    # falta CheckError

    return True



def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"
# k Check

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True

######################
## Deprecated function
######################
def OldMakeMask(maskimage, catfile, scale, regfile):
    "Create a mask image using ellipses for every Object of catfile. Deprecated"
# k Check  Deprecated
# scale=Kronscale

    checkflag = False
    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Mask for neighboring contamination  \n")

    Rkron = scale * ai * kr

    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    for idx, val in enumerate(n):

        # check if object doesn't has saturated regions
        checkflag = check.CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by user
        # in ds9
        regflag = check.CheckSatReg(xx[idx], yy[idx], regfile, Rkron[
                              idx], theta[idx], e[idx])

        if (checkflag == False) and (regflag == False):

            print("Creating ellipse mask for obj {}  \n".format(n[idx]))

            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                           idx], sxmin[idx], sxmax[idx], symin[idx], symax[idx])

        elif(checkflag == True) or (regflag == True):

            print(
                "Skipping object {}, one or more pixels are saturated \n".format(n[idx]))

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeMask(maskimage, catfile, scale, offset, regfile):
    "Create a mask image using ellipses for every Object of catfile. Now includes offset"
# k Check

    checkflag = False
    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Masks for neighbour contamination \n")

    Rkron = scale * ai * kr + offset
    skron = ai * kr #+ offset  # small elipse to check saturation



    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1


    mask = skron < 1
    if mask.any():
        skron[mask] = 1



    hdu = fits.open(maskimage)
    img = hdu[0].data

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = check.CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
#        regflag = check.CheckSatReg(xx[idx], yy[idx], regfile, Rkron[
#            idx], theta[idx], e[idx])  # old with Rkron radius

        regflag = check.CheckSatReg(xx[idx], yy[idx], regfile, skron[
            idx], theta[idx], e[idx])


        if (checkflag == False) and (regflag == False):

            print ("Creating ellipse annuli for object {}  \n".format(n[idx]))
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):

            print ("Skipping object {}, one or more pixels are saturated \n".format(n[idx]))




    if os.path.isfile(maskimage):
        print("{} deleted; a new one is created \n".format(maskimage))
        runcmd = "rm {}".format(maskimage)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

# Check

    xmin = np.int(xmin)
    xmax = np.int(xmax)
    ymin = np.int(ymin)
    ymax = np.int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist < dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


def MakeSatBox(maskimage, region, val, ncol, nrow):
    "Create a mask for saturated regions"
    "Regions must be in DS9 box regions format"

# k Check

#	fileflag=1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    with open(region) as f_in:

        next(f_in)
        next(f_in)
        next(f_in)

        # All lines including the blank ones
        lines = (line.rstrip() for line in f_in)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            (box, info) = line.split('(')

            if(box == "box"):

                (xpos, ypos, xlong, ylong, trash) = info.split(',')

                xpos = float(xpos)
                ypos = float(ypos)
                xlong = float(xlong)
                ylong = float(ylong)

                xlo = (xpos - xlong / 2)
                xhi = (xpos + xlong / 2)

                ylo = (ypos - ylong / 2)
                yhi = (ypos + ylong / 2)

                xlo = int(xlo)
                xhi = int(xhi)

                ylo = int(ylo)
                yhi = int(yhi)

                if (xlo < 1):

                    xlo = 1

                if (xhi > ncol):

                    xhi = ncol

                if (ylo < 1):

                    ylo = 1

                if (yhi > nrow):

                    yhi = nrow

                img[ylo - 1:yhi, xlo - 1:xhi] = val



    if os.path.isfile(maskimage):
        print("{} deleted; a new one is created \n".format(maskimage))
        runcmd = "rm {}".format(maskimage)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True



############################
#### Deprecated function ####
############################
def MakeSkyMask(maskimage, catfile, scale, offset, regfile):
    "Create a mask image using ellipses for sky computation. Old: use MakeMask instead"
# replaced by MakeMask
# k Check

    check = False
    flagsat = 4  # flag value when object is saturated (or close to)
#    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Masks for sky \n")

    Rsky = scale * ai * kr + offset

    mask = Rsky < 1
    if mask.any():
        Rsky[mask] = 1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = check.CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by user
        # in ds9
        regflag = check.CheckSatReg(xx[idx], yy[idx], regfile, Rsky[
                              idx], theta[idx], e[idx])

        if (checkflag == False) and (regflag == False):

            print("Creating ellipse annuli for object {}  \n".format(n[idx]))
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rsky[idx], theta[idx], e[
                           idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):

            print(
                "Skipping object {}, one or more pixels are saturated \n".format(n[idx]))

    if os.path.isfile(maskimage):
        print("{} deleted; a new one is created \n".format(maskimage))
        runcmd = "rm {}".format(maskimage)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)



    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def SetZero(image, num, xlo, xhi, ylo, yhi):
    "remove object from mask image"
    # put all pixels
    # indicated in file to zero
    # in the image provided


# just in case:
    (ncol, nrow) = GetAxis(image)
    mask = xlo < 1
    if mask.any():
        if isinstance(xlo,np.ndarray):
            xlo[mask] = 1
        else:
            xlo = 1

    mask = xhi > ncol
    if mask.any():
        if isinstance(xhi,np.ndarray):
            xhi[mask] = ncol
        else:
            xhi = ncol

    mask = ylo < 1
    if mask.any():
        if isinstance(ylo,np.ndarray):
            ylo[mask] = 1
        else:
            ylo = 1

    mask = yhi > nrow
    if mask.any():
        if isinstance(xhi,np.ndarray):
            yhi[mask] = nrow
        else:
            yhi = nrow


    hdu = fits.open(image)
    dat = hdu[0].data


    mask = dat[ylo - 1:yhi, xlo - 1:xhi] == num

    if mask.any():
        dat[ylo - 1:yhi, xlo - 1:xhi][mask] = 0


    if os.path.isfile(image):
        print("{} updated; object mask number {} was removed \n".format(image,num))
        runcmd = "rm {}".format(image)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu[0].data = dat
    hdu.writeto(image, overwrite=True)
    hdu.close()

    return True


#not checked
def CountPix(imgfile,xmin,xmax,ymin,ymax,val):
    "Counts number of pixels in image region with value equal to val"

    hdu = fits.open(imgfile)
    img = hdu[0].data[ymin - 1:ymax, xmin - 1:xmax]

    mask = img == val

    pixs=len(img[mask])

    hdu.close()


    return pixs
