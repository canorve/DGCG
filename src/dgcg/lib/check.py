import os.path
import os
from astropy.io import fits
import numpy as np
import sys
import scipy
import scipy.special
import subprocess as sp


#from dgcg.lib.image import GetAxis # call with the whole path in the future


#def CheckKron(err):
    # K check
#    if err != 0 and err != 256:
#        s = "Error when calling function.  Error number {0} \n".format(err)
#        sys.exit(s)
#    elif err == 256:
#        print("Error number 256 when calling function\n ")
#        return True


def IsString(val):
    # k Check
    # Valida lo que se ingrese como una cadena
    val = str(val)
    if val.isalnum():
        return True
    else:
        return False


def IsInteger(val):
    # k  Check
    # Valida lo que se ingrese como un enero
    val = str(val)
    if val.isdigit():
        return True
    else:
        return False


def IsFloat(val):
    # k Check
    # Valida lo que se ingrese como un flotante
    val = str(val)
    val = val.lstrip("+-")
    if not val.find(".") == -1:
        if not val.isalnum():
            (n1, n2) = val.split(".")
        if n1.isdigit() and n2.isdigit():
            return True
        else:
            return False
    else:
        return False


def CheckFlag(val, check):
    "Check for flag contained in val, returns True if found "

    flag = False
    mod = 1
    maxx=128


    while (mod != 0):

        res = int(val / maxx)

        if (maxx == check and res == 1):

            flag = True

        mod = val % maxx

        val = mod
        maxx = maxx / 2

    return flag


def CheckOverlap(xpos, ypos, R, theta, q, xpos2, ypos2, R2, theta2, q2):
    "Check the distance of two ellipses. returns True if they overlap"
# check the distance of
# 2 ellipses. returns 1
# if they overlap.

# k Check

    flag = False

    theta = theta + 90  # converting from GALFIT to Sextractor position
    theta2 = theta2 + 90

    bim = q * R
    bim2 = q2 * R2

    theta = theta * np.pi / 180  # Rads!!!
    theta2 = theta2 * np.pi / 180  # Rads!!!
    dx = xpos2 - xpos
    dy = ypos2 - ypos

    dx2 = xpos - xpos2
    dy2 = ypos - ypos2

    dist = np.sqrt(dx**2 + dy**2)

    landa = np.arctan2(dy, dx)

    if landa < 0:
        landa = landa + 2 * np.pi

    landa2 = np.arctan2(dy, dx)

    if landa2 < 0:
        landa2 = landa2 + 2 * np.pi

    landa = landa - theta
    landa2 = landa2 - theta2

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)
    angle2 = np.arctan2(np.sin(landa2) / bim2, np.cos(landa2) / R2)

    xell = xpos + R * np.cos(angle) * np.cos(theta) - \
        bim * np.sin(angle) * np.sin(theta)
    yell = ypos + R * np.cos(angle) * np.sin(theta) + \
        bim * np.sin(angle) * np.cos(theta)

    xell2 = xpos2 + R2 * np.cos(angle2) * np.cos(theta2) - \
        bim2 * np.sin(angle2) * np.sin(theta2)
    yell2 = ypos2 + R2 * np.cos(angle2) * np.sin(theta2) + \
        bim2 * np.sin(angle2) * np.cos(theta2)

    dell = np.sqrt((xell - xpos)**2 + (yell - ypos)**2)
    dell2 = np.sqrt((xell2 - xpos2)**2 + (yell2 - ypos2)**2)

    distell = dell + dell2

    if dist <= distell:
        flag = True

    return flag



def CheckSatReg(x, y, filein, R, theta, ell):
    "Check if object is inside of saturated region. returns True if at least one pixel is inside"

    q = (1 - ell)

    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    flag = False
#    fileflag = 1

    with open(filein) as f_in:

        # All lines including the blank ones
        lines = (line.rstrip() for line in f_in)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines


        for line in lines:


            if (line.find("(") == -1):
                next(lines)
                continue

            (box, info) = line.split('(')

            if(box == "box"):

                (xpos, ypos, xlong, ylong, trash) = info.split(',')

                xpos = float(xpos)
                ypos = float(ypos)
                xlong = float(xlong)
                ylong = float(ylong)

                xlo = xpos - xlong / 2
                xhi = xpos + xlong / 2

                ylo = ypos - ylong / 2
                yhi = ypos + ylong / 2

## center
                dx = xpos - x
                dy = ypos - y

                landa = np.arctan2(dy, dx)

                if landa < 0:
                    landa = landa + 2 * np.pi

                landa = landa - theta

                angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                xell = x + R * np.cos(angle) * np.cos(theta) - \
                    bim * np.sin(angle) * np.sin(theta)
                yell = y + R * np.cos(angle) * np.sin(theta) + \
                    bim * np.sin(angle) * np.cos(theta)

# top left

                dx = xlo - x
                dy = yhi - y

                landa = np.arctan2(dy, dx)

                if landa < 0:
                    landa = landa + 2 * np.pi

                landa = landa - theta

                angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                xelltl = x + R * np.cos(angle) * np.cos(theta) - \
                    bim * np.sin(angle) * np.sin(theta)
                yelltl = y + R * np.cos(angle) * np.sin(theta) + \
                    bim * np.sin(angle) * np.cos(theta)

# top right

                dx = xhi - x
                dy = yhi - y

                landa = np.arctan2(dy, dx)

                if landa < 0:
                    landa = landa + 2 * np.pi

                landa = landa - theta

                angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                xelltr = x + R * np.cos(angle) * np.cos(theta) - \
                    bim * np.sin(angle) * np.sin(theta)
                yelltr = y + R * np.cos(angle) * np.sin(theta) + \
                    bim * np.sin(angle) * np.cos(theta)

# bottom left
                dx = xlo - x
                dy = ylo - y

                landa = np.arctan2(dy, dx)

                if landa < 0:
                    landa = landa + 2 * np.pi

                landa = landa - theta

                angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                xellbl = x + R * np.cos(angle) * np.cos(theta) - \
                    bim * np.sin(angle) * np.sin(theta)
                yellbl = y + R * np.cos(angle) * np.sin(theta) + \
                    bim * np.sin(angle) * np.cos(theta)

# bottom right

                dx = xhi - x
                dy = yhi - y

                landa = np.arctan2(dy, dx)

                if landa < 0:
                    landa = landa + 2 * np.pi

                landa = landa - theta

                angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                xellbr = x + R * np.cos(angle) * np.cos(theta) - \
                    bim * np.sin(angle) * np.sin(theta)
                yellbr = y + R * np.cos(angle) * np.sin(theta) + \
                    bim * np.sin(angle) * np.cos(theta)


                if ((xell > xlo and xell < xhi) and (yell > ylo and yell < yhi)):
                    flag = True
                    break

                if ((xelltl > xlo and xelltl < xhi) and (yelltl > ylo and yelltl < yhi)):
                    flag = True
                    break

                if ((xelltr > xlo and xelltr < xhi) and (yelltr > ylo and yelltr < yhi)):
                    flag = True
                    break

                if ((xellbl > xlo and xellbl < xhi) and (yellbl > ylo and yellbl < yhi)):
                    flag = True
                    break

                if ((xellbr > xlo and xellbr < xhi) and (yellbr > ylo and yellbr < yhi)):
                    flag = True
                    break



    return flag



def CheckSatReg2(x,y,filein):
   "Check if object is inside of saturated region. returns True if at least one pixel is inside"
## check if object is inside of
## saturaded region as indicated by ds9 box region
## returns True if object center is in saturaded region

   flag = False

   with open(filein) as f_in:

       lines = (line.rstrip() for line in f_in) # All lines including the blank ones
       lines = (line.split('#', 1)[0] for line in lines) # remove comments
       lines = (line.rstrip() for line in lines)   # remove lines containing only comments
       lines = (line for line in lines if line) # Non-blank lines

       for line in lines:

           if (line != "image"):

               (box,info)=line.split('(')

               if(box == "box"):

                   (xpos,ypos,xlong,ylong,trash)=info.split(',')

                   xpos=float(xpos)
                   ypos=float(ypos)
                   xlong=float(xlong)
                   ylong=float(ylong)


                   xlo = xpos - xlong/2
                   xhi = xpos + xlong/2

                   ylo = ypos - ylong/2
                   yhi = ypos + ylong/2

                   if ( (x > xlo and x < xhi) and (y > ylo and y < yhi) ):
                       flag=True
                       break


   return flag








def CheckKron(xpos, ypos, x, y, R, theta, q):
    "check if position is inside of the Kron Ellipse saturaded region returns True if object center is in Ellipse region"

    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    flag = False

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    if landa < 0:
        landa = landa + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    if dist < dell:
        flag = True

    return flag




def CheckSaneValues(parvar):
    "Check for sane values in parameters file"


    errmsg="Can't found {}; exiting... \n".format(parvar.Img)
    assert os.path.isfile(parvar.Img), errmsg
    (sizex, sizey) = GetAxis(parvar.Img)


    errmsg="Can't found {}; exiting... \n".format(parvar.SexCat)
    assert os.path.isfile(parvar.SexCat), errmsg


    if not os.path.isfile(parvar.SigImg) and not parvar.SigImg != "none":
        print("Can't found {}; but thats OK... \n".format(parvar.SigImg))

#    errmsg="Can't found {}; exiting... \n".format(parvar.SigImg)
#    assert os.path.isfile(parvar.SigImg), errmsg

    if not os.path.isdir(parvar.PsfDir):
        print("Can't found {}; but thats OK...\n".format(parvar.PsfDir))

#    errmsg="Can't found {}; exiting... \n".format(parvar.PsfDir)
#    assert os.path.isdir(parvar.PsfDir), errmsg

    if not (IsFloat(parvar.MagZpt) or IsInteger(parvar.MagZpt)):
        print("MagZpt = {} is not a number; exiting...\n".format(parvar.MagZpt))
        sys.exit()

    if not (IsFloat(parvar.PlateScale) or IsInteger(parvar.PlateScale)):
        print("PlateScale = {} is not a number; exiting...\n".format(parvar.PlateScale))
        sys.exit()

    if ((parvar.FitFunc != "BD") and (parvar.FitFunc != "sersic")):
        print("FitFunc must be 'BD' or 'sersic' only \n")
        sys.exit()

#    if not (IsFloat(parvar.GalBot) or IsInteger(parvar.GalBot)):
#        print("GalBot = {} is not a number; exiting...\n".format(parvar.GalBot))
#        sys.exit()

    if not (IsFloat(parvar.GalClas) or IsInteger(parvar.GalClas)):
        print("GalClas = {} is not a number; exiting...\n".format(parvar.GalClas))
        sys.exit()

    if not (IsInteger(parvar.ConvBox)):
        print("ConvBox = {} is not a integer; exiting...\n".format(parvar.ConvBox))
        sys.exit()


    if not (IsFloat(parvar.FitBox) or IsInteger(parvar.FitBox)):
        print("FitBox = {} is not a number; exiting...\n".format(parvar.FitBox))
        sys.exit()


    if not (IsFloat(parvar.MagDiff) or IsInteger(parvar.MagDiff)):
        print("MagDiff = {} is not a number; exiting...\n".format(parvar.MagDiff))
        sys.exit()


    if not (IsFloat(parvar.KronScale) or IsInteger(parvar.KronScale)):
        print("KronScale = {} is not a number; exiting...\n".format(parvar.KronScale))
        sys.exit()



    if not (IsFloat(parvar.SkyScale) or IsInteger(parvar.SkyScale)):
        print("SkyScale = {} is not a number; exiting...\n".format(parvar.SkyScale))
        sys.exit()

    if not (IsInteger(parvar.Offset)):
        print("Offset = {} is not a integer; exiting...\n".format(parvar.Offset))
        sys.exit()

    if not (IsInteger(parvar.SkyWidth)):
        print("SkyWidth = {} is not a integer; exiting...\n".format(parvar.SkyWidth))
        sys.exit()

    if not (IsFloat(parvar.NSer) or IsInteger(parvar.NSer)):
        print("NSer = {} is not a number; exiting...\n".format(parvar.NSer))
        sys.exit()

    if not (IsInteger(parvar.MaxFit)):
        print("MaxFit = {} is not a integer; exiting...\n".format(parvar.MaxFit))
        sys.exit()

#    if not (IsFloat(parvar.MagMin) or IsInteger(parvar.MagMin)):
#        print("MagMin = {} is not a number; exiting...\n".format(parvar.MagMin))
#        sys.exit()

    if not (IsFloat(parvar.MagMax) or IsInteger(parvar.MagMax)):
        print("MagMax = {} is not a number; exiting...\n".format(parvar.MagMax))
        sys.exit()

    if not (IsInteger(parvar.FlagSex)):
        print("FlagSex = {} is not a integer; exiting...\n".format(parvar.FlagSex))
        sys.exit()

    if not os.path.isfile(parvar.ConsFile):
        print("Can't found {}; but that's OK... \n".format(parvar.ConsFile))

    if not (parvar.Region == True or parvar.Region == False):
        print("Region must be a boolean; exiting...\n")
        sys.exit()

    if not (IsInteger(parvar.Bxmin)):
        print("Bxmin = {} is not a integer; exiting...\n".format(parvar.Bxmin))
        sys.exit()
    elif (parvar.Bxmin < 1 or parvar.Bxmin > sizex) and (parvar.Region == True):
        print("Boundary number is outside image; exiting...\n")
        sys.exit()
    if not (IsInteger(parvar.Bxmax)):
        print("Bxmax = {} is not a integer; exiting...\n".format(parvar.Bxmax))
        sys.exit()
    elif (parvar.Bxmax < 1 or parvar.Bxmax > sizex) and (parvar.Region == True):
        print("Boundary number is outside image; exiting...\n")
        sys.exit()

    if not (IsInteger(parvar.Bymin)):
        print("Bymin = {} is not a integer; exiting...\n".format(parvar.Bymin))
        sys.exit()
    elif (parvar.Bymin < 1 or parvar.Bymin > sizey) and (parvar.Region == True):
        print("Boundary number is outside image; exiting...\n")
        sys.exit()

    if not (IsInteger(parvar.Bymax)):
        print("Bymax = {} is not a integer; exiting...\n".format(parvar.Bymax))
        sys.exit()
    elif (parvar.Bymax < 1 or parvar.Bymax > sizey) and (parvar.Region == True):
        print("Boundary number is outside image; exiting...\n")
        sys.exit()

    if not (IsInteger(parvar.Split)):
        print("Split = {} is not a integer; exiting...\n".format(parvar.Split))
        sys.exit()

#    if not (parvar.AutoSatRegion == True or parvar.AutoSatRegion == False):
#        print("AutoSatRegion must be a boolean; exiting...\n")
#        sys.exit()

    if not (IsFloat(parvar.SatRegionScale) or IsInteger(parvar.SatRegionScale)):
        print("SatRegionScale = {} is not a number; exiting...\n".format(parvar.SatRegionScale))
        sys.exit()

    if not os.path.isfile(parvar.Ds9SatReg) and not parvar.Ds9SatReg != "none":
        print("Can't found {}; but thats OK... \n".format(parvar.Ds9SatReg))

    if not (IsInteger(parvar.Ds9OutNum)):
        print("Ds9OutNum = {} is not a integer; exiting...\n".format(parvar.Ds9OutNum))
        sys.exit()

    if not (parvar.SimulFit == True or parvar.SimulFit == False):
        print("SimulFit must be a boolean; exiting...\n")
        sys.exit()

    if not (parvar.Erase == True or parvar.Erase == False):
        print("Erase must be a boolean; exiting...\n")
        sys.exit()

    if not (parvar.Nice == True or parvar.Nice == False):
        print("Nice must be a boolean; exiting...\n")
        sys.exit()

    if not (parvar.Overwrite == True or parvar.Overwrite == False):
        print("Overwrite must be a boolean; exiting...\n")
        sys.exit()

    if not (parvar.Execute == 0 or parvar.Execute == 1 or parvar.Execute == 2 or parvar.Execute == 3):
        print("Execute must be 0, 1, 2, 3 only; exiting...\n")
        sys.exit()

    if parvar.LogFile == "fit.log":
        print("LogFile (FileOut) can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.Ds9SatReg == "fit.log":
        print("Ds9SatReg can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.SegFile == "fit.log":
        print("SegFile can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.SkyFile == "fit.log":
        print("SkyFile can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.Ds9OutName == "fit.log":
        print("Ds9OutName can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.Ds9FitReg == "fit.log":
        print("Ds9FitReg can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.BoxOut == "fit.log":
        print("BoxOut can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.BoxSkyOut == "fit.log":
        print("BoxSkyOut can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.SexSort == "fit.log":
        print("SexSort can't be named fit.log; exiting...\n")
        sys.exit()

    if parvar.SexArSort == "fit.log":
        print("SexArSort can't be named fit.log; exiting...\n")
        sys.exit()

# Check if Sextractor catalog has sane values

    sexflag = False

    try:
        N, Alpha, Delta, X, Y, Mg, Kr, Fluxr, Isoa, Ai, E, Theta, Bkgd, Idx, Flg = np.genfromtxt(
            parvar.SexCat, delimiter="", unpack=True)

        N = N.astype(int)
        Flg = Flg.astype(int)

    except:
        print("Unexpected error at reading Sex file:", sys.exc_info()[0])
        raise

    idran = np.random.randint(0, len(N) - 1)

    if not (IsInteger(N[idran])):
        sexflag = True
        print("Error in column 1 \n")

    if not (IsFloat(Alpha[idran]) or IsInteger(Alpha[idran])):
        sexflag = True
        print("Error in column 2 \n")

    if not (IsFloat(Delta[idran]) or IsInteger(Delta[idran])):
        sexflag = True
        print("Error in column 3 \n")

    if not (IsFloat(X[idran]) or IsInteger(X[idran])):
        sexflag = True
        print("Error in column 4 \n")

    if not (IsFloat(Y[idran]) or IsInteger(Y[idran])):
        sexflag = True
        print("Error in column 5 \n")

    if not (IsFloat(Mg[idran]) or IsInteger(Mg[idran])):
        sexflag = True
        print("Error in column 6 \n")

    if not (IsFloat(Kr[idran]) or IsInteger(Kr[idran])):
        sexflag = True
        print("Error in column 7 \n")

    if not (IsFloat(Fluxr[idran]) or IsInteger(Fluxr[idran])):
        sexflag = True
        print("Error in column 8 \n")

    if not (IsFloat(Isoa[idran]) or IsInteger(Isoa[idran])):
        sexflag = True
        print("Error in column 9 \n")

    if not (IsFloat(Ai[idran]) or IsInteger(Ai[idran])):
        sexflag = True
        print("Error in column 10 \n")

    if not (IsFloat(E[idran]) or IsInteger(E[idran])):
        sexflag = True
        print("Error in column 11 \n")
    elif E[idran] > 1 or E[idran] < 0:
        sexflag = True
        print("Error in column 11 \n")

    if not (IsFloat(Theta[idran]) or IsInteger(Theta[idran])):
        sexflag = True
        print("Error in column 12 \n")

    if not (IsFloat(Bkgd[idran]) or IsInteger(Bkgd[idran])):
        sexflag = True
        print("Error in column 13 \n")

    if not (IsFloat(Idx[idran]) or IsInteger(Idx[idran])):
        sexflag = True
        print("Error in column 14 \n")

    elif Idx[idran] > 1 or Idx[idran] < 0:
        sexflag = True
        print("Error in column 14 \n")

    if (not IsInteger(Flg[idran])):
        sexflag = True
        print("Error in column 15 \n")

    if sexflag:

        print("Something is wrong in the Sextractor catalog please check: \n")
        print("  1 NUMBER                 Running object number             \n")
        print(
            "  2 ALPHA_J2000            Right ascension of barycenter (J2000)   [deg]   \n")
        print(
            "  3 DELTA_J2000            Declination of barycenter (J2000)       [deg]   \n")
        print(
            "  4 X_IMAGE                Object position along x                [pixel] \n")
        print(
            "  5 Y_IMAGE                Object position along y                 [pixel] \n")
        print(
            "  6 MAG_APER               Fixed aperture magnitude vector         [mag]   \n")
        print("  7 KRON_RADIUS            Kron apertures in units of A or B                 \n")
        print(
            "  8 FLUX_RADIUS            Fraction-of-light radii                 [pixel] \n")
        print(
            "  9 ISOAREA_IMAGE          Isophotal area above Analysis threshold [pixel**2]\n")
        print(
            "  10 A_IMAGE               Profile RMS along major axis            [pixel] \n")
        print("  11 ELLIPTICITY           1 - B_IMAGE/A_IMAGE                              \n")
        print(
            "  12 THETA_IMAGE           Position angle (CCW/x)                  [deg]   \n")
        print(
            "  13 BACKGROUND            Background at centroid position         [count] \n")
        print("  14 CLASS_STAR            S/G classifier output                            \n")
        print("  15 FLAGS                 Extraction flags                                  \n")
        print("       \n")

        print("Sextractor Catalog don't have 15 columns or some value is not a sane value   \n ")

        sys.exit()


    try:  # is GALFIT installed?
    # pipe output to /dev/null for silence
        null = open("/dev/null", "w")
        sp.Popen("galfit", stdout=null, stderr=null)
        null.close()
    except OSError:
        print("galfit not found")
        sys.exit()
