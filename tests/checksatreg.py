import os.path
import os
from astropy.io import fits
import numpy as np
import sys
import scipy
import scipy.special
import subprocess as sp



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
#                print("hola")
                next(lines)
                continue
            
#                next(lines)
#               next(lines)


            print(line)
            (box, info) = line.split("(")

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


