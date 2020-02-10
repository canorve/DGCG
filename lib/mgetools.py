#! /usr/bin/env python3


import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path
import scipy
import scipy.special
import matplotlib.pyplot as plt
import mimetypes

from mgefit.sectors_photometry import sectors_photometry
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,NullFormatter,
                               AutoMinorLocator,LogLocator,LinearLocator,AutoLocator)

def EllipSec(galfile,qarg,parg,flagsub,flagpix,ranx,rany,dpi,flagout):


    #class for saving parameters
    params=InputParams()
    #
    params.galfile= galfile

    params.flaglogx=False

    params.flagq=True
    params.qarg=qarg

    params.flagpa=True
    params.parg=parg

    params.flagsub=flagsub

    params.flagpix=flagpix

    params.flagranx[0]=True
    params.ranx=ranx

    params.flagrany[0]=True
    params.rany=rany
    
    params.flagrid=False

    params.flagdpi = True
    params.dpival = dpi 

    params.flagnoplot=True
    #params.dplot=dplot

    params.flagout=flagout

    ################## search arguments after the option:
#    if params.flagpa == True:
 #       opt={}
  #      OptionHandle="--pa"
  #      opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
   #     params.parg=np.int(opt['pa'])

 #   if params.flagq == True:
 #       opt={}
 #       OptionHandle="--q"
 #       opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
 #       params.qarg=np.float(opt['q'])

 #   if params.flagranx[0] == True:
 #       opt={}
 #       OptionHandle="--ranx"
 #       opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]

  #      params.rangex=opt["ranx"]
   #     if "-" in params.rangex:
    #        params.flagranx[1] = True
   #         params.ranx=opt['ranx']
    #    else:
   #         params.ranx=np.float(opt['ranx'])

    #if params.flagrany[0]== True:
     #   opt={}
     #   OptionHandle="--rany"
     #   opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]

      #  params.rangey=opt["rany"]
      #  if "-" in params.rangey:
       #     params.flagrany[1] = True
    #        params.rany=opt['rany']
     #   else:
      #      params.rany=np.float(opt['rany'])

#    if params.flagdpi == True:
 #       opt={}
 #       OptionHandle="--dpi"
  #      opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
   #     params.dpival=np.int(opt['dpi'])

    if params.flagnoplot == True:
       params.dplot=False

  #  if params.flagout == True:
  #      opt={}
  #      OptionHandle="--out"
  #      opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
  #      params.output=opt['out']

#    print("angle in multi-plot is measured from the galaxy's major axis ")

    print("analizing GALFIT file: ",params.galfile)
    ########  End of parameter reading #############
    ################################################
    ################################################

    #class for GALFIT's parameters
    galpar=GalfitParams()



    ######################################
    ####### Read Galfit File #############
    #  xc,yc,q,ang,skylevel,scale,file,mgzpt,exptime,mask=ReadGALFITout(params.galfile,galpars)
    ReadGALFITout(params.galfile,galpar)



######################################
######################################

    if params.flagq == True:
        galpar.q=params.qarg

    if params.flagpa == True:
        galpar.ang=params.parg


    str = "q = {} is used ".format(galpar.q)
    print(str)

    str = "pa = {} is used ".format(galpar.ang)
    print(str)

    ##
    str = "dpi = {} for plots ".format(params.dpival)
    print(str)
    ##


    #small modification for dgcg
    #A168-011526.64-p001221.1-410-out.fits
    (tmp)=galpar.outimage.split("-")

    params.namefile=tmp[0] + "-" + tmp[1] +"-"+ tmp[2] +"-"+ tmp[3] 

    # names for the different png

    params.namepng=params.namefile + ".png"
    params.namesec=params.namefile + "-gal.png"
    params.namemod=params.namefile + "-mod.png"
    params.namemul=params.namefile + "-mul.png"
    params.namesub=params.namefile + "-sub.fits"

    params.output=params.namefile + "-sbout.txt"

    if params.flagout == True: 
        msg="surface brightness output file: {} ".format(params.output)
        print(msg)


    # hdu 1 => image   hdu 2 => model

    errmsg="file {} does not exist".format(galpar.outimage)

    assert os.path.isfile(galpar.outimage), errmsg

    hdu = fits.open(galpar.outimage)
    galpar.img = hdu[1].data
    galpar.model = hdu[2].data
    hdu.close()



    #   numsectors=19
    numsectors=15
    minlevel=-100  # minimun value for sky


    limx,limy=EllipSectors(galpar, params, n_sectors=numsectors, minlevel=minlevel)



    params.Comps,params.N=ReadNComp(params.galfile,galpar.xc,galpar.yc)

    print("Number of components = ",params.N)


##############################################
##############################################
##############################################

    if params.dplot:
        plt.pause(1.5)
    plt.savefig(params.namepng,dpi=params.dpival)
    plt.close()



########################################################
################ Multiplots: ###########################
########################################################
    MulEllipSectors(galpar, params, n_sectors=numsectors,  minlevel=minlevel)


    if params.dplot:
        plt.pause(1.5)

    plt.savefig(params.namemul,dpi=params.dpival)
    plt.close()


    if galpar.mask != None:
        os.remove(galpar.mask) # removing temp mask file

##############       #############
##############  END  #############
##############       #############


#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/


### class for parameters
class InputParams:

    flaglogx=False
    flagq=False
    flagpa=False
    flagsub=False
    flagpix=False
    flagranx=[False,False]
    flagrany=[False,False]
    flagnoplot=False
    flagrid=False
    flagdpi=False
    flagout=False
    #init
    qarg=1
    parg=0
    ranx=1
    rany=1
    dplot=True

    dpival=100

    # init sub values
    Comps=np.array([False])
    N=0


    #input file
    galfile= "galfit.01"

    #output file
    output = "sbout.txt"

    namefile="none"
    namepng="none.png"
    namesec="none-gal.png"
    namemod="none-mod.png"
    namemul="none-mul.png"
    namesub="none-sub.fits"



### class for Galfit parameters
class GalfitParams:

    xc=1
    yc=1
    q=1
    ang=0
    skylevel=0
    scale=1
    inputimage="galaxy.fits"
    outimage="galaxy-out.fits"
    maskimage="galaxy-mask.fits"
    mgzpt=25
    exptime=1
    mask="mask.fits"
    xmin=1
    xmax=2
    ymin=1
    ymax=2
    img = np.array([[1,1],[1,1]])
    model = np.array([[1,1],[1,1]])

##### end of classes

def EllipSectors(galpar, params, n_sectors=19, minlevel=0):

    badpixels=galpar.mask

    # removing background:
    galpar.img = galpar.img - galpar.skylevel
    galpar.model = galpar.model - galpar.skylevel

    xradm = []
    ysbm = []
    ysberrm = []

    if badpixels is not None:

        errmsg="file {} does not exist".format(badpixels)
        assert os.path.isfile(badpixels), errmsg

        hdu = fits.open(badpixels)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()
    else:
        maskb=None

    eps=1-galpar.q

    if params.dplot:
        plt.clf()
        print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp=galpar.xc
    xctemp=galpar.yc
    # and angle is different as well:
    angsec=90-galpar.ang
    #    angsec=ang


    ###############################
    #  galaxy:

    g = sectors_photometry(galpar.img, eps, angsec, xctemp, yctemp, minlevel=minlevel,
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namesec)


    ###################################################

    stidxg = np.argsort(g.radius)

    mgerad=g.radius[stidxg]
    mgecount=g.counts[stidxg]
    mgeangle=g.angle[stidxg]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

    #changing to arc sec
    aellarcg=aellabg*galpar.scale


    ####
    # formula according to cappellary mge manual:
    mgesbg= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ####

    stidxq = np.argsort(aellarcg)


    xarcg = aellarcg[stidxq]
    ymgeg = mgesbg[stidxq]

    #############  Function to order SB along X-axis

    xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)

    ################


    ###############################
    #  model:
    m = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=minlevel,
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namemod)
    ###################################################

    stidxm = np.argsort(m.radius)

    mgerad=m.radius[stidxm]
    mgecount=m.counts[stidxm]
    mgeangle=m.angle[stidxm]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabm= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarcm=aellabm*galpar.scale

    # formula according to cappellary mge manual
    mgesbm= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ##

    stidxq = np.argsort(aellarcm)

    xarcm = aellarcm[stidxq]
    ymgem = mgesbm[stidxq]

    ######  Function to order SB along X-axis

    xradm, ysbm, ysberrm    = FindSB(xarcm, ymgem, n_sectors)

    ################ Plotting

    limx,limy,axsec=PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,galpar.scale)

    ###

    if params.flagout == True: 

        OUTFH = open (params.output,"w")

        lineout= "#        sectors_photometry used with pa={}         q={}       \n".format(galpar.ang,galpar.q)
        OUTFH.write(lineout)

        lineout= "#            Galaxy                          Model                \n"
        OUTFH.write(lineout)

        lineout= "#     rad      SB        SBerr       rad       SB        SBerr \n"
        OUTFH.write(lineout)

        lineout= "# (arcsec) (mag/arcsec) (error)   (arcsec) (mag/arcsec)  (error) \n"
        OUTFH.write(lineout)

        for idx, item in enumerate(xradq):
            if idx < len(xradm):
                lineout= "{0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} \n".format(xradq[idx],ysbq[idx],ysberrq[idx],xradm[idx],ysbm[idx],ysberrm[idx])
            else:
                lineout= "{0:.3f} {1:.3f} {2:.3f} \n".format(xradq[idx],ysbq[idx],ysberrq[idx])
            OUTFH.write(lineout)
        OUTFH.close()



    #### Creating Subcomponents images with Galfit

    params.Comps=[]

    if params.flagsub:

        print("running galfit to create subcomponents...")

        rungal = "galfit -o3 {}".format(params.galfile)
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE,
        stderr=sp.PIPE, universal_newlines=True)

        runchg = "mv subcomps.fits {}".format(params.namesub)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE,
        stderr=sp.PIPE, universal_newlines=True)


        # read number of components
        params.Comps,params.N=ReadNComp(params.galfile,galpar.xc,galpar.yc)

        xradq,ysbq,n=SubComp(params.namesub,params.N,params.Comps,galpar.mgzpt,galpar.exptime,galpar.scale,galpar.xc,galpar.yc,galpar.q,galpar.ang,params.flagpix,axsec,skylevel=galpar.skylevel,
              n_sectors=n_sectors, badpixels=badpixels, minlevel=minlevel)



    axsec.legend(loc=1)

    return limx,limy


def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,scale):
    """
    Produces final best-fitting plot

    """

    # subplot for arc sec axis
    fig, axsec = plt.subplots()


    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)


    minrad = np.min(xradq)
    if params.flagranx[1] == False:
        maxrad = np.max(xradq) * params.ranx
    else:
        maxrad = np.max(xradq)

    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt/mincnt)**np.array([-0.05, +1.05])

    if params.flagrany[1] == False:
        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran2 #inverted axis
        yran[1] = yran1


    axsec.set_xlabel("radius ('')")
    axsec.set_ylabel("mag/''")

    axsec.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy")
    axsec.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model")
    if params.flagranx[1] == False:
        axsec.set_xlim(xran)
    else:
        axsec.set_xlim(xmin,xmax)


    if params.flagrany[1] == False:
        axsec.set_ylim(yran)
    else:
        axsec.set_ylim(ymax,ymin) #inverted


    if params.flaglogx == True:

        axsec.set_xscale("log")

        locmaj = LogLocator(base=10,numticks=12)
        axsec.xaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        axsec.xaxis.set_minor_locator(locmin)
        axsec.xaxis.set_minor_formatter(NullFormatter())

    else:
        axsec.xaxis.set_minor_locator(AutoMinorLocator())
        axsec.xaxis.set_major_locator(AutoLocator())


    axsec.tick_params(which='both', width=2)
    axsec.tick_params(which='major', length=7)
    axsec.tick_params(which='minor', length=4, color='r')

    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())

    if params.flagpix == True:

        axpix = axsec.twiny()

        axpix.set_xlabel("(pixels)")
        x1, x2 = axsec.get_xlim()

        axpix.set_xlim(x1/scale, x2/scale)

        axpix.figure.canvas.draw()


        if params.flaglogx == True:
            axpix.set_xscale("log")
            locmaj = LogLocator(base=10,numticks=12)
            axpix.xaxis.set_major_locator(locmaj)

            locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
            axpix.xaxis.set_minor_locator(locmin)
            axpix.xaxis.set_minor_formatter(NullFormatter())
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())
            axpix.xaxis.set_major_locator(AutoLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        axret=axsec
    else:
        axret=axsec

    if params.flagrid == True:
        # Customize the major grid
        axsec.grid(which='major', linestyle='-', linewidth='0.7', color='black')
        # Customize the minor grid
        axsec.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

    return xran,yran,axret



def SubComp(namesub,N,Comps,mgzpt,exptime,scale,xc,yc,q,ang,flagpix,axsec,skylevel=0,
    n_sectors=19, badpixels=None, minlevel=0):

    errmsg="file {} does not exist".format(namesub)
    assert os.path.isfile(namesub), errmsg

    hdu = fits.open(namesub)

    subimgs=[]


    cnt=0  # image =0 do not count
    while(cnt<len(Comps)):

        if Comps[cnt] == True:
            img = hdu[cnt+2].data
            subimgs.append(img)

        cnt=cnt+1

    hdu.close()

    xradm = []
    ysbm = []
    ysberrm = []

    if badpixels is not None:

        errmsg="file {} does not exist".format(badpixels)
        assert os.path.isfile(badpixels), errmsg

        hdu = fits.open(badpixels)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()
    else:
        maskb=None

    eps=1-q


    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp=xc
    xctemp=yc
    # and angle is different as well:
    angsec=90-ang

    ####################

    ab=q
    n=0
    while(n<N):

        subim=subimgs[n]

        scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,
                plot=0, badpixels=maskb, n_sectors=n_sectors)

        ###################################################

        stidxg = np.argsort(scmp.radius)

        mgerad=scmp.radius[stidxg]
        mgecount=scmp.counts[stidxg]
        mgeangle=scmp.angle[stidxg]
        mgeanrad=np.deg2rad(mgeangle)


        aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

        #converting to arcsec

        aellarcg=aellabg*scale


        # formula according to cappellary mge manual

        mgesbg= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(scale**2) + 0.1


        stidxq = np.argsort(aellarcg)

        xarcg = aellarcg[stidxq]
        ymgeg = mgesbg[stidxq]

        xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)


        PlotSub(xradq,ysbq,n,axsec)

        n=n+1


    return  xradq,ysbq,n


def PlotSub(xradq,ysbq,nsub,axsec):
    """
    Produces subcomponent plot

    """

    substr="component "+np.str(nsub+1)

    axsec.plot(xradq, ysbq,'--',color='skyblue',linewidth=4,markersize=0.7,label=substr)


def MulEllipSectors(galpar, params, n_sectors=19, minlevel=0):

    badpixels=galpar.mask

    #it's already done
    #galpar.img = galpar.img - galpar.skylevel
    #galpar.model = galpar.model - galpar.skylevel

    fignum=1

    if badpixels is not None:

        errmsg="file {} does not exist".format(badpixels)
        assert os.path.isfile(badpixels), errmsg

        hdu = fits.open(badpixels)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()
    else:
        maskb=None


    eps=1-galpar.q

    if params.dplot:
        plt.clf()
        print("")

    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)

    ###################
    # I have to switch x and y values because they are different axes for
    # numpy
    xtemp=galpar.xc
    galpar.xc=galpar.yc
    galpar.yc=xtemp

    angsec=90-galpar.ang
    ######################

    sg = sectors_photometry(galpar.img, eps, angsec, galpar.xc, galpar.yc,minlevel=minlevel,
            plot=False, badpixels=maskb, n_sectors=n_sectors)


    sm = sectors_photometry(galpar.model, eps, angsec, galpar.xc, galpar.yc,minlevel=minlevel,
            plot=False, badpixels=maskb, n_sectors=n_sectors)

    ###################################################

    stidx = np.argsort(sg.radius)

    #   galaxy
    mgerad=sg.radius[stidx]

    mgecount=sg.counts[stidx]
    mgeangle=sg.angle[stidx]
    mgeanrad=np.deg2rad(mgeangle)


    # model

    stidx = np.argsort(sm.radius)

    mgemodrad=sm.radius[stidx]

    mgemodcount=sm.counts[stidx]
    mgemodangle=sm.angle[stidx]
    mgemodanrad=np.deg2rad(mgemodangle)


    # converting to pixels

    mgerad=mgerad*galpar.scale
    mgemodrad=mgemodrad*galpar.scale


    # formula according to cappellary mge manual
    # galaxy:
    mgesb= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    # Model:
    mgemodsb= galpar.mgzpt - 2.5*np.log10(mgemodcount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1


    if params.flagsub:

        wtemp=[]
        mgesbsub=[]
        mgeradsub=[]
        mgeanglesub=[]
        rsub=[]
        hdusub = fits.open(params.namesub)

        subimgs=[]

        cnt=0
        while(cnt<len(params.Comps)):

            if params.Comps[cnt] == True:
                imgsub = hdusub[cnt+2].data
                subimgs.append(imgsub)
            cnt=cnt+1

        hdu.close()


    ###############################
    #  galaxy:
        ab=galpar.q
        ni=0
        while(ni<params.N):

            subim=subimgs[ni]

            subcmp = sectors_photometry(subim, eps, angsec, galpar.xc, galpar.yc, minlevel=minlevel,
                    plot=0, badpixels=maskb, n_sectors=n_sectors)

            subidx = np.argsort(subcmp.radius)

            temprad=subcmp.radius[subidx]

            #converting to arcsec
            temprad=temprad*galpar.scale

            mgecountsub=subcmp.counts[subidx]

            tempangle=subcmp.angle[subidx]
            mgeanradsub=np.deg2rad(tempangle)


            # formula according to cappellary mge manual
            tempmge= galpar.mgzpt - 2.5*np.log10(mgecountsub/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1

            mgesbsub.append(tempmge)
            mgeradsub.append(temprad)
            mgeanglesub.append(tempangle)

            ni+=1


    minrad = np.min(mgerad)

    if params.flagranx[1] == False:
        maxrad = np.max(mgerad) * params.ranx
    else:
        maxrad = np.max(mgerad)

    minsb = np.min(mgesb)
    maxsb = np.max(mgesb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = minsb * (maxsb/minsb)**np.array([+1.05,-0.05])


    if params.flagrany[1] == False:

        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran1
        yran[1] = yran2


    sectors = np.unique(mgeangle)
    n = sectors.size
    dn = int(round(n/6.))
    nrows = (n-1)//dn + 1 # integer division

    plt.clf()

    fig, axsec = plt.subplots(nrows, 2, sharex=True, sharey='col', num=fignum)
    fig.subplots_adjust(hspace=0.01)


    if params.flagpix:
        axpix = axsec[0,0].twiny()
        axpix2 = axsec[0,1].twiny()

    fig.text(0.04, 0.5, 'Surface brightness', va='center', rotation='vertical')
    fig.text(0.96, 0.5, 'error (%)', va='center', rotation='vertical')

    axsec[-1, 0].set_xlabel("radius ('')")
    axsec[-1, 1].set_xlabel("radius ('')")

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        if params.flagpix:
            axpix.set_xscale("log")
            axpix2.set_xscale("log")

    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())

    axsec[-1, 0].tick_params(which='both', width=2)
    axsec[-1, 0].tick_params(which='major', length=7)
    axsec[-1, 0].tick_params(which='minor', length=4, color='r')

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())
    axsec[-1, 1].tick_params(which='both', width=2)
    axsec[-1, 1].tick_params(which='major', length=7)
    axsec[-1, 1].tick_params(which='minor', length=4, color='r')


    #    row = 7 # old values
    row = nrows -1
    for j in range(0, n, dn):
        w = np.nonzero(mgeangle == sectors[j])[0]
        w = w[np.argsort(mgerad[w])]
        r = mgerad[w]

        wmod = np.nonzero(mgemodangle == sectors[j])[0]
        wmod = wmod[np.argsort(mgemodrad[wmod])]

        if (len(mgemodrad) < len(mgerad)):
            r2 = mgemodrad[wmod]
        else:
            wmod=w
            r2 = mgemodrad[wmod]


        txtang= sectors[j]
        txt = "$%.f^\circ$" % txtang

        if params.flagranx[1] == False:
            axsec[row, 0].set_xlim(xran)
        else:
            axsec[row, 0].set_xlim(xmin,xmax)

        if params.flagrany[1] == False:
            axsec[row, 0].set_ylim(yran)
        else:
            axsec[row, 0].set_ylim(ymax,ymin) #inverted

        if params.flaglogx == False:
            axsec[row, 0].plot(r, mgesb[w], 'C3o')

            axsec[row, 0].plot(r2, mgemodsb[wmod], 'C0-', linewidth=2)

        else:
            axsec[row, 0].semilogx(r, mgesb[w], 'C3o')

            axsec[row, 0].semilogx(r2, mgemodsb[wmod], 'C0-', linewidth=2)

        if params.flagrid == True:
        # Customize the major grid
            axsec[row,0].grid(which='major', linestyle='-', linewidth='0.7', color='black')
        # Customize the minor grid
            axsec[row,0].grid(which='minor', linestyle=':', linewidth='0.5', color='black')

        #  axsec[row,0].grid(True)

        if params.flagsub == True:
            ii=0
            while(ii<params.N):

                wtemp = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
                wtemp = wtemp[np.argsort(mgeradsub[ii][wtemp])]

                rtemp = mgeradsub[ii][wtemp]


                if params.flaglogx == False:

                    axsec[row, 0].plot(rtemp, mgesbsub[ii][wtemp],'--',color='skyblue', linewidth=2)

                else:

                    axsec[row, 0].semilogx(rtemp, mgesbsub[ii][wtemp], '--',color='skyblue', linewidth=2)


                ii+=1

        axsec[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 0].transAxes)

        if (len(mgemodrad) < len(mgerad)):
            sberr=1-mgemodsb[wmod]/mgesb[wmod]
            axsec[row, 1].plot(r2, sberr*100, 'C0o')
        else:
            sberr=1-mgemodsb[w]/mgesb[w]
            axsec[row, 1].plot(r, sberr*100, 'C0o')


        axsec[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        axsec[row, 1].yaxis.tick_right()
        axsec[row, 1].yaxis.set_label_position("right")
        axsec[row, 1].set_ylim([-19.5, 20])
        # axsec[row, 1].set_ylim([-20, 20])


        if params.flagranx[1] == False:
            axsec[row, 1].set_xlim(xran)
        else:
            axsec[row, 1].set_xlim(xmin,xmax)


        axsec[row, 0].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 0].tick_params(which='both', width=2)
        axsec[row, 0].tick_params(which='major', length=7)
        axsec[row, 0].tick_params(which='minor', length=4, color='r')

        axsec[row, 1].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 1].tick_params(which='both', width=2)
        axsec[row, 1].tick_params(which='major', length=7)
        axsec[row, 1].tick_params(which='minor', length=4, color='r')


        row -= 1

    if params.flagpix == True:
        axpix.set_xlabel("(pixels)")

        ##x1, x2 = axsec[7,0].get_xlim()
                
        if params.flagranx[1] == False:
            x1= xran[0]
            x2= xran[1]
        else:
            x1=xmin
            x2=xmax

        axpix.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix.figure.canvas.draw()

        axpix2.set_xlabel("(pixels)")
        axpix2.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix2.figure.canvas.draw()

        ##
        if params.flaglogx == True:
            axpix.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        if params.flaglogx == True:
            axpix2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix2.xaxis.set_minor_locator(AutoMinorLocator())
        axpix2.tick_params(which='both', width=2)
        axpix2.tick_params(which='major', length=7)
        axpix2.tick_params(which='minor', length=4, color='r')


def FindSB(xarcq, ymgeq, numsectors):
    # the xarcq array must be ordered
    # use mag instead of counts

    xradq=[]
    ysbq=[]
    ysberrq=[]
    xradq=np.array(xradq)
    ysbq=np.array(ysbq)
    ysberrq=np.array(ysberrq)

    numsave=0
    tot=xarcq.size
    count=0
    for i in range(tot,0,-1):

        lima=i-numsectors
        limb=i

        valstd=np.std(xarcq[lima:limb])
        if valstd < 0.1:
            numsave=count
            break
        count=count+1
    init=numsave%numsectors

    n=init

    num=np.int((xarcq.size-init)/numsectors)
    n=xarcq.size-init
    for i in range(num,0,-1):

        lima=n-numsectors
        limb=n

        xradq=np.append(xradq,np.mean(xarcq[lima:limb]))
        ysbq=np.append(ysbq,np.mean(ymgeq[lima:limb]))
        ysberrq=np.append(ysberrq,np.std(ymgeq[lima:limb]))

        n=n-numsectors

    return xradq, ysbq, ysberrq


def ReadGALFITout(inputf,galpar):

    flagfirst = True

    maskimage = ""
    #    skylevel=0

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    while index < len(lines):

        #================================================================================
        # IMAGE and GALFIT CONTROL PARAMETERS
        #A) tempfits/A2399-3-2.fits      # Input data image (FITS file)
        #B) A2399-215615.96-m073822.7-337-out.fits      # Output data image block
        #C) tempfits/none-3-2.fits      # Sigma image name (made from data if blank or "none")
        #D) psfs/PSF-1309-721.fits          # Input PSF image and (optional) diffusion kernel
        #E) 1                   # PSF fine sampling factor relative to data
        #F) mask-337            # Bad pixel mask (FITS image or ASCII coord list)
        #G) constraints         # File with parameter constraints (ASCII file)
        #H) 129  809  265  809  # Image region to fit (xmin xmax ymin ymax)
        #I) 60     60           # Size of the convolution box (x y)
        #J) 21.672              # Magnitude photometric zeropoint
        #K) 0.680  0.680        # Plate scale (dx dy)   [arcsec per pixel]
        #O) regular             # Display type (regular, curses, both)
        #P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        line = lines[index]
        (tmp) = line.split()
        if tmp[0] == "A)":     # input image
            galpar.inputimage=tmp[1]

        if tmp[0] == "B)":     # out image
            galpar.outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            galpar.maskimage=tmp[1]

        if tmp[0] == "H)":     # region fit box
            galpar.xmin=int(tmp[1])
            galpar.xmax=int(tmp[2])
            galpar.ymin=int(tmp[3])
            galpar.ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            galpar.mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            galpar.scale=float(tmp[1])


        # first component
        if tmp[0] == "0)" and flagfirst == True:     # plate scale

            flagfirst=False

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    galpar.xc=float(tmp[1])
                    galpar.yc=float(tmp[2])

                if tmp[0] == "9)":    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)": # position angle
                    galpar.ang=float(tmp[1])

        # sersic component
        if tmp[0] == "0)" and tmp[1] == "sersic":     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcser=float(tmp[1])
                    ycser=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcser)**2+(galpar.yc-ycser)**2)



                if tmp[0] == "9)" and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)" and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])

        # second component exponential model
        if tmp[0] == "0)" and tmp[1] == "expdisk" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcexp=float(tmp[1])
                    ycexp=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcexp)**2+(galpar.yc-ycexp)**2)


                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])


        if tmp[0] == "0)" and tmp[1] == "gaussian":

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgauss=float(tmp[1])
                    ycgauss=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcgauss)**2+(galpar.yc-ycgauss)**2)

                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])



        if tmp[0] == "0)" and tmp[1] == "sky" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":    # axis ratio
                    galpar.skylevel=float(tmp[1])

        index += 1

    GalfitFile.close()


    errmsg="file {} does not exist".format(galpar.inputimage)
    assert os.path.isfile(galpar.inputimage), errmsg

    galpar.exptime=GetExpTime(galpar.inputimage)


    #errmsg="xc and yc are unknown "
    #assert ("xc" in locals()) and ("yc" in locals())  , errmsg

    print("center is at xc, yc = ",galpar.xc,galpar.yc)

    # correcting coordinates
    galpar.xc=galpar.xc-galpar.xmin+1
    galpar.yc=galpar.yc-galpar.ymin+1


    if os.path.isfile(galpar.maskimage):
        mime=mimetypes.guess_type(galpar.maskimage)

        flagbm = not(mime[0] == "text/plain")

        errmsg="Sorry the mask file: {}  must be binary, not ASCII ".format(maskimage)
        assert flagbm is True, errmsg


        galpar.mask="tempmask.fits"
        GetFits(galpar.maskimage, galpar.mask, galpar.xmin, galpar.xmax, galpar.ymin, galpar.ymax)

    else:
        errmsg="Mask file does not exist"
        print(errmsg)
        galpar.mask=None

    #return xc,yc,q,pa,skylevel,scale,outimage,mgzpt,exptime,mask


def ReadNComp(inputf,X,Y):
    ## search and count model components

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    Ncomp=0
    Comps=[]
    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()


        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

            X=X+xmin-1
            Y=Y+ymin-1

        if tmp[0] == "0)":

            if (tmp[1] != "sky"):
                while (tmp[0] != "Z)"):

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":   # center
                        xc=float(tmp[1])
                        yc=float(tmp[2])

                        dist = np.sqrt((xc-X)**2+(yc-Y)**2)
                        if (dist < 3):
                            Comps=np.append(Comps,True)
                            Ncomp=Ncomp + 1
                        else:
                            Comps=np.append(Comps,False)

                    if tmp[0] == "9)":    # axis ratio
                        q=float(tmp[1])

                    if tmp[0] == "10)": # position angle
                        pa=float(tmp[1])
            else:
                Comps=np.append(Comps,False)

        index += 1

    GalfitFile.close()


    return Comps,Ncomp

def GetExpTime(Image):
    "Get exposition time from the image"

    hdu = fits.open(Image)
    exptime = hdu[0].header.get("EXPTIME",1) # return 1 if not found

    hdu.close()
    return exptime


def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"

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


