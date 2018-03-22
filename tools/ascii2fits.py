#!/usr/bin/python3.5

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




def main():

    if len(sys.argv[1:]) != 2:
        print ('Missing arguments')
        print ("Usage:\n %s [InFile] [FitsFile] " % (sys.argv[0]))

        sys.exit()


    InFile   = sys.argv[1]
    FileFits = sys.argv[2]



    Ascii2Table(InFile, FileFits)



def Ascii2Table(infile, filefits):
    "Converts an ascii table file to another fits table file"

    if os.path.exists(filefits):
        print("{}: filename exist! it will be overwritten ".format(filefits))
        os.remove(filefits)

    ban = 1


#>>> from astropy.table import Table
#>>> t = Table([[1, 2], [4, 5], [7, 8]], names=('a', 'b', 'c'))
#>>> t.write('table.fits', format='fits')

#>>> from astropy.io import fits
#>>> import numpy as np
#>>> c1 = fits.Column(name='a', array=np.array([1, 2]), format='K')
#>>> c2 = fits.Column(name='b', array=np.array([4, 5]), format='K')
#>>> c3 = fits.Column(name='c', array=np.array([7, 8]), format='K')
#>>> t = fits.BinTableHDU.from_columns([c1, c2, c3])
#>>> t.writeto('table.fits')

    (name,ext,dgtrash)=infile.split(".")

    if ext == "ser":
        fitfunc = "Sersic"

    elif ext == "bd":
        fitfunc = "BD"
    else:
        print("Ascii2Table: file extention not recognized\n")
        sys.exit()


    if (fitfunc == "BD"):
# 79 cols

        (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, PosXExp, ErPosXExp, PosYExp, ErPosYExp, MagExp, ErMagExp, RsExp, ErRsExp, AxisExp, ErAxisExp, PaExp, ErPaExp, MeanMeExp, ErMeanMeExp, MeExp, ErMeExp, MsExp, ErMsExp, Sky, ErSky, MagTotal, ErMagTotal, BulgeTotal, ErBt, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit) = np.genfromtxt(infile, delimiter="", unpack=True)


        col1 = fits.Column(name='SExtractor Num Object', format='J', array=Num)
        col2 = fits.Column(name='SExtractor RA', format='E', array=RA)
        col3 = fits.Column(name='SExtractor DEC', format='E', array=Dec)
        col4 = fits.Column(name='SExtractor Posx', format='E', array=XPos)
        col5 = fits.Column(name='SExtractor Posy', format='E', array=YPos)
        col6 = fits.Column(name='SExtractor Mag', format='E', array=Mag)
        col7 = fits.Column(name='SExtractor KronRadius', format='E', array=Kron)
        col8 = fits.Column(name='SExtractor FluxRadius', format='E', array=FluxRad)
        col9 = fits.Column(name='SExtractor IsoArea', format='E', array=IsoArea)
        col10 = fits.Column(name='SExtractor AImage', format='E', array=AIm)

        col11 = fits.Column(name='SExtractor AxisRatio', format='E', array=AR)
        col12 = fits.Column(name='SExtractor PositionAngle', format='E', array=Angle)
        col13 = fits.Column(name='SExtractor Background', format='E', array=Background)
        col14 = fits.Column(name='SExtractor ClassStar', format='E', array=Class)
        col15 = fits.Column(name='SExtractor Flag', format='J', array=Flag)
        col16 = fits.Column(name='Xmin', format='E', array=XMin)
        col17 = fits.Column(name='Xmax', format='E', array=XMax)
        col18 = fits.Column(name='Ymin', format='E', array=YMin)
        col19 = fits.Column(name='Ymax', format='E', array=YMax)
        col20 = fits.Column(name='SersicInput', format='E', array=Sersic)

        col21 = fits.Column(name='SkyKronRadius', format='E', array=RSky)
        col22 = fits.Column(name='KronRadius', format='E', array=RKron)
        col23 = fits.Column(name='SkyFlag', format='J', array=SkyFlag)
        col24 = fits.Column(name='OffX', format='J', array=OFFX)
        col25 = fits.Column(name='OffY', format='J', array=OFFY)
#        col26 = fits.Column(name='OutputFlag', format='J', array=lflag)
        col26 = fits.Column(name='Input SubImage', format='20A', array=InputImage)
        col27 = fits.Column(name='Output Image', format='40A', array=OutIm)
        col28 = fits.Column(name='Object File', format='10A', array=Objx)
        col29 = fits.Column(name='Num Object', format='J', array=PPPnum)
        col30 = fits.Column(name='Restart File', format='10A', array=Restart)


        col31 = fits.Column(name='Input Image', format='10A', array=Cluster)
        col32 = fits.Column(name='RA', format='30A', array=ra)
        col33 = fits.Column(name='DEC', format='30A', array=dec)
        col34 = fits.Column(name='PosxSer', format='E', array=PosXSer)
        col35 = fits.Column(name='PosxSerErr', format='E', array=ErPosXSer)
        col36 = fits.Column(name='PosySer', format='E', array=PosYSer)
        col37 = fits.Column(name='PosySerErr', format='E', array=ErPosYSer)
        col38 = fits.Column(name='MagSer', format='E', array=MagSer)
        col39 = fits.Column(name='MagSerErr', format='E', array=ErMagSer)
        col40 = fits.Column(name='ReSer', format='E', array=ReSer)
        col41 = fits.Column(name='ReSerErr', format='E', array=ErReSer)
        col42 = fits.Column(name='IndexSer', format='E', array=NSer)

        col43 = fits.Column(name='IndexSerErr', format='E', array=ErNSer)
        col44 = fits.Column(name='AxisRatioSer', format='E', array=AxisSer)
        col45 = fits.Column(name='AxisRatioSerErr', format='E', array=ErAxisSer)
        col46 = fits.Column(name='PositionAngleSer', format='E', array=PaSer)
        col47 = fits.Column(name='PositionAngleSerErr', format='E', array=ErPaSer)
        col48 = fits.Column(name='KSer', format='E', array=KSer)
        col49 = fits.Column(name='MeanMeSer', format='E', array=MeanMeSer)
        col50 = fits.Column(name='ErMeanMeSer', format='E', array=ErMeanMeSer)

        col51 = fits.Column(name='MeSer', format='E', array=MeSer)
        col52 = fits.Column(name='ErMeSer', format='E', array=ErMeSer)
        col53 = fits.Column(name='PosXExp', format='E', array=PosXExp)
        col54 = fits.Column(name='PosxExpErr', format='E', array=ErPosXExp)
        col55 = fits.Column(name='PosyExp', format='E', array=PosYExp)
        col56 = fits.Column(name='PosyExpErr', format='E', array=ErPosYExp)
        col57 = fits.Column(name='MagExp', format='E', array=MagExp)
        col58 = fits.Column(name='MagExpErr', format='E', array=ErMagExp)
        col59 = fits.Column(name='RsExp', format='E', array=RsExp)
        col60 = fits.Column(name='RsExpErr', format='E', array=ErRsExp)


        col61 = fits.Column(name='AxisRatioExp', format='E', array=AxisExp)
        col62 = fits.Column(name='AxisRatioExpErr', format='E', array=ErAxisExp)
        col63 = fits.Column(name='PositionAngleExp', format='E', array=PaExp)
        col64 = fits.Column(name='PositionAngleExpErr', format='E', array=ErPaExp)
        col65 = fits.Column(name='MeanMeExp', format='E', array=MeanMeExp)
        col66 = fits.Column(name='ErMeanMeExp', format='E', array=ErMeanMeExp)
        col67 = fits.Column(name='MeExp', format='E', array=MeExp)
        col68 = fits.Column(name='ErMeExp', format='E', array=ErMeExp)
        col69 = fits.Column(name='MsExp', format='E', array=MsExp)
        col70 = fits.Column(name='ErMsExp', format='E', array=ErMsExp)

        col71 = fits.Column(name='Sky', format='E', array=Sky)
        col72 = fits.Column(name='SkyErr', format='E', array=ErSky)
        col73 = fits.Column(name='MagTotal', format='E', array=MagTotal)
        col74 = fits.Column(name='MagTotalErr', format='E', array=ErMagTotal)
        col75 = fits.Column(name='BulgeTotal', format='E', array=BulgeTotal)
        col76 = fits.Column(name='BulgeTotalErr', format='E', array=ErBt)
        col77 = fits.Column(name='ChiNu', format='E', array=ChiNu)
        col78 = fits.Column(name='Tidal', format='E', array=Tidal)
        col79 = fits.Column(name='LocalChiNu', format='E', array=ObjChiNu)
        col80 = fits.Column(name='Bumpiness', format='E', array=Bump)

        col81 = fits.Column(name='Signal to Noise Ratio', format='E', array=SNR)
        col82 = fits.Column(name='Num Degrees Of Freedom', format='E', array=NDof)
        col83 = fits.Column(name='FitFlag', format='E', array=FitFlag)
        col84 = fits.Column(name='ReFit Flag', format='E', array=ReFit)


        t=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,
                         col41, col42, col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, col62, col63, col64, col65, col66, col67, col68, col69, col70, col71, col72, col73, col74, col75, col76, col77, col78, col79, col80, col81, col82, col83, col84])



    else:
#64 cols

        (Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, AR, Angle, Background, Class, Flag, XMin, XMax, YMin, YMax, Sersic, RSky, RKron, SkyFlag, OFFX, OFFY, InputImage, OutIm, Objx, PPPnum, Restart, Cluster, ra, dec, PosXSer, ErPosXSer, PosYSer, ErPosYSer, MagSer, ErMagSer, ReSer, ErReSer, NSer, ErNSer, AxisSer, ErAxisSer, PaSer, ErPaSer, KSer, MeanMeSer, ErMeanMeSer, MeSer, ErMeSer, Sky, ErSky, ChiNu, Tidal, ObjChiNu, Bump, SNR, NDof, FitFlag, ReFit)= np.genfromtxt(infile, delimiter="", unpack=True)


        col1 = fits.Column(name='SExtractor Num Object', format='J', array=Num)
        col2 = fits.Column(name='SExtractor RA', format='E', array=RA)
        col3 = fits.Column(name='SExtractor DEC', format='E', array=Dec)
        col4 = fits.Column(name='SExtractor Posx', format='E', array=XPos)
        col5 = fits.Column(name='SExtractor Posy', format='E', array=YPos)
        col6 = fits.Column(name='SExtractor Mag', format='E', array=Mag)
        col7 = fits.Column(name='SExtractor KronRadius', format='E', array=Kron)
        col8 = fits.Column(name='SExtractor FluxRadius', format='E', array=FluxRad)
        col9 = fits.Column(name='SExtractor IsoArea', format='E', array=IsoArea)
        col10 = fits.Column(name='SExtractor AImage', format='E', array=AIm)

        col11 = fits.Column(name='SExtractor AxisRatio', format='E', array=AR)
        col12 = fits.Column(name='SExtractor PositionAngle', format='E', array=Angle)
        col13 = fits.Column(name='SExtractor Background', format='E', array=Background)
        col14 = fits.Column(name='SExtractor ClassStar', format='E', array=Class)
        col15 = fits.Column(name='SExtractor Flag', format='J', array=Flag)
        col16 = fits.Column(name='Xmin', format='E', array=XMin)
        col17 = fits.Column(name='Xmax', format='E', array=XMax)
        col18 = fits.Column(name='Ymin', format='E', array=YMin)
        col19 = fits.Column(name='Ymax', format='E', array=YMax)
        col20 = fits.Column(name='SersicInput', format='E', array=Sersic)


        col21 = fits.Column(name='SkyKronRadius', format='E', array=RSky)
        col22 = fits.Column(name='KronRadius', format='E', array=RKron)
        col23 = fits.Column(name='SkyFlag', format='J', array=SkyFlag)
        col24 = fits.Column(name='OffX', format='J', array=OFFX)
        col25 = fits.Column(name='OffY', format='J', array=OFFY)
#        col26 = fits.Column(name='OutputFlag', format='J', array=lflag)
        col26 = fits.Column(name='Input SubImage', format='20A', array=InputImage)
        col27 = fits.Column(name='Output Image', format='40A', array=OutIm)
        col28 = fits.Column(name='Object File', format='10A', array=Objx)
        col29 = fits.Column(name='Num Object', format='J', array=PPPnum)
        col30 = fits.Column(name='Restart File', format='10A', array=Restart)


        col31 = fits.Column(name='Input Image', format='10A', array=Cluster)
        col32 = fits.Column(name='RA', format='30A', array=ra)
        col33 = fits.Column(name='DEC', format='30A', array=dec)
        col34 = fits.Column(name='PosxSer', format='E', array=PosXSer)
        col35 = fits.Column(name='PosxSerErr', format='E', array=ErPosXSer)
        col36 = fits.Column(name='PosySer', format='E', array=PosYSer)
        col37 = fits.Column(name='PosySerErr', format='E', array=ErPosYSer)
        col38 = fits.Column(name='MagSer', format='E', array=MagSer)
        col39 = fits.Column(name='MagSerErr', format='E', array=ErMagSer)
        col40 = fits.Column(name='ReSer', format='E', array=ReSer)

        col41 = fits.Column(name='ReSerErr', format='E', array=ErReSer)
        col42 = fits.Column(name='IndexSer', format='E', array=NSer)
        col43 = fits.Column(name='IndexSerErr', format='E', array=ErNSer)
        col44 = fits.Column(name='AxisRatioSer', format='E', array=AxisSer)
        col45 = fits.Column(name='AxisRatioSerErr', format='E', array=ErAxisSer)
        col46 = fits.Column(name='PositionAngleSer', format='E', array=PaSer)
        col47 = fits.Column(name='PositionAngleSerErr', format='E', array=ErPaSer)
        col48 = fits.Column(name='KSer', format='E', array=KSer)
        col49 = fits.Column(name='MeanMeSer', format='E', array=MeanMeSer)
        col50 = fits.Column(name='ErMeanMeSer', format='E', array=ErMeanMeSer)


        col51 = fits.Column(name='MeSer', format='E', array=MeSer)
        col52 = fits.Column(name='ErMeSer', format='E', array=ErMeSer)
        col53 = fits.Column(name='Sky', format='E', array=Sky)
        col54 = fits.Column(name='SkyErr', format='E', array=ErSky)
        col55 = fits.Column(name='ChiNu', format='E', array=ChiNu)
        col56 = fits.Column(name='Tidal', format='E', array=Tidal)
        col57 = fits.Column(name='LocalChiNu', format='E', array=ObjChiNu)
        col58 = fits.Column(name='Bumpiness', format='E', array=Bump)
        col59 = fits.Column(name='Signal to Noise Ratio', format='E', array=SNR)
        col60 = fits.Column(name='Num Degrees Of Freedom', format='E', array=NDof)

        col61 = fits.Column(name='FitFlag', format='E', array=FitFlag)
        col62 = fits.Column(name='ReFit Flag', format='E', array=ReFit)


        t=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36, col37, col38, col39, col40,
                         col41, col42, col43, col44, col45, col46, col47, col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59, col60, col61, col62])




# add comments to header

#tbhdu.header.add_comment('Flags value 0: good fit')
#tbhdu.header.add_comment('Flags value 1: ran into constraints during the fit')
#tbhdu.header.add_comment('Flags value 2: parameters values makes non sense ')
#tbhdu.header.add_comment('Flags value 4: need extra component (bar) ')
#tbhdu.header.add_comment('Flags value 8: bad fit')
#tbhdu.header.add_comment('Flags value 16: cD galaxy ')
#tbhdu.header.add_comment('Flags value 32: fit models are not compatible with galaxy. Example: galaxy distorsioned or need model to fit spiral arms.')
#tbhdu.header.add_comment('Flags value 64: Delete it. Probably a star, galaxy too faint or has saturated pixels')

    t.writeto(filefits)  # writing output

    return True







#end of program
if __name__ == '__main__':
    main()
