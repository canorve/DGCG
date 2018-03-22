#!/usr/bin/python


import numpy as np
import sys
import os



RunDir="Directory"

ReRunDir="ReDirectory"


ReFit=[True,False,False,False,False,True,True,False,False,True,False,True]
Num=[34,44,78,999,5434,3,7,98,12,55,36,12]

ReFit=np.array(ReFit)
Num=np.array(Num)

maskrun= ReFit == True

ind=np.where(maskrun == True)

indx=ind[0]


print ("ind",ind)

print ("indx",indx)


for idx in enumerate(indx):

    print ("idx",idx)


    objid = Num[idx[1]]

    parfile =  "obj" + "-" + str(objid)

    pmsg="copying file {} for refit ".format(parfile)
    print(pmsg)

    dirparfile =  RunDir + "/" + parfile

    redirparfile =  ReRunDir + "/" + parfile

    print(dirparfile)
    print(redirparfile)




#        runcmd = "cp {} {}".format(dirparfile,redirparfile)
#        errcp = sp.run([runcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

print("copying files to refit done..")
