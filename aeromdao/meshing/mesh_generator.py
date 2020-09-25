#!/usr/bin/env python3
"""
This script reads a coarse airfoil profile, refines the profile using a spline function, and outputs it as surface mesh. Then it generates a 3D volume mesh with nSpan layers in the z direction using pyHyp, available at https://github.com/mdolab/pyhyp.
Note: the airfoil data should be seperated into PS and SS surfaces, they should start from the LE and ends at TE. We use blunt TE so truncate the PS and SS data at about 99.8% of the chord.
"""

from pyhyp import pyHyp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from pyspline import *
from scipy import interpolate

def mesh_foil_2D(mesh_path, lower_surf, upper_surf, layers, wall_dist, march_dist=30.0, te_points=7):
    ZSpan=0.1                              # width in the z direction
    nSpan=2                                # how many points in z

    # TE parameters
    NpTE = te_points                               # number of points for blunt TE

    # 3D
    NpExtrude=layers                           # how many points to extrude for the 3D volume mesh in y
    yWall= wall_dist                            # first layer mesh length
    marchDist= march_dist                         # march distance for extruding

    ########## user input ################

    # read profiles

    # PS parameters
    Alpha1PS = 1.2                                                              # clustering from the LE
    Alpha2PS = 1.2                                                              # clustering from the TE

    # fPS=open(airfoilProfilePS,'r')
    # linesPS=fPS.readlines()
    # fPS.close()

    # xs = [ list(map(float,line.split()))[0] for line in linesPS ]
    xPS = lower_surf[:,0]               
    yPS = lower_surf[:,1]
    
    dX1PS = xPS[1] if xPS[0]==0. else xPS[0]                                       # first dx from the LE
    dX2PS = xPS[-1] - xPS[-2]                                                        # first dx from the TE
    dXMaxPS = round(max([ abs(j - i) for i, j in zip(xPS[:-1], xPS[1:])]), 5)     # max dx for PS

    zPS=[]

    for i in range(len(xPS)):
        zPS.append(0.0)

    # SS parameters             
    Alpha1SS = 1.2  # clustering from the LE                                 
    Alpha2SS = 1.2  # clustering from the TE  



    # xs = [ list(map(float,line.split()))[0] for line in linesSS ]
    xSS = upper_surf[:,0]
    ySS = upper_surf[:,1]

    dX1SS = xSS[1] if xSS[0]==0. else xSS[0]  # first dx from the LE
    dX2SS = xSS[-1] - xSS[-2]                # first dx from the TE
    dXMaxSS = round(max([ abs(j - i) for i, j in zip(xSS[:-1], xSS[1:])]), 5)     # max dx for SS

    zSS=[]

    for i in range(len(xSS)):
        zSS.append(0.0)

    #------------PS
    # first compute how many stretching points do we need for both ends
    tmp=dX1PS
    for i in range(1000):
        if tmp>=dXMaxPS:
            nStretch1PS=i
            break
        else:
            tmp = tmp*Alpha1PS
    #print (nStretch1PS)
    tmp=dX2PS
    for i in range(1000):
        if tmp>=dXMaxPS:
            nStretch2PS=i
            break
        else:
            tmp = tmp*Alpha2PS


    # now compute how much length does these two stretching end and the constant portion have
    xLPSConst=xPS[-1]
    xLPS1=0
    xLPS2=0
    for i in range(nStretch1PS):
        xLPS1 += dX1PS*(Alpha1PS**i)
        xLPSConst -= dX1PS*(Alpha1PS**i)
    for i in range(nStretch2PS):
        xLPS2 += dX2PS*(Alpha2PS**i)
        xLPSConst -= dX2PS*(Alpha2PS**i)

    # update dXMax, we just want to make sure these three portion add up
    nXConstPS = int(xLPSConst/dXMaxPS)
    dXMaxPS=xLPSConst/nXConstPS

    # now we can add these three portions together
    xInterpPS=[0]
    tmp=dX1PS
    for i in range(nStretch1PS):
        xInterpPS.append(xInterpPS[-1]+tmp)
        tmp = tmp*Alpha1PS
    for i in range(nXConstPS):
        xInterpPS.append(xInterpPS[-1]+dXMaxPS)
    tmp=dX2PS*(Alpha2PS**(nStretch2PS-1))
    for i in range(nStretch2PS):
        xInterpPS.append(xInterpPS[-1]+tmp)
        tmp /= Alpha2PS
    # Finally, we interpolate the refined stretch stuff
    
    c1PS = pySpline.Curve(x=xPS, y=yPS, z=zPS, k=3)
    XPS = c1PS(xInterpPS)
    c2PS = pySpline.Curve(X=XPS, k=3)
    x1PS = c1PS.X[:,0]
    # print(x1PS)
    y1PS = c1PS.X[:,1]
    # print(y1PS)
    #------------SS
    # first compute how many stretching points do we need for both ends
    tmp=dX1SS
    for i in range(1000):
        if tmp>=dXMaxSS:
            nStretch1SS=i
            break
        else:
            tmp = tmp*Alpha1SS
    #print (nStretch1SS)
    tmp=dX2SS
    for i in range(1000):
        if tmp>=dXMaxSS:
            nStretch2SS=i
            break
        else:
            tmp = tmp*Alpha2SS
    #print (nStretch2SS)

    # now compute how much length does these two stretching end and the constant portion have
    xLSSConst=xSS[-1]
    xLSS1=0
    xLSS2=0
    for i in range(nStretch1SS):
        xLSS1 += dX1SS*(Alpha1SS**i)
        xLSSConst -= dX1SS*(Alpha1SS**i)
    for i in range(nStretch2SS):
        xLSS2 += dX2SS*(Alpha2SS**i)
        xLSSConst -= dX2SS*(Alpha2SS**i)
    #print xLSS1,xLSS2,xLSSConst
    # update dXMax, we just want to make sure these three portion add up
    nXConstSS = int(xLSSConst/dXMaxSS)
    dXMaxSS=xLSSConst/nXConstSS
    #print(dXMaxSS)

    # now we can add these three portions together
    xInterpSS=[0]
    tmp=dX1SS
    for i in range(nStretch1SS):
        xInterpSS.append(xInterpSS[-1]+tmp)
        tmp = tmp*Alpha1SS
    for i in range(nXConstSS):
        xInterpSS.append(xInterpSS[-1]+dXMaxSS)
    tmp=dX2SS*(Alpha2SS**(nStretch2SS-1))
    for i in range(nStretch2SS):
        xInterpSS.append(xInterpSS[-1]+tmp)
        tmp /= Alpha2SS
    # print xInterpSS
    # Finally, we interpolate the refined stretch stuff
    # print([ (x, y) for (x,y) in zip(xSS, ySS)])
    c1SS = pySpline.Curve(x=xSS, y=ySS, z=zSS, k=3)
    XSS = c1SS(xInterpSS)
    # print([ (x, y, z) for (x,y,z) in XSS ])
    c2SS = pySpline.Curve(X=XSS, k=3)
    x1SS=c1SS.X[:,0]
    y1SS=c1SS.X[:,1]
    # print(y1PS[-1], y1SS[-1])
    # Since the TE is open we need to close it. Close it multiple linear segments.
    delta_y = np.linspace(y1PS[-1], y1SS[-1],NpTE,'d')
    # print(delta_y)
    delta_y = delta_y[1:]
    delta_x=np.ones_like(delta_y,'d')
    for i in range(len(delta_x)):
        delta_x[i]=x1SS[-1]

    x1SS_Flip=x1SS[::-1] # reverse the array #np.flip(x1SS,axis=0)
    xAll=np.append(x1SS_Flip,x1PS[1:])
    xAll=np.append(xAll,delta_x)

    y1SS_Flip=y1SS[::-1] # reverse the array # np.flip(y1SS,axis=0)
    yAll=np.append(y1SS_Flip,y1PS[1:])
    yAll=np.append(yAll,delta_y)

    # print mesh statistics
    print('nPoints for PS: ',nStretch1PS+nStretch2PS+nXConstPS)
    print('nPoints for SS: ',nStretch1SS+nStretch2SS+nXConstSS)
    print('nPoints for TE: ',NpTE)
    print('nPoints Total: ',nStretch1PS+nStretch2PS+nXConstPS+nStretch1SS+nStretch2SS+nXConstSS+NpTE)
    print('Mesh cells: ',(nStretch1PS+nStretch2PS+nXConstPS+nStretch1SS+nStretch2SS+nXConstSS+NpTE-1)*(NpExtrude-1)*(nSpan-1))

    # Write the plot3d input file:
    f = open(f'{mesh_path}.xyz', 'w')
    f.write('1\n')
    f.write('%d %d %d\n'%(len(xAll), nSpan, 1))
    for iDim in range(3):
        for z in np.linspace(0.0,ZSpan,nSpan):
            for i in range(len(xAll)):
                if iDim == 0:
                    f.write('%20.16f\n'%xAll[i])
                elif iDim == 1:
                    f.write('%20.16f\n'%yAll[i])
                else:
                    f.write('%20.16f\n'%z)
    f.close()

    options= {
        # ---------------------------
        #        Input Parameters
        # ---------------------------
        'inputFile':f'{mesh_path}.xyz',
        'unattachedEdgesAreSymmetry':False,
        'outerFaceBC':'farField',
        'autoConnect':True,
        'BC':{1:{'jLow':'zSymm',
                'jHigh':'zSymm'}},
        'families':'wall',

        # ---------------------------
        #        Grid Parameters
        # ---------------------------
        'N': NpExtrude,
        's0':yWall,
        'marchDist':marchDist,
        
        # ---------------------------
        #   Pseudo Grid Parameters
        # ---------------------------
        'ps0':-1,
        'pGridRatio':-1,
        'cMax':1.0,

        # ---------------------------
        #   Smoothing parameters
        # ---------------------------
        'epsE': 2.0,
        'epsI': 4.0,
        'theta': 2.0,
        'volCoef': .20,
        'volBlend': 0.0005,
        'volSmoothIter': 20,
    }


    hyp = pyHyp(options=options)
    hyp.run()
    hyp.writeCGNS(f'{mesh_path}-{format(NpExtrude)}.cgns')
    hyp.writePlot3D(f'{mesh_path}-{format(NpExtrude)}.xyz')

    return mesh_path

# fig1 = plt.figure(1, figsize=(20,20), dpi=300)
# plt.plot(xPS,yPS,'-k',linewidth=1)
# plt.plot(xAll,yAll,'ro',markersize=2)
# plt.plot(xSS,ySS,'-k',linewidth=1)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.xlim([-0.02,1.02])
# plt.show()
# fig1.savefig('airfoil_plot.png', bbox_inches='tight')   # save the figure to file
