
import numpy, scipy, math, sys
import keparray
from math import modf, cos, sin, radians, exp
from scipy import ndimage, interpolate
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate
from scipy.interpolate import RectBivariateSpline, interp2d
from keparray import rebin2D
from numpy import square, nansum, shape, array, empty, zeros, absolute, size
from sys import stdout, exit


def PRF(params,*args):

    # arguments

    DATx = args[0]
    DATy = args[1]
    DATimg = args[2]
    DATerr = args[3]
    nsrc = args[4]
    splineInterpolation = args[5]
    col = args[6]
    row = args[7]

    # parameters

    f = empty((nsrc))
    x = empty((nsrc))
    y = empty((nsrc))
    for i in range(nsrc):
        f[i] = params[i]
        x[i] = params[nsrc+i]
        y[i] = params[nsrc*2+i]

    # calculate PRF model binned to the detector pixel size

    PRFfit = PRF2DET(f,x,y,DATx,DATy,1.0,1.0,0.0,splineInterpolation)

    # calculate the sum squared difference between data and model

    #    PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))
    PRFres = nansum(square(DATimg - PRFfit))

    # keep the fit centered

    if max(abs(col - x[0]),abs(row - y[0])) > 10.0:
        PRFres = 1.0e300

    return PRFres


def PRF2DET(flux,OBJx,OBJy,DATx,DATy,wx,wy,a,splineInterpolation):

# trigonometry

    cosa = cos(radians(a))
    sina = sin(radians(a))

# where in the pixel is the source position?

# LUGER: The following lines look horribly, horribly bugged.
# FRCx and FRCy get rewritten each iteration of the for loop,
# so only the last source makes it into PRFfit.

    PRFfit = zeros((size(DATy),size(DATx)))
    for i in range(len(flux)):
        FRCx,INTx = modf(OBJx[i])
        FRCy,INTy = modf(OBJy[i])
        if FRCx > 0.5:
            FRCx -= 1.0
            INTx += 1.0
        if FRCy > 0.5:
            FRCy -= 1.0
            INTy += 1.0
        FRCx = -FRCx
        FRCy = -FRCy

# constuct model PRF in detector coordinates

        for (j,y) in enumerate(DATy):
            for (k,x) in enumerate(DATx):
                xx = x - INTx + FRCx
                yy = y - INTy + FRCy
                dx = xx * cosa - yy * sina
                dy = xx * sina + yy * cosa
                PRFfit[j,k] += PRFfit[j,k] + splineInterpolation(dy*wy,dx*wx) * flux[i]

    return PRFfit
