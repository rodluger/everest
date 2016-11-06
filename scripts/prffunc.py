
import numpy, scipy, math, sys
from math import modf, cos, sin, radians, exp
from scipy import ndimage, interpolate
from scipy.ndimage import interpolation
from scipy.ndimage.interpolation import shift, rotate
from scipy.interpolate import RectBivariateSpline, interp2d
from numpy import square, nansum, shape, array, empty, zeros, absolute, size
from sys import stdout, exit
import pyfits


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

    PRFres = nansum(square(DATimg - PRFfit) / square(DATerr))
    # PRFres = nansum(square(DATimg - PRFfit))

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

def openfits(file,mode,logfile,verbose):

    status = 0
    try:
        struct = pyfits.open(file,mode=mode)
    except:
        message = 'ERROR -- prffunc.openfits: cannot open ' + file + ' as a FITS file'
        struct = None
        status = 1
    return struct, status

def closefits(struct,logfile,verbose):

    status = 0
    try:
        struct.close()
    except:
        message = 'ERROR -- prffunc.closefits: cannot close HDU structure'
        status = 1
    return status

def get(file,hdu,keyword,logfile,verbose):

    status = 0
    try:
        value = hdu.header[keyword]
    except:
        message = 'ERROR -- prffunc.get: Cannot read keyword ' + keyword
        message += ' in file ' + file
        status = 1
        value = None
    return value, status

# get WCS keywords
def getWCSp(file,struct,logfile,verbose):

    status = 0
    crpix1p = 0.0
    crpix2p = 0.0
    crval1p = 0.0
    crval2p = 0.0
    cdelt1p = 0.0
    cdelt2p = 0.0
    try:
        crpix1p, status = get(file,struct,'CRPIX1P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CRPIX1P in file ' + file
        status = 1
    try:
        crpix2p, status = get(file,struct,'CRPIX2P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CRPIX2P in file ' + file
        status = 1
    try:
        crval1p, status = get(file,struct,'CRVAL1P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CRVAL1P in file ' + file
        status = 1
    try:
        crval2p, status = get(file,struct,'CRVAL2P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CRVAL2P in file ' + file
        status = 1
    try:
        cdelt1p, status = get(file,struct,'CDELT1P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CDELT1P in file ' + file
        status = 1
    try:
        cdelt2p, status = get(file,struct,'CDELT2P',logfile,verbose)
    except:
        txt = 'WARNING -- prffunc.getWCSp: Cannot read keyword CDELT2P in file ' + file
        status = 1

    return crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status

def readPRFimage(infile,hdu,logfile,verbose):

    status = 0

# open input file

    prf, status = openfits(infile,'readonly',logfile,verbose)

# read bitmap image

    if status == 0:
        try:
            img = prf[hdu].data
        except:
            txt = 'ERROR -- prffunc.readPRFimage: Cannot read PRF image in ' + infile + '[' + str(hdu) + ']'
            # status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            naxis1 = prf[hdu].header['NAXIS1']
        except:
            txt = 'ERROR -- prffunc.readPRFimage: Cannot read NAXIS1 keyword in ' + infile + '[' + str(hdu) + ']'
            # status = kepmsg.err(logfile,txt,verbose)
    if status == 0:
        try:
            naxis2 = prf[hdu].header['NAXIS2']
        except:
            txt = 'ERROR -- prffunc.readPRFimage: Cannot read NAXIS2 keyword in ' + infile + '[' + str(hdu) + ']'
            # status = kepmsg.err(logfile,txt,verbose)

# read WCS keywords

    if status == 0:
        crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status = \
            getWCSp(infile,prf[hdu],logfile,verbose)

# close input file

    if status == 0:
        status = closefits(prf,logfile,verbose)

    return img, crpix1p, crpix2p, crval1p, crval2p, cdelt1p, cdelt2p, status
