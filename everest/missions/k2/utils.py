#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`utils.py` - Mission auxiliary routines
-----------------------------------------------

`K2`-specific auxiliary routines. These are not
generally called from the top level of the code.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from .pipelines import Pipelines
from ...config import EVEREST_SRC, EVEREST_DAT, EVEREST_DEV
from ...utils import _float
from ...mathutils import Chunks
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        raise Exception('Please install the `pyfits` package.')
from astropy.wcs import WCS
from scipy.interpolate import griddata
from k2plr.api import K2_CAMPAIGNS
import numpy as np
from tempfile import NamedTemporaryFile
from six.moves import urllib
import re
import os
import subprocess
import logging
import k2plr as kplr
log = logging.getLogger(__name__)
kplr_client = kplr.API()

__all__ = ['Campaign', 'GetK2Stars', 'GetK2Campaign', 'Channel',
           'RemoveBackground', 'GetNeighboringChannels', 'GetSources',
           'GetHiResImage', 'GetCustomAperture',
           'StatsPicker', 'SaturationFlux', 'Module', 'Channels']


def _range10_90(x):
    '''
    Returns the 10th-90th percentile range of array :py:obj:`x`.

    '''

    x = np.delete(x, np.where(np.isnan(x)))
    i = np.argsort(x)
    a = int(0.1 * len(x))
    b = int(0.9 * len(x))
    return x[i][b] - x[i][a]


class StatsPicker(object):
    '''
    A class that enables clicking on the individual points on
    the :py:func:`k2.Statistics` scatter plots.

    :param axes: A :py:mod:`matplotlib.pyplot` axis instance or a \
           list of axis instances
    :param x: An array or a list of arrays corresponding to the \
           abscissa of :py:obj:`axes`
    :param y: An array or a list of arrays corresponding to the \
           ordinate of :py:obj:`axes`
    :param array_like epic: A list of EPIC target numbers for \
           each of the plotted points
    :param str model: The name of the current :py:mod:`everest` \
           model. Default `"PLD"`
    :param str compare_to: The name of the model against which the data \
           is being compared. Default `"k2sff"`

    '''

    def __init__(self, axes, x, y, epic, model='PLD', compare_to='k2sff',
                 cadence='lc', campaign=None):
        '''

        '''

        from ...user import DVS
        self.show = DVS
        if not hasattr(axes, '__len__'):
            axes = [axes]
            x = [x]
            y = [y]
        self.axes = axes
        self.x = [np.array(xi) for xi in x]
        self.y = [np.array(yi) for yi in y]
        self.xr = [_range10_90(x) for x in self.x]
        self.yr = [_range10_90(y) for y in self.y]
        self.epic = epic
        self.model = model
        self.compare_to = compare_to
        self.last = None
        self.cadence = cadence
        self.campaign = campaign

    def __call__(self, event):
        '''

        '''

        if event.mouseevent.inaxes:

            # Get the axis instance
            j = np.argmax([id(event.mouseevent.inaxes) == id(ax)
                           for ax in self.axes])

            # Index of nearest point
            i = np.nanargmin(((event.mouseevent.xdata - self.x[j]) /
                              self.xr[j]) ** 2 + (
                (event.mouseevent.ydata - self.y[j]) / self.yr[j]) ** 2)

            # HACK: For some reason, this event is being called twice
            # for every click. This is a silly way around that.
            if self.epic[i] == self.last:
                return
            else:
                self.last = self.epic[i]

            # Show the de-trended data for the model
            log.info('Plotting %s model for %d...' %
                     (self.model, self.epic[i]))
            self.show(self.epic[i], mission='k2',
                      cadence=self.cadence, season=self.campaign)

            # Show the de-trended data for the comparison model
            if self.compare_to.lower() in Pipelines:
                log.info('Plotting %s model for %d...' %
                         (self.compare_to, self.epic[i]))
                cmd = ['python', '-c',
                       'import everest; everest.k2.pipelines.' +
                       'plot(%d, pipeline="%s"%s)' %
                       (self.epic[i], self.compare_to,
                        ", campaign=%d" % self.campaign
                        if self.campaign is not None
                        else "")]
                print(" ".join(cmd))
                subprocess.Popen(cmd)
            elif self.compare_to.lower() == 'kepler':
                pass
            else:
                log.info('Plotting %s model for %d...' %
                         (self.compare_to, self.epic[i]))
                self.show(self.epic[i], mission='k2', model=self.compare_to)


def Campaign(EPIC, **kwargs):
    '''
    Returns the campaign number(s) for a given EPIC target. If target
    is not found, returns :py:obj:`None`.

    :param int EPIC: The EPIC number of the target.

    '''

    campaigns = []
    for campaign, stars in GetK2Stars().items():
        if EPIC in [s[0] for s in stars]:
            campaigns.append(campaign)
    if len(campaigns) == 0:
        return None
    elif len(campaigns) == 1:
        return campaigns[0]
    else:
        return campaigns


def GetK2Stars(clobber=False):
    '''
    Download and return a :py:obj:`dict` of all *K2* stars organized by
    campaign. Saves each campaign to a `.stars` file in the
    `everest/missions/k2/tables` directory.

    :param bool clobber: If :py:obj:`True`, download and overwrite \
           existing files. Default :py:obj:`False`

    .. note:: The keys of the dictionary returned by this function are the \
              (integer) numbers of each campaign. Each item in the \
              :py:obj:`dict` is a list of the targets in the corresponding \
              campaign, and each item in that list is in turn a list of the \
              following: **EPIC number** (:py:class:`int`), \
              **Kp magnitude** (:py:class:`float`), **CCD channel number** \
              (:py:class:`int`), and **short cadence available** \
              (:py:class:`bool`).

    '''

    # Download
    if clobber:
        print("Downloading K2 star list...")
        stars = kplr_client.k2_star_info()
        print("Writing star list to disk...")
        for campaign in stars.keys():
            if not os.path.exists(os.path.join(EVEREST_SRC, 'missions',
                                               'k2', 'tables')):
                os.makedirs(os.path.join(
                    EVEREST_SRC, 'missions', 'k2', 'tables'))
            with open(os.path.join(EVEREST_SRC, 'missions', 'k2', 'tables',
                                   'c%02d.stars' % campaign), 'w') as f:
                for star in stars[campaign]:
                    print(",".join([str(s) for s in star]), file=f)

    # Return
    res = {}
    for campaign in K2_CAMPAIGNS:
        f = os.path.join(EVEREST_SRC, 'missions', 'k2',
                         'tables', 'c%02d.stars' % campaign)
        if os.path.exists(f):
            with open(f, 'r') as file:
                lines = file.readlines()
                if len(lines[0].split(',')) == 4:
                    # EPIC number, Kp magnitude, channel number,
                    # short cadence available?
                    stars = [[int(l.split(',')[0]),
                              _float(l.split(',')[1]),
                              int(l.split(',')[2]),
                              eval(l.split(',')[3])] for l in lines]
                else:
                    stars = [[int(l), np.nan, -1, None] for l in lines]
            res.update({campaign: stars})

    return res


def GetK2Campaign(campaign, clobber=False, split=False,
                  epics_only=False, cadence='lc'):
    '''
    Return all stars in a given *K2* campaign.

    :param campaign: The *K2* campaign number. If this is an :py:class:`int`, \
           returns all targets in that campaign. If a :py:class:`float` in \
           the form :py:obj:`X.Y`, runs the :py:obj:`Y^th` decile of campaign \
           :py:obj:`X`.
    :param bool clobber: If :py:obj:`True`, download and overwrite existing \
           files. Default :py:obj:`False`
    :param bool split: If :py:obj:`True` and :py:obj:`campaign` is an \
           :py:class:`int`, returns each of the subcampaigns as a separate \
           list. Default :py:obj:`False`
    :param bool epics_only: If :py:obj:`True`, returns only the EPIC numbers. \
           If :py:obj:`False`, returns metadata associated with each target. \
           Default :py:obj:`False`
    :param str cadence: Long (:py:obj:`lc`) or short (:py:obj:`sc`) cadence? \
           Default :py:obj:`lc`.

    '''

    all = GetK2Stars(clobber=clobber)
    if int(campaign) in all.keys():
        all = all[int(campaign)]
    else:
        return []

    if cadence == 'sc':
        all = [a for a in all if a[3]]

    if epics_only:
        all = [a[0] for a in all]
    if type(campaign) is int or type(campaign) is np.int64:
        if not split:
            return all
        else:
            all_split = list(Chunks(all, len(all) // 10))

            # HACK: Sometimes we're left with a few targets
            # dangling at the end. Insert them back evenly
            # into the first few subcampaigns.
            if len(all_split) > 10:
                tmp1 = all_split[:10]
                tmp2 = all_split[10:]
                for n in range(len(tmp2)):
                    tmp1[n] = np.append(tmp1[n], tmp2[n])
                all_split = tmp1

            res = []
            for subcampaign in range(10):
                res.append(all_split[subcampaign])

            return res
    elif type(campaign) is float:
        x, y = divmod(campaign, 1)
        campaign = int(x)
        subcampaign = round(y * 10)
        return list(Chunks(all, len(all) // 10))[subcampaign]
    else:
        raise Exception('Argument `subcampaign` must be an `int` ' +
                        'or a `float` in the form `X.Y`')


def Channel(EPIC, campaign=None):
    '''
    Returns the channel number for a given EPIC target.

    '''

    if campaign is None:
        campaign = Campaign(EPIC)
    if hasattr(campaign, '__len__'):
        raise AttributeError(
            "Please choose a campaign/season for this target: %s." % campaign)
    try:
        stars = GetK2Stars()[campaign]
    except KeyError:
        # Not sure what else to do here!
        log.warn("Unknown channel for target. Defaulting to channel 2.")
        return 2
    i = np.argmax([s[0] == EPIC for s in stars])
    return stars[i][2]


def Module(EPIC, campaign=None):
    '''
    Returns the module number for a given EPIC target.

    '''

    channel = Channel(EPIC, campaign=campaign)
    nums = {2: 1, 3: 5, 4: 9, 6: 13, 7: 17, 8: 21, 9: 25,
            10: 29, 11: 33, 12: 37, 13: 41, 14: 45, 15: 49,
            16: 53, 17: 57, 18: 61, 19: 65, 20: 69, 22: 73,
            23: 77, 24: 81}
    for c in [channel, channel - 1, channel - 2, channel - 3]:
        if c in nums.values():
            for mod, chan in nums.items():
                if chan == c:
                    return mod
    return None


def Channels(module):
    '''
    Returns the channels contained in the given K2 module.

    '''

    nums = {2: 1, 3: 5, 4: 9, 6: 13, 7: 17, 8: 21, 9: 25,
            10: 29, 11: 33, 12: 37, 13: 41, 14: 45, 15: 49,
            16: 53, 17: 57, 18: 61, 19: 65, 20: 69, 22: 73,
            23: 77, 24: 81}

    if module in nums:
        return [nums[module], nums[module] + 1,
                nums[module] + 2, nums[module] + 3]
    else:
        return None


def KepMag(EPIC, campaign=None):
    '''
    Returns the *Kepler* magnitude for a given EPIC target.

    '''

    if campaign is None:
        campaign = Campaign(EPIC)
    if hasattr(campaign, '__len__'):
        raise AttributeError(
            "Please choose a campaign/season for this target: %s." % campaign)
    stars = GetK2Stars()[campaign]
    i = np.argmax([s[0] == EPIC for s in stars])
    return stars[i][1]


def RemoveBackground(EPIC, campaign=None):
    '''
    Returns :py:obj:`True` or :py:obj:`False`, indicating whether or not
    to remove the background flux for the target. If ``campaign < 3``,
    returns :py:obj:`True`, otherwise returns :py:obj:`False`.

    '''

    if campaign is None:
        campaign = Campaign(EPIC)
    if hasattr(campaign, '__len__'):
        raise AttributeError(
            "Please choose a campaign/season for this target: %s." % campaign)
    if campaign < 3:
        return True
    else:
        return False


def GetNeighboringChannels(channel):
    '''
    Returns all channels on the same module as :py:obj:`channel`.

    '''

    x = divmod(channel - 1, 4)[1]
    return channel + np.array(range(-x, -x + 4), dtype=int)


def MASTRADec(ra, dec, darcsec, stars_only=False):
    '''
    Detector location retrieval based upon RA and Dec.
    Adapted from `PyKE <http://keplergo.arc.nasa.gov/PyKE.shtml>`_.

    '''

    # coordinate limits
    darcsec /= 3600.0
    ra1 = ra - darcsec / np.cos(dec * np.pi / 180)
    ra2 = ra + darcsec / np.cos(dec * np.pi / 180)
    dec1 = dec - darcsec
    dec2 = dec + darcsec

    # build mast query
    url = 'http://archive.stsci.edu/k2/epic/search.php?'
    url += 'action=Search'
    url += '&k2_ra=' + str(ra1) + '..' + str(ra2)
    url += '&k2_dec=' + str(dec1) + '..' + str(dec2)
    url += '&max_records=10000'
    url += '&selectedColumnsCsv=id,k2_ra,k2_dec,kp'
    url += '&outputformat=CSV'
    if stars_only:
        url += '&ktc_target_type=LC'
        url += '&objtype=star'

    # retrieve results from MAST
    try:
        lines = urllib.request.urlopen(url)
    except:
        log.warn('Unable to retrieve source data from MAST.')
        lines = ''

    # collate nearby sources
    epicid = []
    kepmag = []
    ra = []
    dec = []
    for line in lines:

        line = line.strip().decode('ascii')

        if (len(line) > 0 and 'EPIC' not in line and 'integer' not in line and
                              'no rows found' not in line):

            out = line.split(',')
            r, d = sex2dec(out[1], out[2])
            epicid.append(int(out[0]))
            kepmag.append(float(out[3]))
            ra.append(r)
            dec.append(d)

    epicid = np.array(epicid)
    kepmag = np.array(kepmag)
    ra = np.array(ra)
    dec = np.array(dec)

    return epicid, ra, dec, kepmag


def sex2dec(ra, dec):
    '''
    Convert sexadecimal hours to decimal degrees. Adapted from
    `PyKE <http://keplergo.arc.nasa.gov/PyKE.shtml>`_.

    :param float ra: The right ascension
    :param float dec: The declination

    :returns: The same values, but in decimal degrees

    '''

    ra = re.sub('\s+', '|', ra.strip())
    ra = re.sub(':', '|', ra.strip())
    ra = re.sub(';', '|', ra.strip())
    ra = re.sub(',', '|', ra.strip())
    ra = re.sub('-', '|', ra.strip())
    ra = ra.split('|')
    outra = (float(ra[0]) + float(ra[1]) / 60. + float(ra[2]) / 3600.) * 15.0

    dec = re.sub('\s+', '|', dec.strip())
    dec = re.sub(':', '|', dec.strip())
    dec = re.sub(';', '|', dec.strip())
    dec = re.sub(',', '|', dec.strip())
    dec = dec.split('|')

    if float(dec[0]) > 0.0:
        outdec = float(dec[0]) + float(dec[1]) / 60. + float(dec[2]) / 3600.
    else:
        outdec = float(dec[0]) - float(dec[1]) / 60. - float(dec[2]) / 3600.

    return outra, outdec


def GetSources(ID, darcsec=None, stars_only=False):
    '''
    Grabs the EPIC coordinates from the TPF and searches MAST
    for other EPIC targets within the same aperture.

    :param int ID: The 9-digit :py:obj:`EPIC` number of the target
    :param float darcsec: The search radius in arcseconds. \
           Default is four times the largest dimension of the aperture.
    :param bool stars_only: If :py:obj:`True`, only returns objects \
           explicitly designated as `"stars"` in MAST. Default :py:obj:`False`
    :returns: A list of :py:class:`Source` instances containing \
              other :py:obj:`EPIC` targets within or close to this \
              target's aperture
    '''

    client = kplr.API()
    star = client.k2_star(ID)
    tpf = star.get_target_pixel_files()[0]
    with tpf.open() as f:
        crpix1 = f[2].header['CRPIX1']
        crpix2 = f[2].header['CRPIX2']
        crval1 = f[2].header['CRVAL1']
        crval2 = f[2].header['CRVAL2']
        cdelt1 = f[2].header['CDELT1']
        cdelt2 = f[2].header['CDELT2']
        pc1_1 = f[2].header['PC1_1']
        pc1_2 = f[2].header['PC1_2']
        pc2_1 = f[2].header['PC2_1']
        pc2_2 = f[2].header['PC2_2']
        pc = np.array([[pc1_1, pc1_2], [pc2_1, pc2_2]])
        pc = np.linalg.inv(pc)
        crpix1p = f[2].header['CRPIX1P']
        crpix2p = f[2].header['CRPIX2P']
        crval1p = f[2].header['CRVAL1P']
        crval2p = f[2].header['CRVAL2P']
        cdelt1p = f[2].header['CDELT1P']
        cdelt2p = f[2].header['CDELT2P']
        if darcsec is None:
            darcsec = 4 * max(f[2].data.shape)

    epicid, ra, dec, kepmag = MASTRADec(
        star.k2_ra, star.k2_dec, darcsec, stars_only)
    sources = []
    for i, epic in enumerate(epicid):
        dra = (ra[i] - crval1) * np.cos(np.radians(dec[i])) / cdelt1
        ddec = (dec[i] - crval2) / cdelt2
        sx = pc[0, 0] * dra + pc[0, 1] * ddec + crpix1 + crval1p - 1.0
        sy = pc[1, 0] * dra + pc[1, 1] * ddec + crpix2 + crval2p - 1.0
        sources.append(dict(ID=epic, x=sx, y=sy, mag=kepmag[i],
                            x0=crval1p, y0=crval2p))

    return sources


def GetHiResImage(ID):
    '''
    Queries the Palomar Observatory Sky Survey II catalog to
    obtain a higher resolution optical image of the star with EPIC number
    :py:obj:`ID`.

    '''

    # Get the TPF info
    client = kplr.API()
    star = client.k2_star(ID)
    k2ra = star.k2_ra
    k2dec = star.k2_dec
    tpf = star.get_target_pixel_files()[0]
    with tpf.open() as f:
        k2wcs = WCS(f[2].header)
        shape = np.array(f[1].data.field('FLUX'), dtype='float64')[0].shape

    # Get the POSS URL
    hou = int(k2ra * 24 / 360.)
    min = int(60 * (k2ra * 24 / 360. - hou))
    sec = 60 * (60 * (k2ra * 24 / 360. - hou) - min)
    ra = '%02d+%02d+%.2f' % (hou, min, sec)
    sgn = '' if np.sign(k2dec) >= 0 else '-'
    deg = int(np.abs(k2dec))
    min = int(60 * (np.abs(k2dec) - deg))
    sec = 3600 * (np.abs(k2dec) - deg - min / 60)
    dec = '%s%02d+%02d+%.1f' % (sgn, deg, min, sec)
    url = 'https://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&' + \
          'r=%s&d=%s&e=J2000&h=3&w=3&f=fits&c=none&fov=NONE&v3=' % (ra, dec)

    # Query the server
    r = urllib.request.Request(url)
    handler = urllib.request.urlopen(r)
    code = handler.getcode()
    if int(code) != 200:
        # Unavailable
        return None
    data = handler.read()

    # Atomically write to a temp file
    f = NamedTemporaryFile("wb", delete=False)
    f.write(data)
    f.flush()
    os.fsync(f.fileno())
    f.close()

    # Now open the POSS fits file
    with pyfits.open(f.name) as ff:
        img = ff[0].data

    # Map POSS pixels onto K2 pixels
    xy = np.empty((img.shape[0] * img.shape[1], 2))
    z = np.empty(img.shape[0] * img.shape[1])
    pwcs = WCS(f.name)
    k = 0
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            ra, dec = pwcs.all_pix2world(float(j), float(i), 0)
            xy[k] = k2wcs.all_world2pix(ra, dec, 0)
            z[k] = img[i, j]
            k += 1

    # Resample
    grid_x, grid_y = np.mgrid[-0.5:shape[1] - 0.5:0.1, -0.5:shape[0] - 0.5:0.1]
    resampled = griddata(xy, z, (grid_x, grid_y), method='cubic')

    # Rotate to align with K2 image. Not sure why, but it is necessary
    resampled = np.rot90(resampled)

    return resampled


def GetCustomAperture(data):
    '''
    .. warning:: Routine not yet implemented.

    '''

    raise NotImplementedError('TODO: This routine still needs to be written.')


def SaturationFlux(EPIC, campaign=None, **kwargs):
    '''
    Returns the well depth for the target. If any of the target's pixels
    have flux larger than this value, they are likely to be saturated and
    cause charge bleeding. The well depths were obtained from Table 13
    of the Kepler instrument handbook. We assume an exposure time of 6.02s.

    '''

    channel, well_depth = np.loadtxt(os.path.join(EVEREST_SRC, 'missions',
                                                  'k2',
                                                  'tables', 'well_depth.tsv'),
                                     unpack=True)
    satflx = well_depth[channel == Channel(EPIC, campaign=campaign)][0] / 6.02
    return satflx
