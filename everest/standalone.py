#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`standalone.py` - Standalone de-trending
------------------------------------------------

Provides the :py:func:`DetrendFITS` function for
manual de-trending of user-provided `K2` FITS files.


'''

from __future__ import division, print_function, absolute_import
import os
import shutil
import numpy as np
import everest
from everest.mathutils import Interpolate, SavGol
from everest.utils import AP_COLLAPSED_PIXEL, AP_SATURATED_PIXEL, DataContainer
from everest.config import EVEREST_DAT
from everest.missions.k2.utils import GetHiResImage, GetSources, \
     SaturationFlux, RemoveBackground
from tempfile import NamedTemporaryFile
import matplotlib
from matplotlib.widgets import Slider
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as pl
from scipy.ndimage import zoom
try:
    import pyfits
except ImportError:
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        raise Exception('Please install the `pyfits` package.')
import logging
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
log = logging.getLogger(__name__)


def DetrendFITS(fitsfile, raw=False, season=None, clobber=False, **kwargs):
    """
    De-trend a K2 FITS file using :py:class:`everest.detrender.rPLD`.

    :param str fitsfile: The full path to the FITS file
    :param ndarray aperture: A 2D integer array corresponding to the \
           desired photometric aperture (1 = in aperture, 0 = outside \
           aperture). Default is to interactively select an aperture.
    :param kwargs: Any kwargs accepted by :py:class:`everest.detrender.rPLD`.

    :returns: An :py:class:`everest.Everest` instance.

    """
    # Get info
    EPIC = pyfits.getheader(fitsfile, 0)['KEPLERID']
    if season is None:
        season = pyfits.getheader(fitsfile, 0)['CAMPAIGN']
        if season is None or season == "":
            season = 0
    everestfile = os.path.join(
        everest.missions.k2.TargetDirectory(EPIC, season),
        everest.missions.k2.FITSFile(EPIC, season))

    # De-trend?
    if clobber or not os.path.exists(everestfile):

        # Get raw data
        data = GetData(fitsfile, EPIC, season, clobber=clobber, **kwargs)

        # De-trend
        model = everest.rPLD(EPIC,
                             data=data,
                             season=season, debug=True,
                             clobber=clobber, **kwargs)

        # Publish it
        everest.fits.MakeFITS(model)
        shutil.copyfile(os.path.join(model.dir, model.name + '.pdf'),
                        os.path.join(model.dir,
                                     model._mission.DVSFile(model.ID,
                                                            model.season,
                                                            model.cadence)))

    # Return an Everest instance
    return everest.Everest(EPIC, season=season)


class ApertureSelector(object):
    '''

    '''

    def __init__(self, time, images, title='Aperture'):
        '''

        '''

        self.cadence = 0
        self.time = time
        self.fig, self.ax = pl.subplots(1, figsize=(10, 7))
        self.fig.subplots_adjust(left=0.1, bottom=0.25, top=0.925, right=0.45)
        self.images = images
        self.nt, self.ny, self.nx = self.images.shape
        self.x = np.arange(0, self.nx)
        self.y = np.arange(0, self.ny)
        self.aperture = np.zeros((self.ny, self.nx), dtype=int)
        self.aperture[self.ny // 2 - 2:self.ny // 2 +
                      2][:, self.nx // 2 - 2:self.nx // 2 + 2] = 1
        self.contour = None
        self.last_j = None
        self.last_i = None
        self.title = title

        # Slider
        self.axslider = pl.axes([0.105, 0.2, 0.34, 0.03])
        self.slider = Slider(self.axslider, '', 0,
                             self.nt - 1, valinit=0, valfmt='%d')
        self.slider.valtext.set_x(0.5)
        self.slider.valtext.set_ha('center')
        self.slider.on_changed(self.replot)

        # Background
        self.axbkg = pl.axes([0.105, 0.05, 0.34, 0.125])
        bkg = self.colbkg
        self.bkgplot1, = self.axbkg.plot(self.x, bkg, 'ro')
        self.bkgplot2, = self.axbkg.plot(self.x, bkg, 'r-', alpha=0.3)
        pad = 0.2 * (bkg.max() - bkg.min())
        self.axbkg.set_ylim(bkg.min() - pad, bkg.max() + pad)
        self.axbkg.set_xlim(-0.7, self.nx - 0.3)
        for tick in self.axbkg.get_yticklabels():
            tick.set_fontsize(7)
        self.axbkg.get_yaxis().set_major_formatter(
            FuncFormatter(lambda x, p: '%.2f' % x))
        self.axbkg.set_ylabel('Bkg (%)', fontsize=9)

        # Light curve
        self.axlc = pl.axes([0.5, 0.5, 0.4, 0.425])
        self.lcplot, = self.axlc.plot(
            self.time, self.flux, 'k.', alpha=0.3, ms=3)
        self.axlc.set_xticklabels([])
        self.axlc.yaxis.tick_right()
        self.axlc.set_ylabel('Light curve', fontsize=14)
        self.lcstdtxt = self.axlc.annotate('%.2f ppm' % self.lcstd,
                                           xy=(0.025, 0.975),
                                           xycoords='axes fraction',
                                           ha='left', va='top',
                                           fontsize=12, color='r')

        # Light curve background
        self.axlcbkg = pl.axes([0.5, 0.05, 0.4, 0.425])
        self.lcbkgplot, = self.axlcbkg.plot(
            self.time, self.lcbkg, 'k.', alpha=0.3, ms=3)
        self.axlcbkg.yaxis.tick_right()
        self.axlcbkg.set_ylabel('Background', fontsize=14)
        self.bkgstdtxt = self.axlcbkg.annotate('%.2f ppm' % self.bkgstd,
                                               xy=(0.025, 0.975),
                                               xycoords='axes fraction',
                                               ha='left', va='top',
                                               fontsize=12, color='r')

        # Trackers
        self.tracker1 = self.axlc.axvline(
            self.time[self.cadence], color='r', alpha=0.5, lw=1)
        self.tracker2 = self.axlcbkg.axvline(
            self.time[self.cadence], color='r', alpha=0.5, lw=1)

        # Appearance
        self.fig.canvas.set_window_title('Select an aperture')
        self.ax.axis('off')
        self.ax.set_xlim(-0.7, self.nx - 0.3)
        self.ax.set_ylim(-0.7, self.ny - 0.3)
        self.ax.set_title(title, fontsize=18)

        # Plot the image
        try:
            plasma = pl.get_cmap('plasma')
        except ValueError:
            plasma = pl.get_cmap('Greys')
        plasma.set_bad(alpha=0)
        self.implot = self.ax.imshow(self.images[self.cadence],
                                     aspect='auto', interpolation='nearest',
                                     cmap=plasma, picker=True)
        self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_drag)
        self.fig.canvas.mpl_connect('pick_event', self.mouse_click)

        # Update the contour
        self.update()

        # Enter interactive mode
        pl.show()

    @property
    def colbkg(self):
        '''

        '''

        # Flux in background pixels
        bkg = np.zeros(self.nx)
        for col in range(self.nx):
            b = np.where(self.aperture[:, col] == 0)
            bkg[col] = np.nanmedian(self.images[self.cadence][b, col])
        return 100 * (bkg / np.mean(bkg) - 1.)

    @property
    def lcbkg(self):
        '''

        '''

        binds = np.where(self.aperture ^ 1)
        bkg = np.nanmedian(
            np.array([f[binds] for f in self.images], dtype='float64'), axis=1)
        return bkg.reshape(-1, 1)

    @property
    def flux(self):
        '''

        '''

        ap = np.where(self.aperture & 1)
        fpix2D = np.array([f[ap] for f in self.images], dtype='float64')
        return np.sum(fpix2D - self.lcbkg, axis=1)

    @property
    def lcstd(self):
        '''

        '''

        return everest.k2.CDPP(self.flux)

    @property
    def bkgstd(self):
        '''

        '''

        return everest.k2.CDPP(self.lcbkg)

    def update_bkg(self):
        '''

        '''

        bkg = self.colbkg
        self.bkgplot1.set_ydata(bkg)
        self.bkgplot2.set_ydata(bkg)
        pad = 0.2 * (bkg.max() - bkg.min())
        self.axbkg.set_ylim(bkg.min() - pad, bkg.max() + pad)
        self.axbkg.set_xlim(-0.7, self.nx - 0.3)

    def update_lc(self):
        '''

        '''

        flux = self.flux
        self.lcplot.set_ydata(flux)
        pad = 0.2 * (flux.max() - flux.min())
        self.axlc.set_ylim(flux.min() - pad, flux.max() + pad)
        self.axlc.set_xlim(self.time[0], self.time[-1])
        self.lcstdtxt.set_text('%.2f ppm' % self.lcstd)

    def update_lcbkg(self):
        '''

        '''

        lcbkg = self.lcbkg
        self.lcbkgplot.set_ydata(lcbkg)
        pad = 0.2 * (lcbkg.max() - lcbkg.min())
        self.axlcbkg.set_ylim(lcbkg.min() - pad, lcbkg.max() + pad)
        self.axlcbkg.set_xlim(self.time[0], self.time[-1])
        self.bkgstdtxt.set_text('%.2f ppm' % self.bkgstd)

    def PadWithZeros(self, vector, pad_width, iaxis, kwargs):
        '''

        '''

        vector[:pad_width[0]] = 0
        vector[-pad_width[1]:] = 0
        return vector

    def mouse_drag(self, event):
        '''

        '''

        if event.inaxes == self.ax and event.button == 1:

            # Index of nearest point
            i = np.nanargmin(((event.xdata - self.x) / self.nx) ** 2)
            j = np.nanargmin(((event.ydata - self.y) / self.ny) ** 2)

            if (i == self.last_i) and (j == self.last_j):
                return
            else:
                self.last_i = i
                self.last_j = j

            # Toggle pixel
            if self.aperture[j, i]:
                self.aperture[j, i] = 0
            else:
                self.aperture[j, i] = 1

            # Update the contour
            self.update()

    def mouse_click(self, event):
        '''

        '''

        if event.mouseevent.inaxes == self.ax:

            # Index of nearest point
            i = np.nanargmin(
                ((event.mouseevent.xdata - self.x) / self.nx) ** 2)
            j = np.nanargmin(
                ((event.mouseevent.ydata - self.y) / self.ny) ** 2)
            self.last_i = i
            self.last_j = j

            # Toggle pixel
            if self.aperture[j, i]:
                self.aperture[j, i] = 0
            else:
                self.aperture[j, i] = 1

            # Update the contour
            self.update()

    def update(self):
        '''

        '''

        # Update plot
        contour = np.zeros((self.ny, self.nx))
        contour[np.where(self.aperture)] = 1
        contour = np.lib.pad(contour, 1, self.PadWithZeros)
        highres = zoom(contour, 100, order=0, mode='nearest')
        extent = np.array([-1, self.nx, -1, self.ny])
        if self.contour is not None:
            for coll in self.contour.collections:
                self.ax.collections.remove(coll)
        self.contour = self.ax.contour(highres, levels=[0.5], extent=extent,
                                       origin='lower', colors='r',
                                       linewidths=2)
        self.update_bkg()
        self.update_lc()
        self.update_lcbkg()
        self.fig.canvas.draw()

    def replot(self, val):
        '''

        '''

        # Update plot
        self.cadence = int(val)
        self.implot.set_data(self.images[int(val)])
        self.implot.set_clim(vmin=np.nanmin(
            self.images[int(val)]), vmax=np.nanmax(self.images[int(val)]))
        self.tracker1.set_xdata(
            [self.time[self.cadence], self.time[self.cadence]])
        self.tracker2.set_xdata(
            [self.time[self.cadence], self.time[self.cadence]])
        self.update_bkg()
        self.update_lc()
        self.update_lcbkg()
        self.fig.canvas.draw()


def GetData(fitsfile, EPIC, campaign, clobber=False,
            saturation_tolerance=-0.1,
            bad_bits=[1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17],
            get_hires=False, get_nearby=False,
            aperture=None, **kwargs):
    '''
    Returns a :py:obj:`DataContainer` instance with the
    raw data for the target.

    :param str fitsfile: The full raw target pixel file path
    :param bool clobber: Overwrite existing files? Default :py:obj:`False`
    :param float saturation_tolerance: Target is considered saturated \
           if flux is within this fraction of the pixel well depth. \
           Default -0.1
    :param array_like bad_bits: Flagged :py:obj`QUALITY` bits to consider \
           outliers when computing the model. \
           Default `[1,2,3,4,5,6,7,8,9,11,12,13,14,16,17]`
    :param bool get_hires: Download a high resolution image of the target? \
           Default :py:obj:`True`
    :param bool get_nearby: Retrieve location of nearby sources? \
           Default :py:obj:`True`

    '''

    # Get the npz file name
    filename = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % campaign,
                            ('%09d' % EPIC)[:4] +
                            '00000', ('%09d' % EPIC)[4:],
                            'data.npz')

    # Create the dir
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

    # Check for saved data
    if not os.path.exists(filename) or clobber:

        log.info("Fetching data for target...")

        # Load the tpf
        with pyfits.open(fitsfile) as f:
            qdata = f[1].data

        # Get the header info
        fitsheader = [pyfits.getheader(fitsfile, 0).cards,
                      pyfits.getheader(fitsfile, 1).cards,
                      pyfits.getheader(fitsfile, 2).cards]

        # Get a hi res image of the target
        if get_hires:
            try:
                hires = GetHiResImage(EPIC)
            except ValueError:
                hires = None
        else:
            hires = None

        # Get nearby sources
        if get_nearby:
            try:
                nearby = GetSources(EPIC)
            except ValueError:
                nearby = []
        else:
            nearby = []

        # Get the arrays
        cadn = np.array(qdata.field('CADENCENO'), dtype='int32')
        time = np.array(qdata.field('TIME'), dtype='float64')
        fpix = np.array(qdata.field('FLUX'), dtype='float64')
        fpix_err = np.array(qdata.field('FLUX_ERR'), dtype='float64')
        qual = np.array(qdata.field('QUALITY'), dtype=int)

        # Get rid of NaNs in the time array by interpolating
        naninds = np.where(np.isnan(time))
        time = Interpolate(np.arange(0, len(time)), naninds, time)

        # Get the motion vectors (if available!)
        pc1 = np.array(qdata.field('POS_CORR1'), dtype='float64')
        pc2 = np.array(qdata.field('POS_CORR2'), dtype='float64')
        if not np.all(np.isnan(pc1)) and not np.all(np.isnan(pc2)):
            pc1 = Interpolate(time, np.where(np.isnan(pc1)), pc1)
            pc2 = Interpolate(time, np.where(np.isnan(pc2)), pc2)
        else:
            pc1 = None
            pc2 = None

        # Get the static pixel images for plotting
        pixel_images = [fpix[0], fpix[len(fpix) // 2], fpix[len(fpix) - 1]]

        # Get the aperture interactively
        if aperture is None:
            aperture = ApertureSelector(time[::10], fpix[::10],
                                        title='EPIC %d' % EPIC).aperture
        if np.sum(aperture) == 0:
            raise ValueError("Empty aperture!")

        # Atomically write to disk.
        # http://stackoverflow.com/questions/2333872/
        # atomic-writing-to-file-with-python
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        f = NamedTemporaryFile("wb", delete=False)
        np.savez_compressed(f, cadn=cadn, time=time, fpix=fpix,
                            fpix_err=fpix_err,
                            qual=qual, aperture=aperture,
                            pc1=pc1, pc2=pc2, fitsheader=fitsheader,
                            pixel_images=pixel_images, nearby=nearby,
                            hires=hires)
        f.flush()
        os.fsync(f.fileno())
        f.close()
        shutil.move(f.name, filename)

    # Load
    data = np.load(filename)
    aperture = data['aperture'][()]
    pixel_images = data['pixel_images']
    nearby = data['nearby'][()]
    hires = data['hires'][()]
    fitsheader = data['fitsheader']
    cadn = data['cadn']
    time = data['time']
    fpix = data['fpix']
    fpix_err = data['fpix_err']
    qual = data['qual']
    pc1 = data['pc1']
    pc2 = data['pc2']

    # Compute the saturation flux and the 97.5th percentile
    # flux in each pixel of the aperture. We're going
    # to compare these to decide if the star is saturated.
    satflx = SaturationFlux(EPIC, campaign=campaign) * \
        (1. + saturation_tolerance)
    f97 = np.zeros((fpix.shape[1], fpix.shape[2]))
    for i in range(fpix.shape[1]):
        for j in range(fpix.shape[2]):
            if aperture[i, j]:
                # Let's remove NaNs...
                tmp = np.delete(fpix[:, i, j], np.where(
                    np.isnan(fpix[:, i, j])))
                # ... and really bad outliers...
                if len(tmp):
                    f = SavGol(tmp)
                    med = np.nanmedian(f)
                    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
                    bad = np.where((f > med + 10. * MAD) |
                                   (f < med - 10. * MAD))[0]
                    np.delete(tmp, bad)
                    # ... so we can compute the 97.5th percentile flux
                    i97 = int(0.975 * len(tmp))
                    tmp = tmp[np.argsort(tmp)[i97]]
                    f97[i, j] = tmp

    # Check if any of the pixels are actually saturated
    if np.nanmax(f97) <= satflx:
        log.info("No saturated columns detected.")
        saturated = False
        aperture[np.isnan(fpix[0])] = 0
        ap = np.where(aperture & 1)
        fpix2D = np.array([f[ap] for f in fpix], dtype='float64')
        fpix_err2D = np.array([p[ap] for p in fpix_err], dtype='float64')
    else:
        # We need to collapse the saturated columns
        saturated = True
        ncol = 0
        fpixnew = []
        ferrnew = []
        for j in range(aperture.shape[1]):
            if np.any(f97[:, j] > satflx):
                marked = False
                collapsed = np.zeros(len(fpix[:, 0, 0]))
                collapsed_err2 = np.zeros(len(fpix[:, 0, 0]))
                for i in range(aperture.shape[0]):
                    if aperture[i, j]:
                        if not marked:
                            aperture[i, j] = AP_COLLAPSED_PIXEL
                            marked = True
                        else:
                            aperture[i, j] = AP_SATURATED_PIXEL
                        collapsed += fpix[:, i, j]
                        collapsed_err2 += fpix_err[:, i, j] ** 2
                if np.any(collapsed):
                    fpixnew.append(collapsed)
                    ferrnew.append(np.sqrt(collapsed_err2))
                    ncol += 1
            else:
                for i in range(aperture.shape[0]):
                    if aperture[i, j]:
                        fpixnew.append(fpix[:, i, j])
                        ferrnew.append(fpix_err[:, i, j])
        fpix2D = np.array(fpixnew).T
        fpix_err2D = np.array(ferrnew).T
        log.info("Collapsed %d saturated column(s)." % ncol)

    # Compute the background
    binds = np.where(aperture ^ 1)
    if RemoveBackground(EPIC, campaign=campaign) and (len(binds[0]) > 0):
        bkg = np.nanmedian(np.array([f[binds]
                                     for f in fpix], dtype='float64'), axis=1)
        # Uncertainty of the median:
        # http://davidmlane.com/hyperstat/A106993.html
        bkg_err = 1.253 * np.nanmedian(np.array([e[binds] for e in fpix_err],
                                       dtype='float64'), axis=1) \
            / np.sqrt(len(binds[0]))
        bkg = bkg.reshape(-1, 1)
        bkg_err = bkg_err.reshape(-1, 1)
    else:
        bkg = 0.
        bkg_err = 0.

    # Make everything 2D and remove the background
    fpix = fpix2D - bkg
    fpix_err = np.sqrt(fpix_err2D ** 2 + bkg_err ** 2)
    flux = np.sum(fpix, axis=1)

    # Get NaN data points
    nanmask = np.where(np.isnan(flux) | (flux == 0))[0]

    # Get flagged data points -- we won't train our model on them
    badmask = []
    for b in bad_bits:
        badmask += list(np.where(qual & 2 ** (b - 1))[0])

    # Flag >10 sigma outliers -- same thing.
    tmpmask = np.array(list(set(np.concatenate([badmask, nanmask]))))
    t = np.delete(time, tmpmask)
    f = np.delete(flux, tmpmask)
    f = SavGol(f)
    med = np.nanmedian(f)
    MAD = 1.4826 * np.nanmedian(np.abs(f - med))
    bad = np.where((f > med + 10. * MAD) | (f < med - 10. * MAD))[0]
    badmask.extend([np.argmax(time == t[i]) for i in bad])

    # Campaign 2 hack: the first day or two are screwed up
    if campaign == 2:
        badmask.extend(np.where(time < 2061.5)[0])

    # Finalize the mask
    badmask = np.array(sorted(list(set(badmask))))

    # Interpolate the nans
    fpix = Interpolate(time, nanmask, fpix)
    fpix_err = Interpolate(time, nanmask, fpix_err)

    # Return
    data = DataContainer()
    data.ID = EPIC
    data.campaign = campaign
    data.cadn = cadn
    data.time = time
    data.fpix = fpix
    data.fpix_err = fpix_err
    data.nanmask = nanmask
    data.badmask = badmask
    data.aperture = aperture
    data.aperture_name = 'custom'
    data.apertures = dict(custom=aperture)
    data.quality = qual
    data.Xpos = pc1
    data.Ypos = pc2
    data.meta = fitsheader
    data.mag = fitsheader[0]['KEPMAG'][1]
    if type(data.mag) is pyfits.card.Undefined:
        data.mag = np.nan
    data.pixel_images = pixel_images
    data.nearby = nearby
    data.hires = hires
    data.saturated = saturated
    data.bkg = bkg

    return data
