#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`slurm.py` - SLURM cluster routines
-------------------------------------------

Routines for submitting batch jobs to a cluster.

'''

from __future__ import division, print_function, absolute_import, \
     unicode_literals
from .utils import *
from .k2 import GetData, FITSFile
from ...config import EVEREST_SRC, EVEREST_DAT, EVEREST_DEV
from ...utils import ExceptionHook, FunctionWrapper
from ...pool import Pool
import os
import sys
import subprocess
import numpy as np
import pickle
import traceback
import glob
import logging
log = logging.getLogger(__name__)

# Constants
RED = '\033[0;31m'
GREEN = '\033[0;32m'
BLACK = '\033[0m'
BLUE = '\033[0;34m'


def Download(campaign=0, queue='cca', email=None, walltime=8, **kwargs):
    '''
    Submits a cluster job to the build queue to download all TPFs for a given
    campaign.

    :param int campaign: The `K2` campaign to run
    :param str queue: The name of the queue to submit to. Default `build`
    :param str email: The email to send job status notifications to. \
           Default `None`
    :param int walltime: The number of hours to request. Default `8`

    '''

    # Figure out the subcampaign
    if type(campaign) is int:
        subcampaign = -1
    elif type(campaign) is float:
        x, y = divmod(campaign, 1)
        campaign = int(x)
        subcampaign = round(y * 10)
    # Submit the cluster job
    slurmfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'download.sh')
    str_w = '--time=%d:00:00' % walltime
    str_v = '--export=ALL,EVERESTDAT=%s,CAMPAIGN=%d,SUBCAMPAIGN=%d' % (
        EVEREST_DAT, campaign, subcampaign)
    if subcampaign == -1:
        str_name = 'dwn%02d' % campaign
    else:
        str_name = 'dwn%02d.%d' % (campaign, subcampaign)
    str_out = "--output=%s" % os.path.join(EVEREST_DAT, 'k2', str_name + '.log')
    sbatch_args = ['sbatch',
                   "--partition=%s" % queue,
                   str_v,
                   str_out,
                   "--job-name=%s" % str_name,
                   str_w]
    if email is not None:
        sbatch_args.append(['--mail-user=%s' % email, '--mail-type=END,FAIL'])
    sbatch_args.append(slurmfile)
    # Now we submit the job
    print("Submitting the job...")
    print(" ".join(sbatch_args))
    subprocess.call(sbatch_args)


def _Download(campaign, subcampaign):
    '''
    Download all stars from a given campaign. This is
    called from ``missions/k2/download.sh``

    '''

    # Are we doing a subcampaign?
    if subcampaign != -1:
        campaign = campaign + 0.1 * subcampaign
    # Get all star IDs for this campaign
    stars = [s[0] for s in GetK2Campaign(campaign)]
    nstars = len(stars)
    # Download the TPF data for each one
    for i, EPIC in enumerate(stars):
        print("Downloading data for EPIC %d (%d/%d)..." %
              (EPIC, i + 1, nstars))
        if not os.path.exists(os.path.join(EVEREST_DAT, 'k2',
                                           'c%02d' % int(campaign),
                                           ('%09d' % EPIC)[:4] + '00000',
                                           ('%09d' % EPIC)[4:],
                                           'data.npz')):
            try:
                GetData(EPIC, season=campaign, download_only=True)
            except KeyboardInterrupt:
                sys.exit()
            except:
                # Some targets could be corrupted...
                print("ERROR downloading EPIC %d." % EPIC)
                exctype, value, tb = sys.exc_info()
                for line in traceback.format_exception_only(exctype, value):
                    ln = line.replace('\n', '')
                    print(ln)
                continue


def Run(campaign=0, EPIC=None, nodes=5, ppn=None, walltime=100,
        mpn=None, email=None, queue='cca', **kwargs):
    '''
    Submits a cluster job to compute and plot data for all targets
    in a given campaign.

    :param campaign: The K2 campaign number. If this is an :py:class:`int`, \
           returns all targets in that campaign. If a :py:class:`float` \
           in the form `X.Y`, runs the `Y^th` decile of campaign `X`.
    :param str queue: The queue to submit to. Default `None` (default queue)
    :param str email: The email to send job status notifications to. \
           Default `None`
    :param int walltime: The number of hours to request. Default `100`
    :param int nodes: The number of nodes to request. Default `5`
    :param int ppn: The number of processors per node to request. Default `12`
    :param int mpn: Memory per node in gb to request. Default no setting.

    '''

    # Figure out the subcampaign
    if type(campaign) is int:
        subcampaign = -1
    elif type(campaign) is float:
        x, y = divmod(campaign, 1)
        campaign = int(x)
        subcampaign = round(y * 10)

    # Convert kwargs to string. This is really hacky. Pickle creates an array
    # of bytes, which we must convert into a regular string to pass to the pbs
    # script and then back into python. Decoding the bytes isn't enough, since
    # we have pesky escaped characters such as newlines that don't behave well
    # when passing this string around. My braindead hack is to replace newlines
    # with '%%%', then undo the replacement when reading the kwargs. This works
    # for most cases, but sometimes pickle creates a byte array that can't be
    # decoded into utf-8; this happens when trying to pass numpy arrays around,
    # for instance. This needs to be fixed in the future, but for now we'll
    # restrict the kwargs to be ints, floats, lists, and strings.
    try:
        strkwargs = pickle.dumps(kwargs, 0).decode(
            'utf-8').replace('\n', '%%%')
    except UnicodeDecodeError:
        raise ValueError('Unable to pickle `kwargs`. Currently the `kwargs` ' +
                         'values may only be `int`s, `float`s, `string`s, ' +
                         '`bool`s, or lists of these.')

    # Number of tasks
    if ppn is not None:
        tasks = nodes * ppn
    else:
        tasks = nodes * 40

    # Submit the cluster job
    slurmfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'run.sh')
    str_w = '--time=%d:00:00' % walltime
    str_v = "--export=ALL,EVERESTDAT=%s,NODES=%d," % (EVEREST_DAT, nodes) + \
            "EPIC=%d," % (0 if EPIC is None else EPIC) + \
            "CAMPAIGN=%d,SUBCAMPAIGN=%d,STRKWARGS='%s',TASKS=%d" % \
            (campaign, subcampaign, strkwargs, tasks)

    if EPIC is None:
        if subcampaign == -1:
            str_name = 'c%02d' % campaign
        else:
            str_name = 'c%02d.%d' % (campaign, subcampaign)
    else:
        str_name = 'EPIC%d' % EPIC
    str_out = "--output=%s" % os.path.join(EVEREST_DAT, 'k2', str_name + '.log')
    sbatch_args = ['sbatch',
                   "--partition=%s" % queue,
                   str_v,
                   str_out,
                   "--job-name=%s" % str_name,
                   "--nodes=%d" % nodes,
                   str_w]
    if email is not None:
        sbatch_args.append(['--mail-user=%s' % email, '--mail-type=END,FAIL'])
    if mpn is not None:
        sbatch_args.append('--mem=%dG' % mpn)
    if ppn is not None:
        sbatch_args.append('--ntasks-per-node=%d' % ppn)
    if queue == 'preempt':
        sbatch_args.append('--qos=preempt')
    sbatch_args.append(slurmfile)
    
    # Now we submit the job
    print("Submitting the job...")
    print(" ".join(sbatch_args))
    subprocess.call(sbatch_args)


def _Run(campaign, subcampaign, epic, strkwargs):
    '''
    The actual function that runs a given campaign; this must
    be called from ``missions/k2/run.sh``.

    '''

    # Get kwargs from string
    kwargs = pickle.loads(strkwargs.replace('%%%', '\n').encode('utf-8'))

    # Check the cadence
    cadence = kwargs.get('cadence', 'lc')

    # Model wrapper
    m = FunctionWrapper(EverestModel, season=campaign, **kwargs)

    # Set up our custom exception handler
    sys.excepthook = ExceptionHook

    # Are we running a campaign or a single target?
    if epic == 0:

        # Initialize our multiprocessing pool
        with Pool() as pool:
            # Are we doing a subcampaign?
            if subcampaign != -1:
                campaign = campaign + 0.1 * subcampaign
            # Get all the stars
            stars = GetK2Campaign(campaign, epics_only=True, cadence=cadence)
            # Run
            pool.map(m, stars)

    else:

        m(epic)


def Publish(campaign=0, EPIC=None, nodes=5, ppn=None, walltime=100,
            mpn=None, email=None, queue='cca', **kwargs):
    '''
    Submits a cluster job to generate the FITS files for publication.
    Make sure to run :py:func:`everest.k2.GetCBVs` for this campaign
    beforehand.

    :param campaign: The K2 campaign number. If this is an :py:class:`int`, \
           returns all targets in that campaign. If a :py:class:`float` \
           in the form `X.Y`, runs the `Y^th` decile of campaign `X`.
    :param str queue: The queue to submit to. Default `None` (default queue)
    :param str email: The email to send job status notifications to. \
           Default `None`
    :param int walltime: The number of hours to request. Default `100`
    :param int nodes: The number of nodes to request. Default `5`
    :param int ppn: The number of processors per node to request. Default `12`
    :param int mpn: Memory per node in gb to request. Default no setting.

    '''

    # Figure out the subcampaign
    if type(campaign) is int:
        subcampaign = -1
    elif type(campaign) is float:
        x, y = divmod(campaign, 1)
        campaign = int(x)
        subcampaign = round(y * 10)

    # Convert kwargs to string. This is really hacky. Pickle creates an array
    # of bytes, which we must convert into a regular string to pass to the pbs
    # script and then back into python. Decoding the bytes isn't enough, since
    # we have pesky escaped characters such as newlines that don't behave well
    # when passing this string around. My braindead hack is to replace newlines
    # with '%%%', then undo the replacement when reading the kwargs. This works
    # for most cases, but sometimes pickle creates a byte array that can't be
    # decoded into utf-8; this happens when trying to pass numpy arrays around,
    # for instance. This needs to be fixed in the future, but for now we'll
    # restrict the kwargs to be ints, floats, lists, and strings.
    try:
        strkwargs = pickle.dumps(kwargs, 0).decode(
            'utf-8').replace('\n', '%%%')
    except UnicodeDecodeError:
        raise ValueError('Unable to pickle `kwargs`. Currently the `kwargs` ' +
                         'values may only be `int`s, `float`s, `string`s, ' +
                         '`bool`s, or lists of these.')

    # Number of tasks
    if ppn is not None:
        tasks = nodes * ppn
    else:
        tasks = nodes * 40

    # Submit the cluster job
    slurmfile = os.path.join(EVEREST_SRC, 'missions', 'k2', 'publish.sh')
    str_w = '--time=%d:00:00' % walltime
    str_v = "--export=ALL,EVERESTDAT=%s,NODES=%d," % (EVEREST_DAT, nodes) + \
            "CAMPAIGN=%d,SUBCAMPAIGN=%d,STRKWARGS='%s',TASKS=%d" % \
            (campaign, subcampaign, strkwargs, tasks)

    if subcampaign == -1:
        str_name = 'c%02d' % campaign
    else:
        str_name = 'c%02d.%d' % (campaign, subcampaign)
    str_out = "--output=%s" % os.path.join(EVEREST_DAT, 'k2', str_name + '.log')
    sbatch_args = ['sbatch',
                   "--partition=%s" % queue,
                   str_v,
                   str_out,
                   "--job-name=%s" % str_name,
                   "--nodes=%d" % nodes,
                   str_w]
    if email is not None:
        sbatch_args.append(['--mail-user=%s' % email, '--mail-type=END,FAIL'])
    if mpn is not None:
        sbatch_args.append('--mem=%dG' % mpn)
    if ppn is not None:
        sbatch_args.append('--ntasks-per-node=%d' % ppn)
    if queue == 'preempt':
        sbatch_args.append('--qos=preempt')
    sbatch_args.append(slurmfile)
    
    # Now we submit the job
    print("Submitting the job...")
    print(" ".join(sbatch_args))
    subprocess.call(sbatch_args)


def _Publish(campaign, subcampaign, strkwargs):
    '''
    The actual function that publishes a given campaign; this must
    be called from ``missions/k2/publish.sh``.

    '''

    # Get kwargs from string
    kwargs = pickle.loads(strkwargs.replace('%%%', '\n').encode('utf-8'))

    # Check the cadence
    cadence = kwargs.get('cadence', 'lc')

    # Model wrapper
    m = FunctionWrapper(EverestModel, season=campaign, publish=True, **kwargs)

    # Set up our custom exception handler
    sys.excepthook = ExceptionHook

    # Initialize our multiprocessing pool
    with Pool() as pool:
        # Are we doing a subcampaign?
        if subcampaign != -1:
            campaign = campaign + 0.1 * subcampaign
        # Get all the stars
        stars = GetK2Campaign(campaign, epics_only=True, cadence=cadence)

        # Run
        pool.map(m, stars)


def Status(season=range(18), model='nPLD', purge=False, injection=False,
           cadence='lc', purge_hard=False, **kwargs):
    '''
    Shows the progress of the de-trending runs for the specified campaign(s).

    '''

    # Mission compatibility
    campaign = season

    # Injection?
    if injection:
        return InjectionStatus(campaign=campaign, model=model,
                               purge=purge, **kwargs)

    # Cadence
    if cadence == 'sc':
        model = '%s.sc' % model

    if not hasattr(campaign, '__len__'):
        if type(campaign) is int:
            # Return the subcampaigns
            all_stars = [s for s in GetK2Campaign(
                campaign, split=True, epics_only=True, cadence=cadence)]
            campaign = [campaign + 0.1 * n for n in range(10)]
        else:
            all_stars = [[s for s in GetK2Campaign(
                campaign, epics_only=True, cadence=cadence)]]
            campaign = [campaign]
    else:
        all_stars = [[s for s in GetK2Campaign(
            c, epics_only=True, cadence=cadence)] for c in campaign]

    print("CAMP      TOTAL      DOWNLOADED    PROCESSED      FITS    ERRORS")
    print("----      -----      ----------    ---------      ----    ------")
    for c, stars in zip(campaign, all_stars):
        if len(stars) == 0:
            continue
        down = 0
        proc = 0
        err = 0
        fits = 0
        bad = []
        remain = []
        total = len(stars)
        if os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'c%02d' % c)):
            path = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % c)
            for folder in [f for f in os.listdir(path) if f.endswith('00000')]:
                for subfolder in os.listdir(os.path.join(path, folder)):
                    ID = int(folder[:4] + subfolder)
                    if ID in stars:
                        if os.path.exists(os.path.join(EVEREST_DAT,
                                                       'k2', 'c%02d' % c,
                                                       folder,
                                                       subfolder, 'data.npz')):
                            down += 1
                        if os.path.exists(os.path.join(EVEREST_DAT, 'k2',
                                                       'c%02d' % c, folder,
                                                       subfolder, FITSFile(
                                                         ID, c,
                                                         cadence=cadence))):
                            fits += 1
                        if os.path.exists(os.path.join(EVEREST_DAT, 'k2',
                                                       'c%02d' % c, folder,
                                                       subfolder,
                                                       model + '.npz')):
                            proc += 1
                        elif os.path.exists(os.path.join(EVEREST_DAT, 'k2',
                                                         'c%02d' % c, folder,
                                                         subfolder,
                                                         model + '.err')):
                            err += 1
                            bad.append(folder[:4] + subfolder)
                            if purge:
                                os.remove(os.path.join(
                                    EVEREST_DAT, 'k2', 'c%02d' % c,
                                    folder, subfolder, model + '.err'))
                            if purge_hard:
                                for file in glob.glob(os.path.join(
                                    EVEREST_DAT, 'k2', 'c%02d' % c,
                                    folder, subfolder, '*')):
                                    os.remove(file)
                        else:
                            remain.append(folder[:4] + subfolder)
        if proc == total:
            cc = ct = cp = ce = GREEN
            cd = BLACK if down < total else GREEN
        else:
            cc = BLACK
            ct = BLACK
            cd = BLACK if down < total else BLUE
            cp = BLACK if proc < down or proc == 0 else BLUE
            ce = RED if err > 0 else BLACK
        cf = BLACK if fits < total else GREEN
        if type(c) is int:
            print("%s{:>4d}   \033[0m%s{:>8d}\033[0m%s{:>16d}\033[0m%s{:>13d}\033[0m%s{:>10d}\033[0m%s{:>10d}\033[0m".format(c, total, down, proc, fits, err)
                  % (cc, ct, cd, cp, cf, ce))
        else:
            print("%s{:>4.1f}   \033[0m%s{:>8d}\033[0m%s{:>16d}\033[0m%s{:>13d}\033[0m%s{:>10d}\033[0m%s{:>10d}\033[0m".format(c, total, down, proc, fits, err)
                  % (cc, ct, cd, cp, cf, ce))
        if len(remain) <= 25 and len(remain) > 0 and len(campaign) == 1:
            remain.extend(["         "] * (4 - (len(remain) % 4)))
            print()
            for A, B, C, D in zip(remain[::4], remain[1::4],
                                  remain[2::4], remain[3::4]):
                if A == remain[0]:
                    print("REMAIN:  %s   %s   %s   %s" % (A, B, C, D))
                    print()
                else:
                    print("         %s   %s   %s   %s" % (A, B, C, D))
                    print()
        if len(bad) and len(campaign) == 1:
            bad.extend(["         "] * (4 - (len(bad) % 4)))
            print()
            for A, B, C, D in zip(bad[::4], bad[1::4], bad[2::4], bad[3::4]):
                if A == bad[0]:
                    print("ERRORS:  %s   %s   %s   %s" % (A, B, C, D))
                    print()
                else:
                    print("         %s   %s   %s   %s" % (A, B, C, D))
                    print()


def InjectionStatus(campaign=range(18), model='nPLD', purge=False,
                    depths=[0.01, 0.001, 0.0001], **kwargs):
    '''
    Shows the progress of the injection de-trending runs for
    the specified campaign(s).

    '''

    if not hasattr(campaign, '__len__'):
        if type(campaign) is int:
            # Return the subcampaigns
            all_stars = [s for s in GetK2Campaign(
                campaign, split=True, epics_only=True)]
            campaign = [campaign + 0.1 * n for n in range(10)]
        else:
            all_stars = [[s for s in GetK2Campaign(campaign, epics_only=True)]]
            campaign = [campaign]
    else:
        all_stars = [[s for s in GetK2Campaign(
            c, epics_only=True)] for c in campaign]
    print("CAMP      MASK       DEPTH     TOTAL      DONE     ERRORS")
    print("----      ----       -----     -----      ----     ------")
    for c, stars in zip(campaign, all_stars):
        if len(stars) == 0:
            continue
        done = [[0 for d in depths], [0 for d in depths]]
        err = [[0 for d in depths], [0 for d in depths]]
        total = len(stars)
        if os.path.exists(os.path.join(EVEREST_DAT, 'k2', 'c%02d' % c)):
            path = os.path.join(EVEREST_DAT, 'k2', 'c%02d' % c)
            for folder in os.listdir(path):
                for subfolder in os.listdir(os.path.join(path, folder)):
                    ID = int(folder[:4] + subfolder)
                    for m, mask in enumerate(['U', 'M']):
                        for d, depth in enumerate(depths):
                            if os.path.exists(
                                os.path.join(
                                    EVEREST_DAT, 'k2', 'c%02d' % c, folder,
                                    subfolder, '%s_Inject_%s%g.npz' %
                                    (model, mask, depth))):
                                done[m][d] += 1
                            elif os.path.exists(
                                    os.path.join(
                                        EVEREST_DAT, 'k2', 'c%02d' % c, folder,
                                        subfolder, '%s_Inject_%s%g.err' %
                                        (model, mask, depth))):
                                err[m][d] += 1
        for d, depth in enumerate(depths):
            for m, mask in enumerate(['F', 'T']):
                if done[m][d] == total:
                    color = GREEN
                else:
                    color = BLACK
                if err[m][d] > 0:
                    errcolor = RED
                else:
                    errcolor = ''
                if type(c) is int:
                    print("%s{:>4d}{:>8s}{:>14g}{:>10d}{:>10d}%s{:>9d}\033[0m".format(
                        c, mask, depth, total, done[m][d], err[m][d]) % (color, errcolor))
                else:
                    print("%s{:>4.1f}{:>8s}{:>14g}{:>10d}{:>10d}%s{:>9d}\033[0m".format(
                        c, mask, depth, total, done[m][d], err[m][d]) % (color, errcolor))


def EverestModel(ID, model='nPLD', publish=False, csv=False, **kwargs):
    '''
    A wrapper around an :py:obj:`everest` model for PBS runs.

    '''

    if model != 'Inject':
        from ... import detrender

        # HACK: We need to explicitly mask short cadence planets
        if kwargs.get('cadence', 'lc') == 'sc':
            EPIC, t0, period, duration = \
                np.loadtxt(os.path.join(EVEREST_SRC, 'missions', 'k2',
                                        'tables', 'scmasks.tsv'), unpack=True)
            if ID in EPIC and kwargs.get('planets', None) is None:
                ii = np.where(EPIC == ID)[0]
                planets = []
                for i in ii:
                    planets.append([t0[i], period[i], 1.25 * duration[i]])
                kwargs.update({'planets': planets})

        # Run the model
        m = getattr(detrender, model)(ID, **kwargs)

        # Publish?
        if publish:
            if csv:
                m.publish_csv()
            else:
                m.publish()

    else:
        from ...inject import Inject
        Inject(ID, **kwargs)
    return True
