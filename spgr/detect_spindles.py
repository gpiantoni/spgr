from copy import deepcopy
from logging import getLogger
from multiprocessing import Pool
from pickle import load, dump

from numpy import asarray, hstack
from phypno.detect import DetectSpindle
from phypno.graphoelement import Spindles

from .constants import DATA_PATH, SPINDLE_FOLDER
from .read_data import get_data

lg = getLogger(__name__)

try:
    from lsf import map_lsf
    FORCE_LOCAL = False
except ImportError:
    FORCE_LOCAL = True
    lg.info('Could not import LSF, running local jobs only')


def get_spindles(subj, method='Nir2011', frequency=(None, None),
                 duration=(None, None), reref=None, resample_freq=None,
                 hp_filter=None, lp_filter=None, chan_type=('grid', )):

    spindle_name = ('spindles_{subj}_{method}_{frequency[0]}-{frequency[1]}Hz_'
                    '{duration[0]}-{duration[1]}s_{reref}_{resample_freq}_'
                    '{hp_filter}-{lp_filter}Hz.pkl'.format(**locals()))

    subj_dir = DATA_PATH / subj / SPINDLE_FOLDER
    if not subj_dir.exists():
        subj_dir.mkdir()
    spindle_file = subj_dir / spindle_name
    if spindle_file.exists():
        with open(str(spindle_file), 'rb') as f:
            spindles = load(f)

    else:
        detsp = DetectSpindle(method=method, frequency=frequency,
                              duration=duration)
        data = get_data(subj, 'sleep', chan_type, reref=reref,
                        resample_freq=resample_freq, hp_filter=hp_filter,
                        lp_filter=lp_filter)

        spindles = calc_spindle_values(data, detsp, parallel='')
        spindles = spindles(lambda x: frequency[0] <= x['peak_freq'] <= frequency[1])

        with open(str(spindle_file), 'wb') as f:
            dump(spindles, f)

    return spindles


def get_one_chan(data):
    """Generator that returns one channel at the time.

    Parameters
    ----------
    data : instance of DataTime
        recordings with multiple channels

    Returns
    -------
    instance of DataTime
        recording of only one channel.

    Notes
    -----
    This could also work but it's very slow:

    >>> def get_one_chan(data):
    >>>     for i, chan in enumerate(data.axis['chan'][0]):
    >>>         sel_one_chan = Select(chan=(chan, ))
    >>>         yield sel_one_chan(data)

    """
    for i_chan in range(data.number_of('chan')[0]):
        one_chan = deepcopy(data)
        for i_trl in range(data.number_of('trial')):
            one_chan.axis['chan'][i_trl] = asarray((data.axis['chan'][i_trl][i_chan], ))
            one_chan.data[i_trl] = data.data[i_trl][(i_chan, ), :]

        yield one_chan


def det_sp_in_one_chan(data):
    """Detect spindles in one channel. This is a convenience function for lsf.

    Parameters
    ----------
    data : instance of DataTime
        data of only one channel

    Returns
    -------
    instance of Spindles
        info about spindles, data and thresholds

    """
    spindles = detsp(data)
    return spindles


def calc_spindle_values(data, detsp=None, parallel='lsf'):
    """Detect spindles one channel in parallel with lsf.

    Parameters
    ----------
    data : instance of DataTime
        data with the recordings
    detsp : instance of DetectSpindle
        detection parameters as DetectSpindle
    parallel : str
        run on lsf ('lsf'), as parallel ('pool') or as normal loop ('map')

    Returns
    -------
    dict
        with 'spindles' (list of detected spindles), 'chan' (channel labels),
        'mean' (mean of values used for detection)

    """
    if parallel == 'lsf' and not FORCE_LOCAL:
        all_sp = map_lsf(det_sp_in_one_chan, get_one_chan(data),
                         queue='short',
                         variables={'detsp': detsp})
    elif parallel == 'pool':
        with Pool() as p:
            all_sp = p.map(detsp, get_one_chan(data))
    else:
        all_sp = list(map(detsp, get_one_chan(data)))

    spindles = [item for sublist in all_sp for item in sublist.spindle]
    spindles = sorted(spindles, key=lambda x: x['start_time'])

    sp = Spindles()
    sp.spindle = spindles
    sp.chan_name = hstack([x.chan_name for x in all_sp])
    sp.mean = hstack([x.mean for x in all_sp])
    sp.std = hstack([x.std for x in all_sp])
    sp.det_value = hstack([x.det_value for x in all_sp])
    sp.sel_value = hstack([x.sel_value for x in all_sp])

    return sp
