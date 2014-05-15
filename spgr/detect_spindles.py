from logging import getLogger, INFO, DEBUG
lg = getLogger('spgr')

from copy import deepcopy
from glob import glob
from os.path import join
from pickle import load

from numpy import asarray, hstack

from phypno.detect import DetectSpindle

from lsf import map_lsf

from .read_data import DATA_DIR, REC_FOLDER, STAGES

STAGE = 'sleep'

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


def calc_spindle_values(subj=None, detsp=None, ref_to_avg=None, lsf=True):
    """Detect spindles one channel in parallel with lsf.

    Parameters
    ----------
    subj : str
        code of patient
    detsp : instance of DetectSpindle
        detection parameters as DetectSpindle
    ref_to_avg : bool
        if true, use the data that was referenced to average.
    lsf : bool
        run on lsf or as normal loop

    Returns
    -------
    dict
        with 'spindles' (list of detected spindles), 'chan' (channel labels),
        'mean' (mean of values used for detection)

    """
    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    if ref_to_avg:
        data_file = glob(join(subj_dir, '*_' + STAGE + '_avg.pkl'))[0]
    else:
        data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]

    lg.info('Subj %s, reading data: %s', subj, data_file)
    with open(data_file, 'rb') as f:
        data = load(f)

    if lsf:
        lg.info('Subj %s, ready to submit %d jobs', subj,
                data.number_of('chan')[0])
        all_sp = map_lsf(det_sp_in_one_chan, get_one_chan(data),
                         queue='short',
                         variables={'detsp': detsp})
    else:
        all_sp = []
        for one_chan in get_one_chan(data):
            spindles = detsp(data)
            all_sp.append(spindles)

    spindles = [item for sublist in all_sp for item in sublist.spindle]
    spindles = sorted(spindles, key=lambda x:x['start_time'])

    sp = {'spindles': spindles,
          'chan': hstack([x.chan_name for x in all_sp]),
          'mean': hstack([x.mean for x in all_sp]),
          'std': hstack([x.std for x in all_sp]),
          'det_value': hstack([x.det_value for x in all_sp]),
          'sel_value': hstack([x.sel_value for x in all_sp]),
          }

    return sp

def det_sp_in_one_epoch(one_trial_data=None):
    from phypno.trans import Select, Math, Montage

    calc_std = Math(operator_name='std', axis='time')
    std_per_chan = calc_std(one_trial_data)
    good_chan = where((std_per_chan(trial=0) > .001) & (std_per_chan(trial=0) < thresh))[0]
    normal_chan = Select(chan=one_trial_data.axis['chan'][0][good_chan])

    if ref_to_avg:
        reref = Montage(ref_to_avg=True)
        one_trial_data = reref(one_trial_data)

    one_trial_data = normal_chan(one_trial_data)

    spindles = detsp(one_trial_data)
    return spindles.spindle


def calc_spindle_values_one_epoch(subj=None, detection_options=None, ref_to_avg=None):

    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]  # multiple files

    with open(data_file, 'rb') as f:
        data = load(f)

    detsp = DetectSpindle(**detection_options)

    all_sp = map_lsf(det_sp_in_one_epoch, iter(data),
                     queue='short',
                     variables={'detsp': detsp,
                                'thresh': 300,
                                'ref_to_avg': ref_to_avg},
                     imports={'numpy': 'where'})

    spindles = [x for sp in all_sp for x in sp]
    return spindles
