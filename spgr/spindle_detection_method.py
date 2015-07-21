from copy import deepcopy
from numpy import ones, NaN, swapaxes, where
from scipy.signal import fftconvolve
from pyqtgraph import mkPen
from phypno.trans import Math, Filter

from phypno.detect import DetectSpindle
from phypno.trans import Select
from phypno.viz import Viz1

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        PLOT_COLOR,
                        SPINDLE_OPTIONS,
                        )
from .detect_spindles import get_spindles
from .read_data import get_data

from .log import with_log


PLOT_SIZE = 960, 240

SUBJ = 'EM09'
chan_name = 'GR50'
start_good_time = 44850
end_good_time = 44880


@with_log
def Spindle_Detection_Method(lg, images_dir):

    lg.info('## Method description')
    lg.info('Spindles were detected on each electrode independently, using '
            'previously reported methods with stringent criteria ({method}).'
            'The raw signal was filtered between {freq0} and {freq1} Hz.'
            'Spindle duration had to be between {dur0} and {dur1} s.'
            ''.format(method=SPINDLE_OPTIONS['method'],
                      freq0=SPINDLE_OPTIONS['frequency'][0],
                      freq1=SPINDLE_OPTIONS['frequency'][1],
                      dur0=SPINDLE_OPTIONS['duration'][0],
                      dur1=SPINDLE_OPTIONS['duration'][1]))

    lg.info('Example from {subj} on chan {chan} at {time} s'
            ''.format(subj=SUBJ, chan=chan_name, time=start_good_time))

    data = get_data(SUBJ, 'sleep', CHAN_TYPE, **DATA_OPTIONS)

    good_trial = [i for i, trl in enumerate(data)
                  if start_good_time // 30 * 30 == trl.axis['time'][0][0]][0]

    sel = Select(trial=(good_trial, ), chan=(chan_name, ))
    sel_data = sel(data)

    v = Viz1(color=PLOT_COLOR)
    v.size = PLOT_SIZE
    v.add_data(sel_data, limits_x=(start_good_time, end_good_time),
               limits_y=(-200, 200))
    v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')

    png_file = str(images_dir.joinpath('raw.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('Raw EEG', png_file))

    detsp = DetectSpindle(method=SPINDLE_OPTIONS['method'],
                          frequency=SPINDLE_OPTIONS['frequency'],
                          duration=SPINDLE_OPTIONS['duration'])
    sp = detsp(sel_data)

    filt = Filter(low_cut=detsp.det_butter['freq'][0],
                  high_cut=detsp.det_butter['freq'][1],
                  order=detsp.det_butter['order'], s_freq=data.s_freq)
    x_filt = filt(sel_data)

    def _move_avg(dat, axis):
        """Function to compute moving average."""
        flat = ones(detsp.moving_avg['dur'] * data.s_freq)

        dat = swapaxes(dat, axis, -1)
        for i, x in enumerate(dat):
            dat[i, :] = fftconvolve(x, flat / sum(flat), mode='same')
        return swapaxes(dat, -1, axis)

    hilbert_abs = Math(operator_name=('hilbert', 'abs'), axis='time')
    moving_avg = Math(operator=_move_avg, axis='time')

    dat_det = moving_avg(hilbert_abs(x_filt))

    v = Viz1(color=PLOT_COLOR)
    v.size = PLOT_SIZE
    v.add_data(x_filt, color=mkPen('b', width=1),
               limits_x=(start_good_time, end_good_time), limits_y=(-100, 100))
    v.add_data(dat_det, color=mkPen('r', width=4))
    v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')

    png_file = str(images_dir.joinpath('envelope.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('Envelope', png_file))

    t = dat_det.axis['time'][0]
    good_sp = [x for x in sp.spindle
               if x['start_time'] < t[-1] and x['end_time'] > t[0]]

    detected_spindles = deepcopy(dat_det)

    for one_sp in good_sp:
        i0 = where(t == one_sp['start_time'])[0]
        i1 = where(t == one_sp['end_time'])[0]
        detected_spindles.data[0][0, i0:i1] = NaN

    v = Viz1(color=PLOT_COLOR)
    v.size = PLOT_SIZE
    v.add_data(sel_data, limits_x=(start_good_time, end_good_time),
               limits_y=(-200, 200))
    v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')
    v.add_graphoelement(sp)
    png_file = str(images_dir.joinpath('detected.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('Detected', png_file))

    for ref in ('avg', 15):
        for subj in HEMI_SUBJ:
            print(subj)
            sp = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref,
                              **SPINDLE_OPTIONS)
