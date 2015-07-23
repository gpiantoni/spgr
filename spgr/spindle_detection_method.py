from phypno.detect import DetectSpindle
from phypno.trans import Select
from phypno.viz import Viz1

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        PLOT_COLOR,
                        SPINDLE_OPTIONS,
                        RAW_LIMITS_Y,
                        )
from .detect_spindles import get_spindles
from .read_data import get_data

from .log import with_log


PLOT_SIZE = 960, 480

SUBJ = 'EM09'
chan_names = ['GR' + str(i) for i in range(28, 55)]
chan_names = ['GR35', 'GR50']
start_good_time = 44850
end_good_time = 44880


@with_log
def Spindle_Detection_Method(lg, images_dir):

    lg.info('## Method description')
    lg.info('Spindles were detected on each electrode independently, using '
            'previously reported methods with stringent criteria ({method}). '
            'The raw signal was filtered between {freq0} and {freq1} Hz. '
            'Spindle duration had to be between {dur0} and {dur1} s.'
            ''.format(method=SPINDLE_OPTIONS['method'],
                      freq0=SPINDLE_OPTIONS['frequency'][0],
                      freq1=SPINDLE_OPTIONS['frequency'][1],
                      dur0=SPINDLE_OPTIONS['duration'][0],
                      dur1=SPINDLE_OPTIONS['duration'][1]))

    lg.info('Example from {subj} at {time} s'
            ''.format(subj=SUBJ, time=start_good_time))

    data = get_data(SUBJ, 'sleep', CHAN_TYPE, **DATA_OPTIONS)

    good_trial = [i for i, trl in enumerate(data)
                  if start_good_time // 30 * 30 == trl.axis['time'][0][0]][0]

    sel = Select(trial=(good_trial, ), chan=chan_names)
    sel_data = sel(data)

    detsp = DetectSpindle(method=SPINDLE_OPTIONS['method'],
                          frequency=SPINDLE_OPTIONS['frequency'],
                          duration=SPINDLE_OPTIONS['duration'])
    sp = detsp(sel_data)

    v = Viz1(color=PLOT_COLOR)
    v.size = PLOT_SIZE
    v.add_data(sel_data, limits_x=(start_good_time, end_good_time),
               limits_y=RAW_LIMITS_Y)
    for p in v._plots.values():
        p.setLabels(left='amplitude (μV)',
                    bottom='time (s)')
    v.add_graphoelement(sp)
    png_file = str(images_dir.joinpath('detected.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('Detected', png_file))

    for ref in ('avg', 15):
        for subj in HEMI_SUBJ:
            print(subj)
            sp = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref,
                              **SPINDLE_OPTIONS)
