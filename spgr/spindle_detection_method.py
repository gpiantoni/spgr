from numpy import array

from phypno.detect import DetectSpindle
from phypno.trans import Select, Filter
from phypno.viz import Viz1

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        HIGHLIGHT_COLOR,
                        DPI,
                        SPINDLE_OPTIONS,
                        RAW_LIMITS_Y,
                        LABEL_FONT_SIZE,
                        TICKS_FONT_SIZE,
                        )
from .detect_spindles import get_spindles
from .read_data import get_data

from .log import with_log


PLOT_SIZE = 112, 50

SUBJ = 'EM09'
chan_names = ['GR' + str(i) for i in range(28, 55)]
chan_names = ['GR35', 'GR50']
start_good_time = 44850
end_good_time = 44880
HIGHLIGHT_FILTER = 1, 30


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

    filt = Filter(low_cut=HIGHLIGHT_FILTER[0],
                  high_cut=HIGHLIGHT_FILTER[1],
                  s_freq=sel_data.s_freq)
    sel_data = filt(sel_data)

    v = _make_cooccur_plot(sel_data, sp)
    png_file = str(images_dir.joinpath('detected.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('Detected', png_file))

    for ref in ('avg', 15):
        for subj in HEMI_SUBJ:
            print(subj)
            sp = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref,
                              **SPINDLE_OPTIONS)


def _make_cooccur_plot(sel_data, sp):

    v = Viz1(show=True, size_mm=PLOT_SIZE, dpi=DPI)
    v.add_data(sel_data,
               limits_y=[RAW_LIMITS_Y[0] * 2, RAW_LIMITS_Y[1]*2])
    v.add_graphoelement(sp, color=HIGHLIGHT_COLOR, height=400)

    def tick_frac():
        major_tick_fractions = array([0, 1])
        minor_tick_fractions = array([.3, .6])
        tick_labels = array(['0', '30'])
        return major_tick_fractions, minor_tick_fractions, tick_labels

    for plt in v._fig.plot_widgets:
        plt.xaxis.axis.ticker._get_tick_frac_labels = tick_frac
        plt.view.camera.set_range(x=(44851.5, 44878.5), y=(-200, 200))
        plt.xaxis.axis.major_tick_length = 10
        plt.xaxis.axis.minor_tick_length = 0
        plt.xaxis.axis._text.font_size = TICKS_FONT_SIZE

        plt.yaxis.axis._text.font_size = TICKS_FONT_SIZE
        plt.ylabel.text = 'amplitude (μV)'
        plt.ylabel._text_visual.font_size = LABEL_FONT_SIZE

        plt.margin = 25  # otherwise xtick label overlaps with border
        plt.view.border_color = 'w'

    plt.xlabel.text = 'time (s)'
    plt.xlabel._text_visual.font_size = LABEL_FONT_SIZE

    return v
