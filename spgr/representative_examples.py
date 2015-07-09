from numpy import where
from phypno.trans import Select
from phypno.viz import Viz1
from scipy.signal import periodogram

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        PLOT_COLOR,
                        SPINDLE_OPTIONS,
                        )
from .detect_spindles import get_spindles
from .read_data import get_data
from .spindle_source import get_chan_with_regions

from .log import with_log


REREF = 'avg'
PAD = 2
YLIM = 150
WINDOW = 1

BEST_SPINDLE = 'sigma_ratio'  # 'area_under_curve' or 'sigma_ratio'


@with_log
def Representative_Examples(lg, images_dir):
    """Show some good spindle examples.
    A better way of doing this is by doing fft at spindle, before and after and
    take the spindle with the best ratio.
    """
    lg.info('## Representative Spindles')

    lg.info('### Spindles in each brain region')

    all_pos = []
    all_title = []
    all_png_file = []

    for subj in HEMI_SUBJ:

        data = get_data(subj, 'sleep', CHAN_TYPE, reref='avg', **DATA_OPTIONS)
        best_spindles = find_best_spindles(subj, data)
        chan = data.attr['chan']

        for region, spindle in best_spindles.items():

            spindle_data = find_spindle_data(data, spindle)
            v = Viz1(color=PLOT_COLOR)
            v.add_data(spindle_data, limits_y=(-YLIM, YLIM))
            v.add_graphoelement([spindle, ])

            png_file = str(images_dir.joinpath('{}_{}.png'.format(region, subj)))
            v.save(png_file)

            y_pos = chan(lambda x: x.label == spindle['chan']).return_xyz()[0, 1]
            all_pos.append(y_pos)
            all_title.append('### Representative for {} ({})'.format(region,
                                                                     subj))
            all_png_file.append(png_file)

    # sort by anterior-posterior axis
    for _, title, png_file in sorted(zip(all_pos, all_title, all_png_file)):
        lg.info(title)
        lg.info('![]({})'.format(png_file))


def find_best_spindles(subj, data):

    spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
    chan = get_chan_with_regions(subj, REREF)

    best_spindles = {}

    for one_sp in spindles:
        sp_region = chan(lambda x: x.label == one_sp['chan']).return_attr('region')[0]

        if not sp_region.startswith('ctx-'):
            continue

        sp_region = sp_region[len('ctx-Xh-'):]

        # there are different ways to define the best spindle.
        # the quality of a spindle will be stored in goodness (the higher the better)
        if BEST_SPINDLE == 'area_under_curve':
            one_sp['goodness'] = one_sp['area_under_curve']
        elif BEST_SPINDLE == 'sigma_ratio':
            one_sp['goodness'] = _find_sigma_ratio(one_sp, data)

        if sp_region in best_spindles:
            if one_sp['goodness'] > best_spindles[sp_region]['goodness']:
                best_spindles[sp_region] = one_sp
        else:
            best_spindles[sp_region] = one_sp

    return best_spindles


def find_spindle_data(data, spindle):
    i_trial = _find_trial_with_spindle(data, spindle)

    sel = Select(time=(spindle['start_time'] - PAD,
                       spindle['end_time'] + PAD),
                 trial=(i_trial, ), chan=(spindle['chan'], ))
    return sel(data)


def _find_trial_with_spindle(data, sp):
    for i, t in enumerate(data.axis['time']):
        if sp['start_time'] > t[0] and sp['start_time'] < t[-1]:
            return i


def _find_sigma_ratio(one_sp, data):

    tb_0 = one_sp['peak_time'] - WINDOW / 2
    ib_0 = where(data.axis['time'][0] >= tb_0)[0][0]
    ib_1 = ib_0 + WINDOW * data.s_freq

    ia_0 = ib_0 - 1 * WINDOW * data.s_freq
    ia_1 = ib_0 - 0 * WINDOW * data.s_freq
    ic_0 = ib_0 + 1 * WINDOW * data.s_freq
    ic_1 = ib_0 + 2 * WINDOW * data.s_freq

    ia = data.data[0][0, ia_0:ia_1]
    ib = data.data[0][0, ib_0:ib_1]
    ic = data.data[0][0, ic_0:ic_1]

    sigma_ratio = sigma_power(ib, data.s_freq) / (sigma_power(ia, data.s_freq) +
                                                  sigma_power(ic, data.s_freq))

    return sigma_ratio


def sigma_power(x, s_freq):
    [f, Pxx] = periodogram(x, fs=s_freq, nfft=s_freq)
    SIGMA_FREQ = SPINDLE_OPTIONS['frequency']
    return sum(Pxx[(f >= SIGMA_FREQ[0]) & (f <= SIGMA_FREQ[1])])
