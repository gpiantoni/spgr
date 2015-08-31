from numpy import log, mean, where
from phypno.trans import Filter, Select
from phypno.viz import Viz1
from scipy.signal import periodogram

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        HIGHLIGHT_COLOR,
                        HIGHLIGHT_HEIGHT,
                        PLOT_COLOR,
                        RAW_LIMITS_Y,
                        SPINDLE_OPTIONS,
                        avg_regions,
                        avg_vert,
                        surf_avg,
                        )
from .detect_spindles import get_spindles
from .read_data import get_data
from .spindle_source import get_chan_with_regions

from .log import with_log

REREF = 'avg'
PAD = 1
WINDOW = 1
HIGHLIGHT_FILTER = 2, 30

BEST_SPINDLE = 'sigma_ratio'  # 'area_under_curve' or 'sigma_ratio'


@with_log
def Representative_Examples(lg, images_dir):
    """Show some good spindle examples."""
    lg.info('## Representative Spindles')

    # create dict with empty lists
    regions_with_png = {k: [] for k in avg_regions}

    for subj in HEMI_SUBJ:

        data = get_data(subj, 'sleep', CHAN_TYPE, reref=REREF, **DATA_OPTIONS)
        best_spindles = find_best_spindles(subj, data)

        filt = Filter(low_cut=HIGHLIGHT_FILTER[0],
                      high_cut=HIGHLIGHT_FILTER[1],
                      s_freq=data.s_freq)
        data = filt(data)

        for region, spindle in best_spindles.items():

            spindle_data = find_spindle_data(data, spindle)

            v = Viz1(color=PLOT_COLOR)
            v.add_data(spindle_data, limits_y=RAW_LIMITS_Y)
            v.add_graphoelement([spindle, ], color=HIGHLIGHT_COLOR,
                                height=HIGHLIGHT_HEIGHT)

            png_file = str(images_dir.joinpath('{}_{}.png'.format(region,
                                                                  subj)))
            v.save(png_file)

            regions_with_png[region].append(png_file)

    for region, png_files in regions_with_png.items():
        _, y, z = _find_region_xyz(region)
        lg.info('### {}'.format(region))
        lg.info('y: {}, z: {}'.format(y, z))
        for png_file in png_files:
            lg.info('![]({})'.format(png_file))


def find_best_spindles(subj, data):
    """For each subject and brain region, find the best spindle.

    Parameters
    ----------
    subj : str
        subj code
    data : instance of ChanTimeData
        complete recordings for a subject

    Returns
    -------
    dict
        dictionary, each brain region has the best spindle.
    """
    spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
    chan = get_chan_with_regions(subj, REREF)

    best_spindles = {}

    for one_sp in spindles:
        sp_region = chan(lambda x: x.label == one_sp['chan']).return_attr('region')[0]

        # aparc uses "ctx-Xh" and aparc.a2009s uses "ctx_Xh"
        if not sp_region.startswith('ctx'):
            continue

        sp_region = sp_region[len('ctx-Xh-'):]

        # there are different ways to define the best spindle.
        # the quality of a spindle will be stored in goodness (the higher the better)
        if BEST_SPINDLE == 'area_under_curve':
            one_sp['goodness'] = one_sp['area_under_curve']
        elif BEST_SPINDLE == 'sigma_ratio':
            try:
                one_sp['goodness'] = _find_sigma_ratio(one_sp, data)
            except (IndexError, ZeroDivisionError):  # how is this possible?
                one_sp['goodness'] = - 10000000

        if sp_region in best_spindles:
            if one_sp['goodness'] > best_spindles[sp_region]['goodness']:
                best_spindles[sp_region] = one_sp
        else:
            best_spindles[sp_region] = one_sp

    return best_spindles


def find_spindle_data(data, spindle):
    """Return the data around a spindle, for the channel with a spindle.

    Parameters
    ----------
    data : instance of ChanTimeData
        complete recordings for a subject
    spindle : dict
        parameters of one spindle

    Returns
    -------
    instance of ChanTimeData
        recordings for the channel with a spindle, time interval before and
        after the spindle
    """
    i_trial = _find_trial_with_spindle(data, spindle)

    center_time = (spindle['end_time'] + spindle['start_time']) / 2
    sel = Select(time=(center_time - PAD, center_time + PAD),
                 trial=(i_trial, ), chan=(spindle['chan'], ))
    return sel(data)


def _find_trial_with_spindle(data, sp):
    """Simple function to find the trial containing the spindle."""
    for i, t in enumerate(data.axis['time']):
        if sp['start_time'] > t[0] and sp['start_time'] < t[-1]:
            return i


def _find_sigma_ratio(one_sp, data):
    """For one spindle, get the sigma ratio (amount of spindle power).

    Parameters
    ----------
    one_sp : dict
        parameters of one spindle
    instance of ChanTimeData
        data with spindle data

    Returns
    -------
    float
        amount of spindle power, normalized by other frequencies / intervals
    """
    sp_t = where(data.axis['time'][0] >= one_sp['peak_time'])[0][0]

    ia_0 = sp_t - int(data.s_freq * WINDOW * 1.5)
    ia_1 = sp_t - int(data.s_freq * WINDOW * 0.5)
    ib_0 = sp_t - int(data.s_freq * WINDOW * 0.5)
    ib_1 = sp_t + int(data.s_freq * WINDOW * 0.5)
    ic_0 = sp_t + int(data.s_freq * WINDOW * 0.5)
    ic_1 = sp_t + int(data.s_freq * WINDOW * 1.5)

    ia = data.data[0][0, ia_0:ia_1]
    ib = data.data[0][0, ib_0:ib_1]
    ic = data.data[0][0, ic_0:ic_1]

    pa_0, pa_1, pa_2 = power_in_bands(ia, data.s_freq)
    pb_0, pb_1, pb_2 = power_in_bands(ib, data.s_freq)
    pc_0, pc_1, pc_2 = power_in_bands(ic, data.s_freq)

    sigma_ratio = (pb_1 / pb_2) / ((pa_1 / pa_2 + pc_1 / pc_2) +
                                   (pa_0 + pb_0 + pc_0))

    return sigma_ratio


def power_in_bands(x, s_freq):
    """Power below spindle band, in the spindle band, above spindle band."""
    [f, Pxx] = periodogram(x, fs=s_freq, nfft=s_freq)
    SIGMA_FREQ = SPINDLE_OPTIONS['frequency']

    Pxx = log(Pxx)

    below_sigma = sum(Pxx[(f >= 1) & (f <= SIGMA_FREQ[0])])
    sigma = sum(Pxx[(f >= SIGMA_FREQ[0]) & (f <= SIGMA_FREQ[1])])
    above_sigma = sum(Pxx[(f >= SIGMA_FREQ[1]) & (f <= 40)])
    return below_sigma, sigma, above_sigma


def _find_region_xyz(region):
    vert_in_region = avg_vert == avg_regions.index(region)
    xyz_vert_in_region = surf_avg.vert[vert_in_region, :]
    return mean(xyz_vert_in_region, axis=0)
