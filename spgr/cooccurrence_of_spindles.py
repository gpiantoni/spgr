from numpy import array, sqrt, polyval, where
from scipy.stats import linregress
from vispy.io import write_png
from vispy.plot import Fig
from vispy.scene import Text, LinePlot

from copy import deepcopy
from numpy import array, mean, std
from numpy.random import uniform, seed
from multiprocessing import Pool
from functools import partial

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        HISTOGRAM_WIDTH,
                        SPINDLE_OPTIONS,
                        SURF_PLOT_SIZE
                        )
from .detect_spindles import get_spindles
from .lmer_stats import add_to_dataframe, lmer
from .plot_spindles import plot_surf
from .plot_histogram import make_hist_overlap
from .read_data import get_chan_used_in_analysis, get_data
from .spindle_source import (get_chan_with_regions,
                             get_morph_linear,
                             )
from .stats_on_spindles import (cooccur_likelihood,
                                create_spindle_groups,
                                ratio_spindles_with_chan,
                                count_sp_at_any_time,
                                )

from .log import with_log


@with_log
def Cooccurrence_Histogram(lg, images_dir):

    lg.info('## Histogram of Co-occurrence of Spindles')

    for REREF in ('avg', 15)[:1]:
        for subj in HEMI_SUBJ:
            if REREF == 'avg':
                nchan = 30
            else:
                nchan = 60
            v = make_hist_overlap(subj, reref=REREF,
                                  width=HISTOGRAM_WIDTH, nchan=nchan)
            png_file = str(images_dir.joinpath('hist_{}_{}.png'.format(REREF,
                                               subj)))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(REREF, subj), png_file))


@with_log
def Cooccurrence_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for REREF in ('avg', 15):

        lg.info('### reref {}'.format(REREF))

        all_values = []
        dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

        for subj in HEMI_SUBJ:
            spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)

            chan = get_chan_used_in_analysis(subj, 'sleep', chan_type=CHAN_TYPE, reref=REREF, **DATA_OPTIONS)
            chan = get_chan_with_regions(subj, REREF)

            data = get_data(subj, period_name='sleep', chan_type=CHAN_TYPE, reref=REREF, **DATA_OPTIONS)
            t_range = concatenate(data.axis['time'][:])

            p = count_sp_at_any_time(spindles, t_range)
            t_range_with_sp = t_range[p >= 1]
            p_with_sp = p[p >= 1]
            chan_prob = zeros(chan.n_chan)
            chan = get_chan_with_regions(subj, REREF)

            for i, one_chan in enumerate(chan.chan):
                t_at_sp = zeros(t_range_with_sp.shape, dtype=bool)
                for one_sp in spindles.spindle:
                    if one_sp['chan'] == one_chan.label:
                        t_at_sp = t_at_sp | ((t_range_with_sp >= one_sp['start_time']) & (t_range_with_sp < one_sp['end_time']))

                chan_prob[i] = mean(p_with_sp[t_at_sp]) # / mean(p_with_sp[~t_at_sp])

            morphed = get_morph_linear(subj, chan_prob, reref=REREF)
            all_values.append(morphed)

            add_to_dataframe(dataframe, subj, chan_prob, chan)

        cooccur_coef, _ = lmer(dataframe, lg)

        threshold = 0.01, None
        limits = 3, 6

        v = plot_surf(all_values, threshold=threshold, limits=limits,
                      size_mm=SURF_PLOT_SIZE)
        png_name = 'cooccurrence_map_{}.png'.format(REREF)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(REREF),
                png_file))


def _make_fake_chan_prob(i, start_min, start_max, chan, spindles):
    seed(i)
    fake_spindles = deepcopy(spindles)
    for fake_sp in fake_spindles.spindle:
        dur = fake_sp['end_time'] - fake_sp['start_time']
        fake_sp['start_time'] = uniform(start_min, start_max)
        fake_sp['end_time'] = fake_sp['start_time'] + dur

    fake_spindle_group = create_spindle_groups(fake_spindles)
    fake_chan_prob = ratio_spindles_with_chan(chan, fake_spindle_group)
    return fake_chan_prob


def shuffle_ratio(chan, spindles):
    spindle_group = create_spindle_groups(spindles)
    ref_chan_prob = ratio_spindles_with_chan(chan, spindle_group)

    start_min = min(x['start_time'] for x in spindles.spindle)
    start_max = max(x['start_time'] for x in spindles.spindle)

    _make_fake_partial = partial(_make_fake_chan_prob, start_min=start_min,
                                 start_max=start_max, chan=chan,
                                 spindles=spindles)
    N_RPT = 2000
    with Pool() as p:
        v = list(p.map(_make_fake_partial, range(N_RPT)))

    b = array(v)
    return (ref_chan_prob - mean(b, axis=0)) / std(b, axis=0)


def _cooccur_v_density(dens_coef, cooccur_coef, lg):
    """Compare spindle per minute vs cooccurrence to see if they are linked and
    if there are some exceptions, like some regions that have higher spindle
    cooccurrence than explained only by spindle count.
    """
    regions = array(sorted(dens_coef))
    dens = array([dens_coef[one_region] for one_region in regions])
    cooccur = array([cooccur_coef[one_region] for one_region in regions])

    lreg = linregress(dens, cooccur)
    lg.info('Correlation spindle density v cooccurrence r({}) = {:.3f}, '
            'p = {:.3f}'.format(len(regions) - 2, lreg.rvalue, lreg.pvalue))
    lg.info('Best-fit line: {: 5.3f}y + {: 5.3f}'.format(lreg.slope,
                                                         lreg.intercept))
    estimates = polyval((lreg.slope, lreg.intercept), dens)

    dist_from_line = cooccur - estimates
    deviance = sqrt(sum(dist_from_line ** 2) / (len(dist_from_line) - 2))
    lg.info('deviance: {: 5.3f}'.format(deviance))

    def _find_region_xy(region_name):
        idx = where(regions == region_name)[0][0]
        return array((dens[idx], cooccur[idx]))

    threshold = deviance
    offset = array([0.03, 0])

    fig = Fig(show=True)
    f = fig[0, 0]
    f._configure_2d()
    f.view.add(LinePlot((dens, cooccur), width=0, symbol='*'))
    f.view.add(LinePlot((dens, estimates), width=3))

    for one_region in regions[dist_from_line > threshold]:
        pos = _find_region_xy(one_region) - offset
        tv = Text(text=one_region, pos=pos, anchor_x='right')
        f.view.add(tv)

    for one_region in regions[dist_from_line < -threshold]:
        pos = _find_region_xy(one_region) + offset
        tv = Text(text=one_region, pos=pos, anchor_x='left')
        f.view.add(tv)

    f.view.camera.set_range(x=(0.3, 2), y=(0.05, 0.15))

    f.yaxis.axis._text.font_size = 10
    f.xaxis.axis._text.font_size = 10

    f.xlabel.text = 'Estimated Spindles per Minute'
    f.ylabel.text = 'Estimated Co-occurrence Probability'

    return fig
