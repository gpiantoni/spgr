from numpy import nanmean, NaN, isfinite, log, zeros

from .constants import (COOCCUR_CHAN_LIMITS,
                        HEMI_SUBJ,
                        PARAMETERS,
                        SPINDLE_OPTIONS,
                        SURF_PLOT_SIZE)
from .detect_spindles import get_spindles
from .lmer_stats import add_to_dataframe, lmer
from .plot_spindles import plot_lmer
from .spindle_source import get_chan_with_regions

from .log import with_log


@with_log
def Direction_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for reref in ('avg', 15):

        limits = COOCCUR_CHAN_LIMITS[reref]

        lg.info('### reref {}'.format(reref))

        dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

        for subj in HEMI_SUBJ:
            chan_val = spindle_direct(subj, reref)
            chan = get_chan_with_regions(subj, reref)
            add_to_dataframe(dataframe, subj, chan_val, chan)

        lg.info('\nCorrected at FDR 0.05')
        coef, pvalues = lmer(dataframe, lg)

        v = plot_lmer(coef, pvalues=pvalues, limits=limits,
                      size_mm=SURF_PLOT_SIZE)
        png_name = 'cooccurrence_map_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))


def spindle_direct(subj, reref):

    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
    chan_names = list(spindles.chan_name)
    n_chan = len(spindles.chan_name)

    chan_names = list(spindles.chan_name)
    n_chan = len(spindles.chan_name)

    x = zeros((n_chan, n_chan))

    for sp0 in spindles.spindle:
        for sp1 in spindles.spindle:
            if (sp0['start_time'] > sp1['start_time']) and (sp0['start_time'] < sp1['end_time']):
                i0 = chan_names.index(sp0['chan'])
                i1 = chan_names.index(sp1['chan'])
                x[i0, i1] += 1

    c = log(x / x.T)
    c[~ isfinite(c)] = NaN

    return nanmean(c, axis=0)
