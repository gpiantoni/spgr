from functools import partial
from multiprocessing import Pool
from numpy import nanmean, NaN, isfinite, log, zeros, nansum, abs, array, max, nanmax, nanmin, r_

from .constants import (HEMI_SUBJ,
                        SPINDLE_OPTIONS,
                        SURF_PLOT_SIZE)
from .detect_spindles import get_spindles
from .plot_spindles import plot_lmer
from .spindle_source import get_chan_with_regions, get_regions_with_elec

from .log import with_log


from numpy import triu
from numpy.random import binomial, seed

NULL_PROBABILITY=.5
N_RND = 1000


@with_log
def Direction_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for reref in ('avg', 15):

        lg.info('### reref {}'.format(reref))

        regions = get_regions_with_elec(reref)

        x = zeros((len(regions), len(regions)))

        for subj in HEMI_SUBJ:

            spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
            chan = get_chan_with_regions(subj, reref=reref)
            chan_with_regions = dict(zip(spindles.chan_name,
                                         chan.return_attr('region',
                                                          spindles.chan_name)))

            for sp0 in spindles.spindle:  # this is the lead
                for sp1 in spindles.spindle:  # this is the follower

                    if ((sp1['start_time'] > sp0['start_time']) and
                        (sp1['start_time'] < sp0['end_time'])):

                        region0 = chan_with_regions[sp0['chan']]
                        region1 = chan_with_regions[sp1['chan']]
                        try:
                            i0 = regions.index(region0[7:])
                            i1 = regions.index(region1[7:])
                            x[i0, i1] += 1
                        except ValueError:  # if it's not in the regions of interest
                            pass

        d = sum(x, axis=1) / (sum(x, axis=0) + sum(x, axis=1)) * 100
        coef = dict(zip(regions, d))
        """
        _partial_null = partial(_compute_null_direction, x)

        with Pool() as p:
            n_d = p.map(_partial_null, range(N_RND))
        n_d = r_[n_d]


        lg.info('\nCorrected at Bonferroni 0.05')
        if pvalues[region] < p_threshold:
            lg.info('{:30} coef={:.3f}  p={:.4f}'.format(region, coef[region],
                                                         pvalues[region]))
        """
        v = plot_lmer(coef, limits=(40, 60), size_mm=SURF_PLOT_SIZE)
        png_name = 'direction_map_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))


def _compute_null_direction(x, seed_value):
    s_x = triu(x + x.T)

    seed(seed_value)
    x2 = zeros(x.shape)

    for i0 in range(x.shape[0]):
        for i1 in range(x.shape[1]):
            if i0 <= i1:
                b = binomial(s_x[i0, i1], NULL_PROBABILITY)
                x2[i0, i1] = b
                x2[i1, i0] = s_x[i0, i1] - b

    return sum(x2, axis=1) / (sum(x2, axis=0) + sum(x2, axis=1)) * 100
