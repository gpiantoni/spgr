from .constants import (HEMI_SUBJ,
                        HISTOGRAM_WIDTH,
                        SURF_PLOT_SIZE)
from .lmer_stats import add_to_dataframe, lmer
from .plot_spindles import plot_surf
from .plot_histogram import make_hist_overlap
from .spindle_source import get_morph_linear, get_chan_with_regions
from .stats_on_spindles import count_cooccur_per_chan

from .log import with_log


@with_log
def Cooccurrence_Histogram(lg, images_dir):

    lg.info('## Histogram of Co-occurrence of Spindles')

    for reref in ('avg', 15):
        for subj in HEMI_SUBJ:
            if reref == 'avg':
                nchan = 30
            else:
                nchan = 60
            v = make_hist_overlap(subj, reref=reref,
                                  width=HISTOGRAM_WIDTH, nchan=nchan)
            png_file = str(images_dir.joinpath('hist_{}_{}.png'.format(reref,
                                               subj)))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(reref, subj), png_file))


@with_log
def Cooccurrence_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for reref in ('avg', 15):

        lg.info('### reref {}'.format(reref))

        all_values = []
        dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

        for subj in HEMI_SUBJ:
            chan_val = count_cooccur_per_chan(subj, reref)

            morphed = get_morph_linear(subj, chan_val, reref=reref)
            all_values.append(morphed)

            chan = get_chan_with_regions(subj, reref)
            add_to_dataframe(dataframe, subj, chan_val, chan)

        cooccur_coef, _ = lmer(dataframe, lg)

        threshold = 0.01, None
        limits = 3, 6

        v = plot_surf(all_values, threshold=threshold, limits=limits,
                      size_mm=SURF_PLOT_SIZE)
        png_name = 'cooccurrence_map_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))
