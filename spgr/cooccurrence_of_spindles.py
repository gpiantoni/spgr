from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        HISTOGRAM_WIDTH,
                        PLOT_COLOR,
                        SPINDLE_OPTIONS,
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
                                )

from .log import with_log


@with_log
def Cooccurrence_Histogram(lg, images_dir):

    lg.info('## Histogram of Co-occurrence of Spindles')

    for REREF in ('avg', 15):
        for subj in HEMI_SUBJ:
            if REREF == 'avg':
                nchan = 30
            else:
                nchan = 60
            v = make_hist_overlap(subj, color=PLOT_COLOR, reref=REREF,
                                  width=HISTOGRAM_WIDTH, nchan=nchan)
            png_file = str(images_dir.joinpath('hist_{}_{}.png'.format(REREF,
                                               subj)))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(REREF, subj), png_file))


@with_log
def Cooccurrence_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for REREF in ('avg', 15):
        for NORMALIZATION in ('cooccur01', 'cooccur05', 'cooccur10',
                              'exclusive'):

            lg.info('### reref {}, normalization {}'.format(REREF,
                                                            NORMALIZATION))

            all_values = []
            dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

            for subj in HEMI_SUBJ:
                spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
                spindle_group = create_spindle_groups(spindles)

                chan = get_chan_used_in_analysis(subj, 'sleep',
                                                 chan_type=CHAN_TYPE,
                                                 reref=REREF, **DATA_OPTIONS)

                if NORMALIZATION.startswith('cooccur'):
                    chan_prob = cooccur_likelihood(chan, spindle_group,
                                                   NORMALIZATION[-2:])

                elif NORMALIZATION == 'exclusive':
                    chan_prob = ratio_spindles_with_chan(chan, spindle_group)

                morphed = get_morph_linear(subj, chan_prob, reref=REREF)
                all_values.append(morphed)

                chan = get_chan_with_regions(subj, REREF)
                add_to_dataframe(dataframe, subj, chan_prob, chan)

            lmer(dataframe, lg)

            if NORMALIZATION.startswith('cooccur'):
                threshold = 0.01, 1
                limits = 0, 1
            elif NORMALIZATION == 'exclusive':
                threshold = 0.01, None
                limits = 0, .15

            v = plot_surf(all_values, threshold=threshold, limits=limits)
            png_name = 'cooccurrence_map_{}_{}.png'.format(NORMALIZATION,
                                                           REREF)
            png_file = str(images_dir.joinpath(png_name))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(NORMALIZATION, REREF),
                    png_file))
