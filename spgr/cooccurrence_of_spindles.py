from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        HISTOGRAM_WIDTH,
                        PLOT_COLOR,
                        SPINDLE_OPTIONS,
                        )
from .detect_spindles import get_spindles
from .lmer_stats import lmer
from .plot_spindles import plot_surf
from .plot_histogram import make_hist_overlap
from .read_data import get_chan_used_in_analysis
from .spindle_source import (get_chan_with_regions,
                             get_morph_linear,
                             )
from .stats_on_spindles import (create_spindle_cooccurrence_matrix,
                                create_spindle_groups,
                                mean_spindle_cooccurrence,
                                )

from .log import with_log


@with_log
def Cooccurrence_Histogram(lg, images_dir):

    lg.info('## Histogram of Co-occurrence of Spindles')

    for REREF in ('avg', 15):
        for subj in HEMI_SUBJ:
            if REREF == 'avg':
                nchan = 40
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

    all_values = []

    for REREF in ('avg', 15):
        for NORMALIZED_BY in ('source', 'target'):

            dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}
            for subj in HEMI_SUBJ:
                spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
                spindle_group = create_spindle_groups(spindles)

                chan = get_chan_used_in_analysis(subj, 'sleep',
                                                 chan_type=CHAN_TYPE,
                                                 reref=REREF, **DATA_OPTIONS)
                chan_prob = create_spindle_cooccurrence_matrix(chan,
                                                               spindle_group)
                chan_prob_per_chan = mean_spindle_cooccurrence(chan_prob,
                                                               normalized_by=NORMALIZED_BY)

                morphed = get_morph_linear(subj, chan_prob_per_chan,
                                           reref=REREF)
                all_values.append(morphed)

                chan = get_chan_with_regions(subj, REREF)
                add_to_dataframe(dataframe, subj, chan_prob_per_chan, chan)

            threshold = 0.01, None
            limits = 0, .15

            lmer(dataframe, lg)

            v = plot_surf(all_values, threshold=threshold, limits=limits)
            png_name = 'cooccurrence_map_{}_{}.png'.format(NORMALIZED_BY,
                                                           REREF)
            png_file = str(images_dir.joinpath(png_name))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(NORMALIZED_BY, REREF),
                    png_file))


def add_to_dataframe(df, subj, values, chan):

    for one_chan, one_value in zip(chan.chan, values):
        region = one_chan.attr['region']
        if region[:3] == 'ctx':
            df['subj'].append(subj)
            df['region'].append(region[7:])
            df['elec'].append(one_chan.label)
            df['value'].append(one_value)
