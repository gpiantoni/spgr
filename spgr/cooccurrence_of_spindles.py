from numpy import max, mean, min
from pickle import Pickler, Unpickler

from .constants import (ALL_REREF,
                        COOCCUR_CHAN_LIMITS,
                        GROUP_PATH,
                        HEMI_SUBJ,
                        PARAMETERS,
                        P_CORRECTION,
                        P_THRESHOLD,
                        SURF_PLOT_SIZE)
from .lmer_stats import add_to_dataframe, lmer
from .plot_spindles import plot_lmer
from .plot_histogram import make_hist_overlap
from .spindle_source import get_chan_with_regions
from .stats_on_spindles import count_cooccur_per_chan, get_cooccur_percent

from .log import with_log


@with_log
def Cooccurrence_Histogram(lg, images_dir):

    lg.info('## Histogram of Co-occurrence of Spindles')

    for reref in ALL_REREF:
        lg.info('### reref {}'.format(reref))

        all_p = []
        for subj in HEMI_SUBJ:

            v, p_mean = make_hist_overlap(subj, reref=reref)
            all_p.append(p_mean)
            lg.info('mean channels with co-occuring spindles for {}: {: 6.2f}'
                    ''.format(subj, p_mean))

            png_file = str(images_dir.joinpath('hist_{}_{}.png'.format(reref,
                                               subj)))
            v.save(png_file)
            lg.info('![{}]({})'.format('{} {}'.format(reref, subj), png_file))

        lg.info('Average number of channels with co-occurring spindles:'
                'mean {: 6.2f}, range {: 6.2f} - {: 6.2f}'
                ''.format(mean(all_p), min(all_p), max(all_p)))


@with_log
def Cooccurrence_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for reref in ALL_REREF:

        limits = COOCCUR_CHAN_LIMITS[reref]

        lg.info('### reref {}'.format(reref))

        dataframe_file = GROUP_PATH / ('dataframe_' + reref + '_' +
                                       PARAMETERS['summarize_cooccur'] + '.pkl')

        if dataframe_file.exists():
            with dataframe_file.open('rb') as f:
                dataframe = Unpickler(f).load()
        else:

            dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

            for subj in HEMI_SUBJ:
                chan_val = count_cooccur_per_chan(subj, reref,
                                                  PARAMETERS['summarize_cooccur'])
                chan = get_chan_with_regions(subj, reref)
                add_to_dataframe(dataframe, subj, chan_val, chan)

            with dataframe_file.open('wb') as f:
                Pickler(f).dump(dataframe)

        lg.info('\nCorrected at {} {}'.format(P_CORRECTION, P_THRESHOLD))
        coef, pvalues = lmer(dataframe, lg, adjust=P_CORRECTION,
                             pvalue=P_THRESHOLD)
        coef['insula_2'] = coef['middletemporal_2']  # insula has only 2 datapoints
        v = plot_lmer(coef, pvalues=pvalues, limits=limits,
                      size_mm=SURF_PLOT_SIZE)
        png_name = 'cooccurrence_map_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))


@with_log
def Cooccurrence_Percentile(lg, images_dir):

    lg.info('## Compare Top v Bottom Percentile of Co-occurrence of Spindles')

    for reref in ALL_REREF:
        lg.info('### reref {}'.format(reref))

        for subj in HEMI_SUBJ:
            get_cooccur_percent(subj, reref, lg)

