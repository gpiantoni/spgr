from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        HEMI_SUBJ,
                        SPINDLE_OPTIONS,
                        SURF_PLOT_SIZE,
                        )
from .detect_spindles import get_spindles
from .lmer_stats import add_to_dataframe, lmer
from .plot_spindles import plot_surf
from .read_data import get_data
from .spindle_source import get_morph_linear, get_chan_with_regions

from .log import with_log


@with_log
def Single_Channel_Statistics(lg, images_dir):

    lg.info('## Spindle Properties: Descriptive')

    lg.info('### Reference to Average')
    plot_average_values('avg', lg, images_dir)

    lg.info('### Bipolar Reference')
    plot_average_values(15, lg, images_dir)


def get_spindle_param(subj, param, ref):
    """Param: 'density', 'duration', 'peak_freq', 'peak_val'"""
    spindles = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref,
                            **SPINDLE_OPTIONS)

    if param == 'density':
        data = get_data(subj, 'sleep', chan_type=CHAN_TYPE, reref=ref,
                        **DATA_OPTIONS)

        dat_count = spindles.to_data('count')
        n_min = (data.number_of('trial') *
                 (data.axis['time'][0][-1] - data.axis['time'][0][0])) / 60
        values = dat_count.data[0] / n_min

    elif param == 'duration':
        for one_spindle in spindles.spindle:
            one_spindle['duration'] = (one_spindle['end_time'] -
                                       one_spindle['start_time'])
        data = spindles.to_data('duration')
        values = data.data[0]

    else:
        data = spindles.to_data(param)
        values = data.data[0]

    return values


def plot_average_values(REREF, lg, images_dir):
    for param in ('density', 'duration', 'peak_val', 'peak_freq'):

        if param == 'density':
            threshold = 0.01, None
            limits = 0, 2
        elif param == 'duration':
            threshold = 0.5, 2
            limits = 1, 1.5
        elif param == 'peak_val':
            threshold = 0.001, None
            limits = 0, 70
        elif param == 'peak_freq':
            threshold = SPINDLE_OPTIONS['frequency']
            limits = (SPINDLE_OPTIONS['frequency'][0] + 2,
                      SPINDLE_OPTIONS['frequency'][1] - 2)

        morphed = []

        dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

        for subj in HEMI_SUBJ:
            values = get_spindle_param(subj, param, REREF)
            morphed.append(get_morph_linear(subj, values, REREF))

            chan = get_chan_with_regions(subj, REREF)
            add_to_dataframe(dataframe, subj, values, chan)

        lmer(dataframe, lg)
        v = plot_surf(morphed, threshold=threshold, limits=limits,
                      size_mm=SURF_PLOT_SIZE)
        png_file = str(images_dir.joinpath('{}_{}.png'.format(param, REREF)))
        v.save(png_file)
        lg.info('![{}]({})'.format('{} {}'.format(REREF, param), png_file))
