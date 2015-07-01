from shutil import rmtree

from spgr.constants import (HEMI_SUBJ,
                            CHAN_TYPE,
                            DATA_OPTIONS,
                            SPINDLE_OPTIONS,
                            )
from spgr.read_data import get_data
from spgr.detect_spindles import get_spindles
from spgr.spindle_source import get_morph_linear, get_chan_with_regions
from spgr.plot_spindles import plot_surf
from spgr.lmer_stats import lmer
from spgr.constants import IMAGES_PATH

# from multiprocessing import set_start_method
# set_start_method('forkserver')


from .log import with_log

images_dir = IMAGES_PATH.joinpath('property_descriptive')
try:
    rmtree(str(images_dir))
except FileNotFoundError:
    pass
images_dir.mkdir(parents=True)

PLOT_COLOR = 'kw'


@with_log
def Single_Channel_Statistics(lg):

    lg.info('## Spindle Properties: Descriptive')

    lg.info('### Reference to Average')
    plot_average_values('avg', lg)

    lg.info('### Bipolar Reference')
    plot_average_values(15, lg)


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


def add_to_dataframe(df, subj, values, chan):

    for one_chan, one_value in zip(chan.chan, values):
        region = one_chan.attr['region']
        if region[:3] == 'ctx':
            df['subj'].append(subj)
            df['region'].append(region[7:])
            df['elec'].append(one_chan.label)
            df['value'].append(one_value)


def plot_average_values(REF, lg):
    lg.propagate = False
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
            limits = 11, 13

        morphed = []

        dataframe = {'subj': [], 'region': [], 'elec': [], 'value': []}

        for subj in HEMI_SUBJ:
            values = get_spindle_param(subj, param, REF)
            morphed.append(get_morph_linear(subj, values, REF))

            chan = get_chan_with_regions(subj, REF)
            add_to_dataframe(dataframe, subj, values, chan)

        lg.propagate = True
        lmer(dataframe)
        lg.propagate = False
        v = plot_surf(morphed, threshold=threshold, limits=limits)
        png_file = str(images_dir.joinpath('{}_{}.png'.format(param, REF)))
        v.save(png_file)
        lg.info('![{}]({})'.format('{} {}'.format(REF, param), png_file))

    lg.propagate = True
