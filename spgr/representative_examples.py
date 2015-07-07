from phypno.trans import Select
from phypno.viz import Viz1

from spgr.constants import (CHAN_TYPE,
                            DATA_OPTIONS,
                            HEMI_SUBJ,
                            PLOT_COLOR,
                            SPINDLE_OPTIONS,
                            )
from spgr.detect_spindles import get_spindles
from spgr.read_data import get_data
from spgr.spindle_source import get_chan_with_regions

from .log import with_log


REREF = 'avg'
PAD = 2
YLIM = 150


@with_log
def Representative_Examples(lg, images_dir):
    """Show some good spindle examples.
    A better way of doing this is by doing fft at spindle, before and after and
    take the spindle with the best ratio.
    """
    lg.info('## Representative Spindles')

    lg.info('### Spindles in each brain region')

    all_pos = []
    all_title = []
    all_png_file = []

    for subj in HEMI_SUBJ:
        best_spindles = find_best_spindles(subj)

        data = get_data(subj, 'sleep', CHAN_TYPE, reref='avg', **DATA_OPTIONS)
        chan = data.attr['chan']

        for region, spindle in best_spindles.items():

            i_trial = find_trial_with_spindle(data, spindle)

            sel = Select(time=(spindle['start_time'] - PAD,
                               spindle['end_time'] + PAD),
                         trial=(i_trial, ), chan=(spindle['chan'], ))
            spindle_data = sel(data)

            v = Viz1(color=PLOT_COLOR)
            v.add_data(spindle_data, limits_y=(-YLIM, YLIM))
            v.add_graphoelement([spindle, ])

            png_file = str(images_dir.joinpath('{}_{}.png'.format(region, subj)))
            v.save(png_file)

            y_pos = chan(lambda x: x.label == spindle['chan']).return_xyz()[0, 1]
            all_pos.append(y_pos)
            all_title.append('### Representative for {} ({})'.format(region,
                                                                     subj))
            all_png_file.append(png_file)

    # sort by anterior-posterior axis
    for _, title, png_file in sorted(zip(all_pos, all_title, all_png_file)):
        lg.info(title)
        lg.info('![]({})'.format(png_file))


def find_best_spindles(subj):

    spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
    chan = get_chan_with_regions(subj, REREF)

    best_spindles = {}

    for one_sp in spindles:
        sp_region = chan(lambda x: x.label == one_sp['chan']).return_attr('region')[0]

        if not sp_region.startswith('ctx-'):
            continue

        sp_region = sp_region[len('ctx-Xh-'):]

        if sp_region in best_spindles:
            if one_sp['area_under_curve'] > best_spindles[sp_region]['area_under_curve']:
                best_spindles[sp_region] = one_sp
        else:
            best_spindles[sp_region] = one_sp

    return best_spindles


def find_trial_with_spindle(data, sp):
    for i, t in enumerate(data.axis['time']):
        if sp['start_time'] > t[0] and sp['start_time'] < t[-1]:
            return i
