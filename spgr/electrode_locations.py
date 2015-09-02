from numpy import ones
from phypno.attr import Freesurfer
from phypno.viz import Viz3

from .constants import (CHAN_COLOR,
                        CHAN_TYPE,
                        DATA_OPTIONS,
                        DPI,
                        FS_FOLDER,
                        HEMI_SUBJ,
                        REC_PATH,
                        SKIN_COLOR,
                        SURF_PLOT_SIZE,
                        SINGLE_SUBJ_SURF,
                        AVERAGE_BW_SURF,
                        surf_avg)
from .spindle_source import get_morph_linear
from .plot_spindles import plot_surf
from .read_data import get_chan_used_in_analysis

from .log import with_log


@with_log
def Electrode_Locations(lg, images_dir):

    lg.info('## Locations for all the subjects')
    for subj, hemi in HEMI_SUBJ.items():
        chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                         **DATA_OPTIONS)
        fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
        surf = getattr(fs.read_brain(), hemi)

        v = Viz3(show=False, dpi=DPI, size_mm=SINGLE_SUBJ_SURF)
        v.add_chan(chan, color=CHAN_COLOR)
        v.add_surf(surf, color=SKIN_COLOR)
        v._plt.camera.scale_factor = 150

        png_file = str(images_dir.joinpath(subj + '.png'))
        v.save(png_file)
        lg.info('![{}]({})'.format(subj, png_file))

    lg.info('## Coverage')
    morphed = []
    for subj in HEMI_SUBJ:

        # set channel values to 1
        chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                         **DATA_OPTIONS)
        values = ones(chan.n_chan)

        morphed.append(get_morph_linear(subj, values, reref='avg'))

    v = plot_surf(morphed, limits=(0, 1), extra_smoothing=False,
                  size_mm=SURF_PLOT_SIZE)

    png_file = str(images_dir.joinpath('coverage_average.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('coverage', png_file))

    lg.info('## Average surface')
    v = Viz3(size_mm=AVERAGE_BW_SURF, dpi=DPI)
    v.add_surf(surf_avg, color=(1, 1, 1))
    png_file = str(images_dir.joinpath('fs_avg.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('surface average', png_file))
