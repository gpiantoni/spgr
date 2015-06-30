from shutil import rmtree

from phypno.attr import Freesurfer

from spgr.constants import (IMAGES_PATH,
                            REC_PATH,
                            FS_FOLDER,
                            HEMI_SUBJ,
                            CHAN_TYPE,
                            DATA_OPTIONS,
                            )
from spgr.read_data import get_chan_used_in_analysis


from phypno.viz import Viz3
from numpy import ones
from spgr.spindle_source import get_morph_linear
from spgr.plot_spindles import plot_surf

from .log import with_log


SKIN_COLOR = (239, 208, 207, 240)


images_dir = IMAGES_PATH.joinpath('electrode_locations')
try:
    rmtree(str(images_dir))
except FileNotFoundError:
    pass
images_dir.mkdir(parents=True)

PLOT_COLOR = 'kw'


@with_log
def Electrode_Locations(lg):

    lg.info('## Locations for all the subjects')

    for subj, hemi in HEMI_SUBJ.items():
        chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, **DATA_OPTIONS)
        fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
        surf = getattr(fs.read_brain(), hemi)

        v = Viz3(color='kw')
        v.add_chan(chan, color=(200, 20, 20, 255))
        v.add_surf(surf, color=SKIN_COLOR)

        png_file = str(images_dir.joinpath(subj + '.png'))
        v.save(png_file)
        lg.info('![{}]({})'.format(subj, png_file))


    lg.info('## Coverage')

    morphed = []

    for subj in HEMI_SUBJ:

        # set channel values to 1
        chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, **DATA_OPTIONS)
        values = ones(chan.n_chan)

        morphed.append(get_morph_linear(subj, values, reref='avg'))

    v = plot_surf(morphed, limits=(0, 1), extra_smoothing=False)

    png_file = str(images_dir.joinpath('coverage_average.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('coverage', png_file))
