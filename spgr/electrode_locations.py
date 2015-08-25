from numpy import empty, NaN, ones
from phypno.attr import Freesurfer
from phypno.viz import Viz3

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        FS_FOLDER,
                        HEMI_SUBJ,
                        PLOT_COLOR,
                        REC_PATH,
                        SKIN_COLOR,
                        avg_vert,
                        avg_color,
                        avg_regions,
                        surf_avg,
                        )
from .lobes import freesurfer_color_code, LOBE_COLOR_LIMITS
from .spindle_source import get_morph_linear, get_regions_with_elec
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

        v = Viz3(color='kw')
        v.add_chan(chan, color=(200, 20, 20, 255))
        v.add_surf(surf, color=SKIN_COLOR)

        png_file = str(images_dir.joinpath(subj + '.png'))
        v.save(png_file)
        lg.info('![{}]({})'.format(subj, png_file))

    lg.info('## Coverage')
    """
    morphed = []
    for subj in HEMI_SUBJ:

        # set channel values to 1
        chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                         **DATA_OPTIONS)
        values = ones(chan.n_chan)

        morphed.append(get_morph_linear(subj, values, reref='avg'))

    v = plot_surf(morphed, limits=(0, 1), extra_smoothing=False)

    png_file = str(images_dir.joinpath('coverage_average.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('coverage', png_file))
    """

    lg.info('## Freesurfer Regions')
    regions_with_elec = get_regions_with_elec()

    colors = empty((surf_avg.vert.shape[0], 4))
    for i, vert in enumerate(avg_vert):
        if avg_regions[avg_vert[i]] in regions_with_elec:
            colors[i, :] = avg_color[vert, :4]
        else:
            colors[i, :] = NaN

    v = Viz3(color=PLOT_COLOR)
    v.add_surf(surf_avg, color=colors)
    png_file = str(images_dir.joinpath('fs_regions.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('regions', png_file))

    lg.info('## Lobes')
    freesurfer_code, color_code = freesurfer_color_code()

    color_avg = avg_vert.copy()
    for i0, i1 in zip(freesurfer_code, color_code):
        color_avg[avg_vert == i0] = i1

    v = Viz3(color=PLOT_COLOR)
    v.add_surf(surf_avg, values=color_avg, colormap='jet',
               limits_c=LOBE_COLOR_LIMITS)
    png_file = str(images_dir.joinpath('fs_lobes.png'))
    v.save(png_file)
    lg.info('![{}]({})'.format('lobes', png_file))
