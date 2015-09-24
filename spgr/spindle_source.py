from logging import getLogger
from pickle import load, dump

from phypno import Data
from phypno.attr import Freesurfer
from phypno.attr.chan import assign_region_to_channels
from phypno.source import Linear, Morph

from .constants import (REC_PATH,
                        FS_FOLDER,
                        DEFAULT_HEMI,
                        HEMI_SUBJ,
                        CHAN_TYPE,
                        DATA_OPTIONS,
                        GROUP_PATH,
                        MORPH_SMOOTHING,
                        SMOOTHING_STD,
                        SMOOTHING_THRESHOLD,
                        PARAMETERS,
                        )
from .read_data import get_chan_used_in_analysis


lg = getLogger(__name__)

STORED_PATH = GROUP_PATH.joinpath('saved_regions')
if not STORED_PATH.exists():
    STORED_PATH.mkdir()


def get_morph_linear(subj, values, reref, to_surf='fsaverage'):
    lg.debug('Projecting values for {}'.format(subj))

    chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, reref=reref,
                                     **DATA_OPTIONS)

    data = Data(values, chan=chan.return_label())
    morphed_data = _reflect_to_avg(subj, data, chan, to_surf)

    return morphed_data


def _reflect_to_avg(subj, data, chan, to_surf):
    """Simplest and inaccurate way to project onto the hemisphere of interest.
    We just reflect the channels from the wrong side to the side of interest.
    """
    if HEMI_SUBJ[subj] != DEFAULT_HEMI:
        for one_chan in chan.chan:
            one_chan.xyz *= (-1, 1, 1)

    fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
    brain = fs.read_brain()
    surf = getattr(brain, DEFAULT_HEMI)

    linear_filename = ('linear_chan{:03d}_std{:03d}_thr{:03d}_{}.pkl'
                       ''.format(chan.n_chan, SMOOTHING_STD,
                                 SMOOTHING_THRESHOLD, subj))
    linear_file = STORED_PATH.joinpath(linear_filename)
    if linear_file.exists():
        with open(str(linear_file), 'rb') as f:
            l = load(f)
    else:
        l = Linear(surf, chan, std=SMOOTHING_STD,
                   threshold=SMOOTHING_THRESHOLD)
        with open(str(linear_file), 'wb') as f:
            dump(l, f)

    m = Morph(surf, to_surf=to_surf, smooth=MORPH_SMOOTHING)
    morphed_data = m(l(data))

    return morphed_data


def get_chan_with_regions(subj, reref, parc_type=None):

    if parc_type is None:
        parc_type = PARAMETERS['PARC_TYPE']

    orig_chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                          reref=reref, **DATA_OPTIONS)

    region_filename = ('{}_chan{:03d}_{}.pkl'
                       ''.format(parc_type,
                                 orig_chan.n_chan, subj))
    region_file = STORED_PATH.joinpath(region_filename)
    if region_file.exists():
        with open(str(region_file), 'rb') as f:
            chan = load(f)
    else:
        fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
        chan = assign_region_to_channels(orig_chan, fs, parc_type=parc_type,
                                         exclude_regions=('Unknown', ))
        with open(str(region_file), 'wb') as f:
            dump(chan, f)

    return chan


def get_regions_with_elec(reref='avg'):
    """Return the list of channels with at least one electrode.

    It loops over every subject, assigns the region to each electrode, and then
    takes the unique region names.

    Notes
    -----
    Region names start with "ctx_?h" for aparc.a2009s at least. So we remove
    that part of the name. By using "if" we also get rid of regions such as
    white matter or other weird regions.
    """
    all_regions = []
    for subj in HEMI_SUBJ:
        chan = get_chan_with_regions(subj, reref)
        all_regions.extend(chan.return_attr('region'))

    region_names = set(x[7:] for x in all_regions if x[:4] == 'ctx_')
    return list(region_names)
