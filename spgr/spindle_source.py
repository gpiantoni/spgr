from logging import getLogger

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
                        SMOOTHING_STD,
                        SMOOTHING_THRESHOLD,
                        PARAMETERS,
                        )
from .read_data import get_chan_used_in_analysis


lg = getLogger(__name__)

SMOOTHING = {'std': None,
             'threshold': None
             }
LINEAR = dict.fromkeys(list(HEMI_SUBJ))
REGIONS = dict.fromkeys(list(HEMI_SUBJ))


def get_morph_linear(subj, values, reref):
    lg.debug('Projecting values for {}'.format(subj))

    chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, reref=reref,
                                     **DATA_OPTIONS)

    data = Data(values, chan=chan.return_label())
    morphed_data = _reflect_to_avg(subj, data, chan)

    return morphed_data


def _reflect_to_avg(subj, data, chan):
    """Simplest and inaccurate way to project onto the hemisphere of interest.
    We just reflect the channels from the wrong side to the side of interest.
    """
    if HEMI_SUBJ[subj] != DEFAULT_HEMI:
        for one_chan in chan.chan:
            one_chan.xyz *= (-1, 1, 1)

    fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
    brain = fs.read_brain()
    surf = getattr(brain, DEFAULT_HEMI)

    same_smoothing = (SMOOTHING['std'] == SMOOTHING_STD and
                      SMOOTHING['threshold'] == SMOOTHING_THRESHOLD)
    # check if 1. was precomputed, 2. same n of channels, 3. same parameters
    if LINEAR[subj] and len(LINEAR[subj].chan) == chan.n_chan and same_smoothing:
        l = LINEAR[subj]

    else:
        l = Linear(surf, chan, std=SMOOTHING_STD,
                   threshold=SMOOTHING_THRESHOLD)

        # we assign to the values of SMOOTHING, not to SMOOTHING
        SMOOTHING['std'] = SMOOTHING_STD
        SMOOTHING['threshold'] = SMOOTHING_THRESHOLD
        LINEAR[subj] = l

    m = Morph(surf)
    morphed_data = m(l(data))

    return morphed_data


def get_chan_with_regions(subj, reref):

    orig_chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                          reref=reref, **DATA_OPTIONS)

    if REGIONS[subj] and REGIONS[subj].n_chan == orig_chan.n_chan:
        chan = REGIONS[subj]
    else:
        fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
        chan = assign_region_to_channels(orig_chan, fs,
                                         parc_type=PARAMETERS['PARC_TYPE'],
                                         exclude_regions=('Unknown', ))
        REGIONS[subj] = chan

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
        chan = get_chan_with_regions(subj, 'avg')
        all_regions.append(chan.return_attr('region'))

    region_names = set(x[7:] for x in all_regions if x[:4] == 'ctx_')
    lg.debug('Number of regions with at least one elec: {}'
             ''.format(len(region_names)))
    return list(region_names)
