from logging import getLogger

from phypno import Data
from phypno.attr import Freesurfer
from phypno.source import Linear, Morph

from .constants import (REC_PATH,
                        FS_FOLDER,
                        DEFAULT_HEMI,
                        HEMI_SUBJ,
                        CHAN_TYPE,
                        DATA_OPTIONS,
                        SMOOTHING_STD,
                        SMOOTHING_THRESHOLD,
                        )
from .read_data import get_chan_used_in_analysis


lg = getLogger(__name__)

SMOOTHING = {'std': None,
             'threshold': None
             }
LINEAR = dict.fromkeys(list(HEMI_SUBJ))


def get_morph_linear(subj, values, reref):
    lg.info('Projecting values for {}'.format(subj))

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
