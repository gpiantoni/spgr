from collections import Counter
from logging import getLogger
from pickle import load, dump
from re import split

from numpy import array, max, mean, min

from phypno import Data
from phypno.attr import Freesurfer
from phypno.attr.chan import assign_region_to_channels
from phypno.source import Linear, Morph

from .constants import (APARC_FOLDER,
                        PROJ_FOLDER,
                        DATA_PATH,
                        REC_PATH,
                        FS_FOLDER,
                        DEFAULT_HEMI,
                        HEMI_SUBJ,
                        CHAN_TYPE,
                        DATA_OPTIONS,
                        MORPH_SMOOTHING,
                        SMOOTHING_STD,
                        SMOOTHING_THRESHOLD,
                        PARAMETERS,
                        )
from .read_data import get_chan_used_in_analysis


lg = getLogger(__name__)


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
    subj_dir = DATA_PATH / subj / PROJ_FOLDER
    if not subj_dir.exists():
        subj_dir.mkdir()
    linear_file = subj_dir / linear_filename
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


def rejected_chan(lg, all_subj, good_chan, all_chan):

    n_good_chan = array([x.n_chan for x in good_chan])
    lg.info('channels used in analysis: mean %.2f (range %d-%d)',
            mean(n_good_chan), min(n_good_chan), max(n_good_chan))

    gr_chan = []
    for chan in all_chan:
        gr_chan.append(chan(lambda x: 'gr' in x.label.lower()))
    n_gr_chan = array([x.n_chan for x in gr_chan])
    lg.info('channels with a grid: mean %.2f (range %d-%d)',
            mean(n_gr_chan), min(n_gr_chan), max(n_gr_chan))

    bad_chan = []
    for chan0, chan1 in zip(good_chan, gr_chan):
        good_labels = chan0.return_label()
        bad_chan.append(chan1(lambda x: x.label not in good_labels))
    n_bad_chan = array([x.n_chan for x in bad_chan])
    lg.info('bad channels with a grid: mean %.2f (range %d-%d)',
            mean(n_bad_chan), min(n_bad_chan), max(n_bad_chan))

    perc_bad_chan = n_bad_chan / n_gr_chan * 100
    lg.info('percent bad channels with a grid: mean %.2f (range %.2f-%.2f)',
            mean(perc_bad_chan), min(perc_bad_chan), max(perc_bad_chan))

    all_regions = []
    for subj, chan in zip(all_subj, bad_chan):
        chan = _assign_labels(subj, chan, PARAMETERS['PARC_TYPE'])
        regions = chan.return_attr('region')
        for region in regions:
            all_regions.append(split('[-_]', region)[2])

    for region, n_elec in Counter(all_regions).most_common():
        lg.info(region + ': ' + str(n_elec))


def get_chan_with_regions(subj, reref, parc_type=None):

    if parc_type is None:
        parc_type = PARAMETERS['PARC_TYPE']

    orig_chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE,
                                          reref=reref, **DATA_OPTIONS)

    region_filename = ('{}_chan{:03d}_{}.pkl'
                       ''.format(parc_type,
                                 orig_chan.n_chan, subj))
    subj_dir = DATA_PATH / subj / APARC_FOLDER
    if not subj_dir.exists():
        subj_dir.mkdir()

    region_file = subj_dir / region_filename
    if region_file.exists():
        with open(str(region_file), 'rb') as f:
            chan = load(f)
    else:
        chan = _assign_labels(subj, orig_chan, parc_type)
        with open(str(region_file), 'wb') as f:
            dump(chan, f)

    return chan


def _assign_labels(subj, chan, parc_type):
    """Subfunction to assign region labels to chan

    Parameters
    ----------
    subj : str
        subject code
    chan : instance of Channels
        the channels to compute the labels of
    parc_type : str
        the type of parcellation

    Returns
    -------
    instance of Channels
        channels with the regions
    """
    if parc_type.startswith('aparc.laus'):
        fs_lut_name = 'aparc.annot.' + parc_type.split('.')[-1] + '.ctab'
        fs_lut = str(REC_PATH / subj / FS_FOLDER / 'label' / fs_lut_name)
    else:
        fs_lut = None
    fs = Freesurfer(str(REC_PATH / subj / FS_FOLDER), fs_lut=fs_lut)
    chan = assign_region_to_channels(chan, fs, parc_type=parc_type,
                                     exclude_regions=('Unknown', ))
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

    now it sorts alphabetically but there should be a better way to organize
    it.
    """
    all_regions = []
    for subj in HEMI_SUBJ:
        chan = get_chan_with_regions(subj, reref)
        all_regions.extend(chan.return_attr('region'))

    region_names = set(x[7:] for x in all_regions if x[:3] == 'ctx')
    return sorted(list(region_names))
