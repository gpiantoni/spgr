from glob import glob
from hashlib import md5
from json import load as json_load
from logging import getLogger, DEBUG
from os import listdir, makedirs, sep
from os.path import isdir, join, exists, expanduser, splitext, basename
from pickle import dump, load

from phypno import Dataset
from phypno.attr import Annotations, Channels
from phypno.trans import Filter, Montage, Resample
from phypno.viz import Viz1


lg = getLogger(__file__)
lg.setLevel(DEBUG)

PROJECT = 'spgr'
HOME = expanduser('~')

REC_DIR = join(HOME, 'recordings')
XLTEK_PATH = 'eeg/raw/xltek'
SCORE_PATH = 'doc/scores'
ELEC_PATH = 'doc/elec'
FS_PATH = 'mri/proc/freesurfer'

DATA_DIR = join(HOME, 'projects', PROJECT, 'subjects')
GROUP_DIR = join(HOME, 'projects', PROJECT, 'group')

REC_FOLDER = 'rec'
HP_FILTER = 1
LP_FILTER = 40
RESAMPLE_FREQ = 256
REREF = ''
REC_NAME = ('{subj}_{period}_{chan_types}_hp{hp:03d}_lp{lp:03d}{reref}'
            '_rs{resample:03d}.pkl')


def select_scores(stages, duration, all_subj, choose='max'):
    """Select subjects/nights depending on some criteria.

    Parameters
    ----------
    stages : tuple of str
        list of stages to analyze
    duration : int
        minimal duration to include an epoch, in seconds
    all_subj : list of str
        list of all the subjects
    choose : str, optional
        criterion to choose best night ('max': period with longest duration in
        stage, 'latest': the last possible period)

    Returns
    -------
    dict
        where key is the subject code and the value is the path to scoring file
        with the longest duration in the stage(s) of interest.

    """
    scores = {}
    for subj in all_subj:
        try:
            scores[subj] = _select_scores_per_subj(stages, duration, subj,
                                                   choose)
            lg.info(subj + ' has ' + scores[subj])
        except IndexError:
            lg.debug(subj + ': no scored recordings')
        except ValueError:
            lg.info(subj + ': duration not sufficient')

    return scores


def save_data(subj, score_file, period_name, stages, chan_type=(),
              hp_filter=HP_FILTER, lp_filter=LP_FILTER,
              resample_freq=RESAMPLE_FREQ, reref=REREF, to_plot=False):
    """Save recordings for one subject, based on some parameters

    Parameters
    ----------
    subj : str
        subject code
    score_file : path to file
        file with sleep scoring and annotations
    period_name : str
        name of the period, used to idenfity it (not important)
    stages : tuple of str
        list of stages to analyze
    hp_filter : float, optional
        high-pass filter
    lp_filter : float, optional
        low-pass filter
    resample_freq : int, optional
        frequency to resample to
    reref : None or str, optional
        if data should be rereferenced to the average ('avg')
    chan_type : tuple
        tuple of str to select channel groups, among 'depth', 'grid', 'scalp'
    to_plot : bool
        if you want a plot of an epoch in the middle of the recordings

    Notes
    -----
    It saves the data to disk. The name of the file contains the parameters.
    """
    lg.info(subj + ' ' + score_file)
    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    if not isdir(subj_dir):
        makedirs(subj_dir)

    xltek_file = join(REC_DIR, subj, XLTEK_PATH, score_file[:-11])
    d = Dataset(xltek_file)
    lg.info('Sampling Frequency {} Hz '.format(d.header['s_freq']))

    score_dir = join(REC_DIR, subj, SCORE_PATH)
    score = Annotations(join(score_dir, score_file))

    chan_file = join(REC_DIR, subj, ELEC_PATH,
                     basename(score_file).replace('_scores.xml',
                                                  '_channels.json'))

    selected_chan = _select_channels(chan_file, chan_type)
    lg.info('N Channels {} '.format(len(selected_chan)))
    lg.debug('Channels: ' + ', '.join(selected_chan))

    start_time = [x['start'] for x in score.epochs if x['stage'] in stages]
    end_time = [x['end'] for x in score.epochs if x['stage'] in stages]
    duration = sum([x1 - x0 for x0, x1 in zip(start_time, end_time)])
    lg.info('Duration: {0: 3.1f} min'.format(duration / 60))

    data = d.read_data(begtime=start_time, endtime=end_time,
                       chan=selected_chan)

    # first re-reference, and then filter
    if reref == 'avg':
        to_avg = Montage(ref_to_avg=True)
        data = to_avg(data)
        reref = '_avg'
    else:
        reref = ''

    if hp_filter is not None:
        hp_filt = Filter(low_cut=hp_filter, s_freq=data.s_freq)
        data = hp_filt(data)
    else:
        hp_filter = 0

    if lp_filter is not None:
        lp_filt = Filter(high_cut=lp_filter, s_freq=data.s_freq)
        data = lp_filt(data)
    else:
        lp_filter = 0

    if resample_freq is not None:
        res = Resample(s_freq=resample_freq)
        data = res(data)
    else:
        resample_freq = 0

    v = None
    if to_plot:
        TRIAL = int(data.number_of('trial') / 3)
        v = Viz1()  # size depends on n of channels
        v.add_data(data, trial=TRIAL, limits_y=(-100, 100))


    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               reref=reref, resample=resample_freq,
                               chan_types='-'.join(chan_type))

    with open(join(subj_dir, pkl_file), 'wb') as f:
        dump(data, f)

    with open(join(subj_dir, splitext(pkl_file)[0] + '_chan.txt'), 'w') as f:
        f.write('\n'.join(data.axis['chan'][0]))

    return v


def list_subj(period_name, chan_type=(), hp_filter=HP_FILTER,
              lp_filter=LP_FILTER, resample_freq=RESAMPLE_FREQ, reref=REREF):
    """Return list of subjects matching some parameters.

    Parameters
    ----------
    period_name : str
        period of interest
    chan_type : tuple of str
        list of channel groups of interest (such as 'grid', 'depth', 'scalp')
    hp_filter : float, optional
        high-pass filter cutoff
    lp_filter : float, optional
        low-pass filter cutoff
    resample_freq : int, optional
        frequency used for resampling
    reref : None or str, optional
        if data should be rereferenced to the average ('avg')

    Returns
    -------
    list of str
        list of subjects with matching parameters
    """
    subj = '*'
    if reref == 'avg':
        reref = '_avg'
    else:
        reref = ''
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq, reref=reref,
                               chan_types=''.join(chan_type))
    matching_files = glob(join(DATA_DIR, subj, REC_FOLDER, pkl_file))
    all_subj = sorted([x.split(sep)[-3] for x in matching_files])
    lg.info('SUBJECTS: ' + ', '.join(all_subj))

    return all_subj


def get_data(subj, period_name, chan_type=(), hp_filter=HP_FILTER,
             lp_filter=LP_FILTER, resample_freq=RESAMPLE_FREQ, reref=REREF):
    """Get the data for one subject quickly.

    Parameters
    ----------
    subj : str
        patient code


    Returns
    -------
    instance of DataTime

    """
    if hp_filter is None:
        hp_filter = 0
    if lp_filter is None:
        hp_filter = 0
    if resample_freq is None:
        resample_freq = 0

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq,
                               reref='',  # no reref, we do it online here
                               chan_types='-'.join(chan_type))

    lg.info('Subj %s, reading data: %s', subj, pkl_file)
    with open(join(subj_dir, pkl_file), 'rb') as f:
        data = load(f)

    chan = get_chan_used_in_analysis(subj, 'sleep', chan_type, reref='',
                                     resample_freq=resample_freq,
                                     hp_filter=hp_filter,
                                     lp_filter=lp_filter )[1]
    data.attr['chan'] = chan

    if reref == 'avg':
        montage = Montage(ref_to_avg=True)
    elif isinstance(reref, int):
        montage = Montage(bipolar=reref)
    else:
        montage = lambda x: x

    data = montage(data)

    return data


def get_hash(hp_filter=1, lp_filter=40, resample_freq=256, reref='',
             chan_type=()):
    """Get a (hopefully) unique code for each set of parameters. The risk of
    coincidence is really small.
    """
    s = REC_NAME[16:].format(hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                             resample=resample_freq, reref=reref,
                             chan_types='-'.join(chan_type))
    return md5(s.encode('utf-8')).hexdigest()[:4]


def get_chan_used_in_analysis(subj, period_name, chan_type=(),
                              hp_filter=HP_FILTER, lp_filter=LP_FILTER,
                              resample_freq=RESAMPLE_FREQ, reref=REREF):
    """Read quickly the channels used in the analysis

    Parameters
    ----------
    all_subj : str
        subject code to read

    Returns
    -------
    list of str
        list of selected channels in the actual dataset
    instance of Channels
        the channels for the patient of interest (but they might not have all
        the channels).
    """
    if hp_filter is None:
        hp_filter = 0
    if lp_filter is None:
        hp_filter = 0
    if resample_freq is None:
        resample_freq = 0
    if reref == 'avg':
        reref = '_avg'
    else:
        reref = ''

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq, reref=reref,
                               chan_types='-'.join(chan_type))

    good_chan = []
    with open(join(subj_dir, splitext(pkl_file)[0] + '_chan.txt'), 'r') as f:
        for one_chan in f:
            good_chan.extend(one_chan.splitlines())

    SESS = 'A'

    chan_file = join(REC_DIR, subj, 'doc', 'elec',
                     subj + '_elec_pos-names_sess' + SESS + '.csv')
    if not exists(chan_file):
        lg.debug('reading not-renamed electrodes')
        chan_file = join(REC_DIR, subj, 'doc', 'elec',
                         subj + '_elec_pos-adjusted_sess' + SESS + '.csv')

    try:
        chan = Channels(chan_file)

        chosen_chan = chan(lambda x: x.label in good_chan)
        lg.info('%s analysis chan %d, with location %d',
                subj, len(good_chan), chosen_chan.n_chan)
    except FileNotFoundError:
        lg.info('No channels for xltek datasets for %s', subj)
        chosen_chan = None

    return good_chan, chosen_chan


def _select_scores_per_subj(stages, duration, subj, choose='max'):
    """For each subject, find the dataset with the longest stages of interest

    Parameters
    ----------
    stages : tuple of str
        list of stages to analyze
    duration : int
        minimal duration to include an epoch, in seconds
    subj : str
        subject code
    choose : str, optional
        criterion to choose best night ('max': period with longest duration in
        stage, 'latest': the last possible period)

    Returns
    -------
    path to file
        scoring file with the longest duration in the stage(s) of interest.
    """
    score_dir = join(REC_DIR, subj, SCORE_PATH)
    all_scores = listdir(score_dir)

    if not all_scores:
        raise IndexError('No scores files for ' + subj)

    time_in_period = {}
    for one_score in all_scores:
        score = Annotations(join(score_dir, one_score))
        time_in_stages = sum(score.time_in_stage(x) for x in stages)
        lg.debug('    %s has % 5.1f min', one_score, time_in_stages / 60)
        if time_in_stages >= duration:
            time_in_period[one_score] = time_in_stages

    if choose == 'max':
        chosen_period = max(time_in_period, key=lambda k: time_in_period[k])
    elif choose == 'latest':
        # use the name to find the last one
        chosen_period = sorted(time_in_period.keys())[-1]

    return chosen_period


def _select_channels(chan_file, chan_type):
    """Read selected channels based on manual selection

    Parameters
    ----------
    chan_name : path to file
        file with the manually selected channels
    chan_type : tuple
        tuple of str to select channel groups, among 'depth', 'grid', 'scalp'

    Returns
    -------
    list of str
        list of channel labels belonging to the selected group.

    """
    with open(chan_file, 'r') as outfile:
        groups = json_load(outfile)

    chan = []
    for chan_grp in groups:
        if chan_grp['name'] in chan_type:
            chan.extend(chan_grp['chan_to_plot'])

    if not chan:
        raise ValueError(', '.join(chan_type) + ' not in selection')

    return chan
