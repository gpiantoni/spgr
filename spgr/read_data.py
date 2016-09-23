from datetime import timedelta
from glob import glob
from json import load as json_load
from logging import getLogger
from os import sep
from os.path import join, splitext
from pickle import dump, load

from phypno import Dataset
from phypno.attr import Annotations, Channels
from phypno.trans import Filter, Montage, Resample
from phypno.trans.montage import create_bipolar_chan

from .constants import (DATA_PATH,
                        DATA_OPTIONS,
                        REC_PATH,
                        SCORES_PATH,
                        REC_FOLDER,
                        ELEC_FOLDER,
                        XLTEK_FOLDER,
                        SESSION)

lg = getLogger(__name__)


REREF = ''
REC_NAME = ('{subj}_{period}_{chan_types}_hp{hp:03d}_lp{lp:03d}_'
            'rs{resample:03d}.pkl')


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
            lg.info(subj + ' has ' + str(scores[subj]))
        except IndexError:
            lg.debug(subj + ': no scored recordings')
        except ValueError:
            lg.info(subj + ': duration not sufficient')

    return scores


def save_data(subj, score_file, period_name, stages, chan_type=(),
              hp_filter=DATA_OPTIONS['hp_filter'],
              lp_filter=DATA_OPTIONS['lp_filter'],
              resample_freq=DATA_OPTIONS['resample_freq']):
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
    chan_type : tuple
        tuple of str to select channel groups, among 'depth', 'grid', 'scalp'

    Returns
    -------
    int
        number of channels
    float
        duration of the selected recordings in min

    Notes
    -----
    It saves the data to disk. The name of the file contains the parameters.
    """
    lg.info(subj + ' ' + str(score_file))
    subj_dir = DATA_PATH.joinpath(subj).joinpath(REC_FOLDER)
    if not subj_dir.is_dir():
        subj_dir.mkdir(parents=True)

    xltek_file = REC_PATH.joinpath(subj).joinpath(XLTEK_FOLDER).joinpath(score_file.name[:-11])
    d = Dataset(str(xltek_file))
    lg.info('Sampling Frequency {} Hz '.format(d.header['s_freq']))

    score = Annotations(str(score_file))
    chan_dir = REC_PATH.joinpath(subj).joinpath(ELEC_FOLDER)
    chan_file = chan_dir.joinpath(score_file.name.replace('_scores.xml',
                                                          '_channels.json'))

    selected_chan = _select_channels(chan_file, chan_type)
    lg.info('N Channels {} '.format(len(selected_chan)))
    lg.debug('Channels: ' + ', '.join(selected_chan))

    start_time = [x['start'] for x in score.epochs if x['stage'] in stages]
    end_time = [x['end'] for x in score.epochs if x['stage'] in stages]


    abs_time = d.header['start_time']
    abs_start = abs_time + timedelta(seconds=start_time[0])
    abs_end = abs_time + timedelta(seconds=end_time[-1])
    duration = sum([x1 - x0 for x0, x1 in zip(start_time, end_time)])
    lg.info('Start Time: ' + str(abs_start))
    lg.info('End Time: ' + str(abs_end))
    lg.info('Duration: {0: 3.1f} min'.format(duration / 60))

    data = d.read_data(begtime=start_time, endtime=end_time,
                       chan=selected_chan)

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

    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq,
                               chan_types='-'.join(chan_type))

    with open(str(subj_dir.joinpath(pkl_file)), 'wb') as f:
        dump(data, f)

    with open(str(subj_dir.joinpath(splitext(pkl_file)[0] +
                                    '_chan.txt')), 'w') as f:
        f.write('\n'.join(data.axis['chan'][0]))

    return len(selected_chan), duration


def list_subj(period_name, chan_type=(), hp_filter=DATA_OPTIONS['hp_filter'],
              lp_filter=DATA_OPTIONS['lp_filter'],
              resample_freq=DATA_OPTIONS['resample_freq']):
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

    Returns
    -------
    list of str
        list of subjects with matching parameters
    """
    subj = '*'
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq,
                               chan_types=''.join(chan_type))
    matching_files = glob(join(DATA_PATH, subj, REC_FOLDER, pkl_file))
    all_subj = sorted([x.split(sep)[-3] for x in matching_files])
    lg.info('SUBJECTS: ' + ', '.join(all_subj))

    return all_subj


def get_data(subj, period_name, chan_type=(),  reref=REREF,
             hp_filter=DATA_OPTIONS['hp_filter'],
             lp_filter=DATA_OPTIONS['lp_filter'],
             resample_freq=DATA_OPTIONS['resample_freq']):
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

    subj_dir = DATA_PATH.joinpath(subj).joinpath(REC_FOLDER)
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq,
                               chan_types='-'.join(chan_type))

    lg.info('Subj %s, reading data: %s', subj, pkl_file)
    with open(str(subj_dir.joinpath(pkl_file)), 'rb') as f:
        data = load(f)

    chan = get_chan_used_in_analysis(subj, 'sleep', chan_type, reref='',
                                     resample_freq=resample_freq,
                                     hp_filter=hp_filter,
                                     lp_filter=lp_filter)
    data.attr['chan'] = chan

    if reref == 'avg':
        montage = Montage(ref_to_avg=True)
    elif isinstance(reref, int):
        montage = Montage(bipolar=reref)
    else:
        montage = lambda x: x

    data = montage(data)
    chan = data.attr['chan']

    return data


def get_chan_used_in_analysis(subj, period_name, chan_type=(), reref=REREF,
                              hp_filter=DATA_OPTIONS['hp_filter'],
                              lp_filter=DATA_OPTIONS['lp_filter'],
                              resample_freq=DATA_OPTIONS['resample_freq'],
                              return_all_chan=False):
    """Read quickly the channels used in the analysis

    Parameters
    ----------
    all_subj : str
        subject code to read
    return_all_chan : bool
        whether it should return one argument (only channels used in the
        analysis) or two arguments (also all the possible channels)

    Returns
    -------
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

    subj_dir = DATA_PATH.joinpath(subj).joinpath(REC_FOLDER)
    pkl_file = REC_NAME.format(subj=subj, period=period_name,
                               hp=int(10 * hp_filter), lp=int(10 * lp_filter),
                               resample=resample_freq, reref='',
                               chan_types='-'.join(chan_type))

    good_chan = []
    with open(str(subj_dir.joinpath(splitext(pkl_file)[0] +
                                    '_chan.txt')), 'r') as f:
        for one_chan in f:
            good_chan.extend(one_chan.splitlines())

    SESS = SESSION(subj)

    chan_dir = REC_PATH.joinpath(subj).joinpath('doc').joinpath('elec')
    chan_file = chan_dir.joinpath(subj + '_elec_pos-names_sess' + SESS +
                                  '.csv')
    if not chan_file.exists():
        lg.debug('reading not-renamed electrodes')
        chan_file = chan_file.joinpath(subj + '_elec_pos-adjusted_sess' +
                                       SESS + '.csv')

    try:
        chan = Channels(str(chan_file))

        chosen_chan = chan(lambda x: x.label in good_chan)
        lg.debug('%s analysis chan %d, with location %d',
                 subj, len(good_chan), chosen_chan.n_chan)
    except FileNotFoundError:
        lg.info('Could not read file %s', chan_file)
        chosen_chan = None

    if isinstance(reref, float) or isinstance(reref, int):
        chosen_chan = create_bipolar_chan(chosen_chan, reref)[0]

    if return_all_chan:
        return chosen_chan, chan
    else:
        return chosen_chan


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
    all_scores = SCORES_PATH.rglob(subj + '_*')

    time_in_period = {}
    for one_score in all_scores:
        score = Annotations(str(one_score))
        time_in_stages = sum(score.time_in_stage(x) for x in stages)
        lg.debug('    %s has % 5.1f min', one_score.stem, time_in_stages / 60)
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
    with open(str(chan_file), 'r') as outfile:
        groups = json_load(outfile)

    chan = []
    for chan_grp in groups:
        if chan_grp['name'] in chan_type:
            chan.extend(chan_grp['chan_to_plot'])

    if not chan:
        raise ValueError(', '.join(chan_type) + ' not in selection')

    return chan
