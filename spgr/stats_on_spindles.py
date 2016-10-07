from logging import getLogger

from numpy import (arange,
                   concatenate,
                   mean,
                   median,
                   percentile,
                   where,
                   zeros)
from scipy.stats import ttest_rel

from .constants import (DATA_OPTIONS,
                        PARAMETERS,
                        SPINDLE_OPTIONS)
from .detect_spindles import get_spindles
from .read_data import keep_time_chan

lg = getLogger('spgr')
PERCENT = PARAMETERS['PERCENTILE']
S_FREQ = DATA_OPTIONS['resample_freq']


def count_sp_at_any_time(sp, t_range):
    """At any given time point, check how many spindles you observe in all the
    channels.

    Parameters
    ----------
    sp : instance of Spindles
        spindles to analyze
    t_range : ndarray vector
        vector of actual time point in the recordings

    Returns
    -------
    ndarray vector
        same size as t_range, for each time point it tells you how many
        spindles there are
    """
    t_in = zeros(t_range.shape, dtype=int)
    for one_sp in sp.spindle:
        t_in += ((t_range >= one_sp['start_time']) &
                 (t_range < one_sp['end_time'])).astype(int)

    return t_in


def count_cooccur_per_chan(subj, reref, summarize='mean'):
    """At any given time point, for each channel count how many other channels
    have a spindle. This creates a distribution (x-axis: number of spindles,
    y-axis: number of channels with the number of spindles) where each point is
    a time point. Then you can summarize the distribution by taking the mean
    or the median.

    Parameters
    ----------
    subj : str
        subject code
    reref : str or float
        'avg' or 15, for average reference or bipolar montage

    Returns
    -------
    ndarray vector
        of n_chan length, where each value represents the mean of the
        distribution of spindles that cooccur.
    """
    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)

    time, chan = keep_time_chan(subj, reref)
    t_range = concatenate(time[:])
    chan_list = chan[0]

    p = count_sp_at_any_time(spindles, t_range)
    t_range_with_sp = t_range[p >= 1]
    p_with_sp = p[p >= 1]
    chan_prob = zeros(len(chan_list))

    for i, one_chan in enumerate(chan_list):
        t_at_sp = zeros(t_range_with_sp.shape, dtype=bool)
        for one_sp in spindles.spindle:
            if one_sp['chan'] == one_chan:
                t_at_sp = (t_at_sp |
                           ((t_range_with_sp >= one_sp['start_time']) &
                            (t_range_with_sp < one_sp['end_time'])))

        if summarize == 'mean':
            chan_prob[i] = mean(p_with_sp[t_at_sp])

        elif summarize == 'median':
            chan_prob[i] = median(p_with_sp[t_at_sp])

    return chan_prob


def get_cooccur_percent(subj, reref, lg):
    """Get values for different spindle parameters for most isolated and most
    cooccurring spindles.

    Parameters
    ----------
    subj : str
        subject code
    reref : str
        re-reference type
    lg : Logger
        logger to write to

    Returns
    -------
    dict
        summary parameters for isolated spindles
    dict
        summary parameters for cooccurring spindles
    """
    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
    t = [x['start_time'] for x in spindles.spindle]
    t_range = arange(min(t), max(t), 1 / S_FREQ)
    p = count_sp_at_any_time(spindles, t_range)

    lg.info('{}'.format(subj))
    df_i = _compute_percent(lg, spindles, t_range, p, 'isolated')
    df_c = _compute_percent(lg, spindles, t_range, p, 'cooccurring')

    return df_i, df_c


def _compute_percent(lg, spindles, t_range, p, sp_type):
    """Compute summary parameters for isolated and cooccurring spindles

    Parameters
    ----------
    lg : Logger
        logger to write to
    spindles : instance of Spindles
        spindles for one specific subject
    t_range : ndarray
        vector with all the possible time points (even those with no spindles,
        it doesn't matter)
    p : ndarray (same length as t_range)
        number of spindles in each time point
    sp_type : str
        'isolated' or 'cooccurring'

    Returns
    -------
    dict
        summary parameters for isolated spindles
    dict
        summary parameters for cooccurring spindles
    """
    p_with = p[p >= 1]
    if sp_type == 'isolated':
        p_pool = where((p <= percentile(p_with, PERCENT)) & p >= 1)[0]
    elif sp_type == 'cooccurring':
        p_pool = where((p >= percentile(p_with, 100 - PERCENT)))[0]
    i_t = t_range[p_pool]

    all_sp = []
    for sp in spindles.spindle:
        if any((sp['start_time'] <= i_t) & (sp['end_time'] >= i_t)):
            all_sp.append(sp)

    lg.info('Number of {} spindles: {}'.format(sp_type, len(all_sp)))

    df = {'freq': mean([x['peak_freq'] for x in all_sp]),
          'ampl': mean([x['peak_val'] for x in all_sp]),
          'dur': mean([x['end_time'] - x['start_time'] for x in all_sp]),
          }

    return df


def print_table_percent(lg, df_i, df_c):
    """Print summary tables for statistics on the top and bottom cooccurring
    spindles.

    Parameters
    ----------
    lg : Logger
        logger to write to
    df_i : list of dict
        list of summary statistics for each subjects in the parameters of
        interest for the most isolated spindles
    df_c : list of dict
        list of summary statistics for each subjects in the parameters of
        interest for the most cooccurring spindles
    """
    params = 'freq', 'ampl', 'dur'

    for param in params:
        i_val = [x[param] for x in df_i]
        c_val = [x[param] for x in df_c]
        tstat, pvalue = ttest_rel(i_val, c_val)
        lg.info('with param{:>5}, {:d}% most isolated spindles: {:.3f}, '
                '{:d}% most cooccurring spindles: {:.3f}\n'
                'paired t({:d}) = {:.3f}, p-value = {:.3f}'
                ''.format(param, PERCENT, mean(i_val), PERCENT, mean(c_val),
                          len(i_val) - 1, tstat, pvalue))
