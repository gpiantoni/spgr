from logging import getLogger

from numpy import (concatenate,
                   mean,
                   zeros)

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        SPINDLE_OPTIONS)
from .detect_spindles import get_spindles
from .read_data import get_data

lg = getLogger('spgr')


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


def count_cooccur_per_chan(subj, reref):
    """At any given time point, for each channel count how many other channels
    have a spindle. This creates a distribution (x-axis: number of spindles,
    y-axis: number of channels with the number of spindles) where each point is
    a time point. Then you can summarize the distribution by taking the mean.

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
    spindles = get_spindles(subj, reref, **SPINDLE_OPTIONS)

    data = get_data(subj, period_name='sleep', chan_type=CHAN_TYPE,
                    reref=reref, **DATA_OPTIONS)
    t_range = concatenate(data.axis['time'][:])
    chan_list = data.axis['chan'][0]

    p = count_sp_at_any_time(spindles, t_range)
    t_range_with_sp = t_range[p >= 1]
    p_with_sp = p[p >= 1]
    chan_prob = zeros(len(chan_list))

    for i, one_chan in enumerate(chan_list):
        t_at_sp = zeros(t_range_with_sp.shape, dtype=bool)
        for one_sp in spindles.spindle:
            if one_sp['chan'] == chan_list:
                t_at_sp = (t_at_sp |
                           ((t_range_with_sp >= one_sp['start_time']) &
                            (t_range_with_sp < one_sp['end_time'])))

        chan_prob[i] = mean(p_with_sp[t_at_sp])

    return chan_prob
