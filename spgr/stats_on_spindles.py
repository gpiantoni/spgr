from logging import getLogger

from numpy import (argsort,
                   asarray,
                   diff,
                   expand_dims,
                   fill_diagonal,
                   NaN,
                   nanmean,
                   r_,
                   sum,
                   triu,
                   where,
                   zeros)

lg = getLogger('spgr')


def estimate_overlap(spindles):
    """Calculate overlap between single-channel spindles.

    Parameters
    ----------
    spindles : instance of phypno.graphoelement.Spindles
        spindles that were detected on each channel individually

    Returns
    -------
    bool matrix of # single-channel spindles
        true if two single-channel spindles occur at the same time. Matrix is
        not symmetric, the values are only in the lower diagonal.
    """
    singlechan_sp_0 = asarray([x['start_time'] for x in spindles.spindle])
    singlechan_sp_1 = asarray([x['end_time'] for x in spindles.spindle])

    idx_sort = argsort(singlechan_sp_0)
    singlechan_sp_0 = singlechan_sp_0[idx_sort]
    singlechan_sp_1 = singlechan_sp_1[idx_sort]

    singlechan_sp_0 = expand_dims(asarray(singlechan_sp_0), 1)
    singlechan_sp_1 = expand_dims(asarray(singlechan_sp_1), 1)
    overlap = ((singlechan_sp_0 <= singlechan_sp_1.T) &
               (singlechan_sp_0.T <= singlechan_sp_1))

    return triu(overlap)


def create_spindle_groups(spindles, param='chan'):
    """Group spindles based on the fact that they occur at the same time.

    Parameters
    ----------
    spindles : instance of phypno.graphoelement.Spindles
        spindles that were detected on each channel individually
    param : str
        take this spindle parameter for each group

    Returns
    -------
    list of list
        each element in list is one spindle group, which is defined by the
        channels in the list.
    """
    # make sure that the list of spindles (as dict) is sorted as well
    spindles.spindle = sorted(spindles.spindle, key=lambda x: x['start_time'])

    group_sep = _determine_spindle_group(spindles)

    sp_groups = []
    i0 = 0
    i1 = 0

    for one_sep in group_sep:
        i1 += one_sep
        sp_groups.append(spindles.spindle[i0:i1])
        i0 += one_sep

    chan_group = []
    for one_group in sp_groups:
        chan_group.append([x[param] for x in one_group])

    chan_group = [x for x in chan_group if len(x) != 0]

    return chan_group


def _determine_spindle_group(spindles):
    """Estimate the spindle groups.

    Parameters
    ----------
    spindles : instance of phypno.graphoelement.Spindles
        spindles that were detected on each channel individually

    Returns
    -------
    vector
        vector with indices that divide the single-spindles in groups.
    """
    overlap = estimate_overlap(spindles)
    overlap_sum = sum(overlap, axis=1)
    group_sep = diff(r_[1, where(overlap_sum == 1)[0]])[1:]

    return group_sep


def cooccur_likelihood(chan, spindle_group, n_other_channels):
    """Create nChan x nChan matrix with the likelihood to see a spindle at the
    same time in two channels, normalized by the total number of spindles of
    each channel.

    Parameters
    ----------
    chan : instance of phypno.attr.Channels
        channels to create the matrix
    spindle_group : list of list
        each element in list is one spindle group, which is defined by the
        channels in the list.
    n_other_channels : str
        number of other channels that need to have a spindle

    Returns
    -------
    nChan vector
        TODO
    """
    n_other_channels = int(n_other_channels)

    all_chan = chan.return_label()

    chan_prob = zeros((chan.n_chan, 2))

    for one_group in spindle_group:
        for chan0 in one_group:
            i0 = all_chan.index(chan0)
            if len(one_group) > n_other_channels:  # 1 -> the same channel
                chan_prob[i0, 1] += 1
            else:
                chan_prob[i0, 0] += 1

    return chan_prob[:, 1] / sum(chan_prob, axis=1)


def ratio_spindles_with_chan(chan, spindle_group):
    chan_prob = zeros((chan.n_chan, chan.n_chan))
    all_chan = chan.return_label()

    for chan0 in all_chan:
        for chan1 in all_chan:
            subgroup = [x for x in spindle_group if chan0 in x]
            sp_ratio = sum(1 for sp_group in subgroup if chan1 in sp_group) / len(subgroup)
            i0 = all_chan.index(chan0)
            i1 = all_chan.index(chan1)
            chan_prob[i0, i1] = sp_ratio

    fill_diagonal(chan_prob, NaN)

    return nanmean(chan_prob, 0)


def mean_spindle_cooccurrence(chan_prob_n, normalized_by='source'):
    """Take the mean of the probability of spindle cooccurrence.

    Parameters
    ----------
    nChan x nChan matrix
        matrix with the normalized co-occurrences of spindles between two
        channels.
    normalized_by : str
        'source' or 'target'

    Returns
    -------
    vector
        vector with the normalized values

    Notes
    -----
    If you normalize by 'source', then it will give you the mean probability of
    finding a spindle in that channel, normalized by the total number of
    spindles in that channel.
    If you normalize by 'target', it gives you the mean probability that one of
    the spindles in the channel of interest is important for the other
    channels.
    """
    if normalized_by == 'source':
        return nanmean(chan_prob_n, axis=0)
    if normalized_by == 'target':
        return nanmean(chan_prob_n, axis=1)


    chan_prob = zeros((chan.n_chan, chan.n_chan))
    all_chan = chan.return_label()
    for chan0 in all_chan:
        for chan1 in all_chan:
            subgroup = [x for x in spindle_group if chan0 in x]
            sp_ratio = sum(1 for sp_group in subgroup if chan1 in sp_group) / len(subgroup)
            i0 = all_chan.index(chan0)
            i1 = all_chan.index(chan1)
            chan_prob[i0, i1] = sp_ratio

