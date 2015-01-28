from numpy import asarray, sum, where, diff, r_, histogram, arange

from .stats_on_spindles import estimate_overlap


def hist_overlap(spindles, width=2, nchan=70):

    s = asarray([x['start_time'] for x in spindles.spindle])
    e = asarray([x['end_time'] for x in spindles.spindle])
    ov = estimate_overlap(s, e)
    x = sum(ov, axis=1)

    v = diff(r_[1, where(x == 1)[0]])

    hist = arange(0, nchan, width)
    return histogram(v, hist)
