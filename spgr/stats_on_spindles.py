from logging import getLogger
lg = getLogger('spgr')

from numpy import expand_dims, argsort, asarray, triu


def estimate_overlap(s, e):
    idx_sort = argsort(s)
    s = s[idx_sort]
    e = e[idx_sort]

    s = expand_dims(asarray(s), 1)
    e = expand_dims(asarray(e), 1)
    overlap = (s <= e.T) & (s.T <= e)
    return triu(overlap)
