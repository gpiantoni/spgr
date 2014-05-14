from logging import getLogger
lg = getLogger('spgr')

from numpy import expand_dims, argsort, asarray, sum, triu
from numpy.random import uniform, seed

from multiprocessing import Pool

from .read_data import MIN_EPOCHS  # it should also be MAX_EPOCHS

DUR_EPOCHS = 30
total_duration = MIN_EPOCHS * DUR_EPOCHS

N_RPT = 1000


def estimate_overlap(s, e):
    idx_sort = argsort(s)
    s = s[idx_sort]
    e = e[idx_sort]

    s = expand_dims(asarray(s), 1)
    e = expand_dims(asarray(e), 1)
    overlap = (s <= e.T) & (s.T <= e)
    return triu(overlap)

def fake_dist(i, s, e):
    seed(i)
    fake_s = uniform(0, total_duration, len(s))
    fake_e = fake_s + (e - s)
    return sum(estimate_overlap(fake_s, fake_e))

def bootstrap_overlap(all_subj, all_spindles):
    for subj, spindles in zip(all_subj, all_spindles):
        s = asarray([x['start_time'] for x in spindles])
        e = asarray([x['end_time'] for x in spindles])
        ov = sum(estimate_overlap(s, e))

        def fake_dist_subj(x):
            return fake_dist(x, s, e)

        dist = list(map(fake_dist_subj, range(N_RPT)))
        lg.info('probability for {0}: {1:.9f}'.format(subj,
                                                      sum(ov <= dist) / N_RPT))
