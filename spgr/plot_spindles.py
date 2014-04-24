from numpy import zeros, sum
from matplotlib.pyplot import subplots


SUBPLOT_ROW = 3
SUBPLOT_COL = 2

SUBPLOT_HEIGHT = 6
SUBPLOT_WIDTH = 6


def hist_values(all_subj, all_spindles,
                get_value, x_lim, nbin):

    f, subp = subplots(SUBPLOT_ROW, SUBPLOT_COL,
                       figsize=(SUBPLOT_HEIGHT * SUBPLOT_ROW,
                                SUBPLOT_WIDTH * SUBPLOT_COL))

    for ax, subj, spindles in zip(subp.flatten(), all_subj, all_spindles):
        value = [get_value(x) for x in spindles]
        value = [x for x in value if x is not None]
        ax.hist(value, range=x_lim,
                bins=(x_lim[1] - x_lim[0]) * nbin,
                align='left')
        ax.set_title(subj)
        ax.set_xlim(x_lim)


def channel_count(all_subj, all_spindles):

    f, subp = subplots(SUBPLOT_ROW, SUBPLOT_COL,
                       figsize=(SUBPLOT_HEIGHT * SUBPLOT_ROW,
                                SUBPLOT_WIDTH * SUBPLOT_COL))

    for ax, subj, spindles in zip(subp.flatten(), all_subj, all_spindles):

        all_chan = set([sp['chan'] for sp in spindles])
        n_spindles = zeros(len(all_chan))
        for i, i_chan in enumerate(all_chan):
            n_spindles[i] = sum([1 for sp in spindles if sp['chan'] == i_chan])

        ax.bar(range(len(n_spindles)), n_spindles)
        ax.set_title(subj)
