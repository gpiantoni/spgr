import warnings

from numpy import zeros, min, max, mean, meshgrid, linspace, NaN, isinf, asarray, sum
from scipy.interpolate import griddata
from matplotlib.pyplot import subplots, colorbar

from .stats_on_spindles import estimate_overlap

SUBPLOT_ROW = 3
SUBPLOT_COL = 2

SUBPLOT_HEIGHT = 6
SUBPLOT_WIDTH = 6

RESOLUTION = 200


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


def topo_values(all_subj, all_chan, all_spindles, get_value, take_mean,
                v_lim=None):

    f, subp = subplots(SUBPLOT_ROW, SUBPLOT_COL,
                       figsize=(SUBPLOT_HEIGHT * SUBPLOT_ROW,
                                SUBPLOT_WIDTH * SUBPLOT_COL))

    for ax, subj, chan, spindles in zip(subp.flatten(), all_subj, all_chan,
                                        all_spindles):

        values = zeros(chan.n_chan)
        n_spindles = zeros(chan.n_chan)
        for i, i_chan in enumerate(chan.chan):
            for sp in spindles:
                if sp['chan'] == i_chan.label:
                    n_spindles[i] += 1
                    values[i] += get_value(sp)

        if take_mean:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                values /= n_spindles
                values[isinf(values)] = 0

        xyz = chan.return_xyz()

        xy = xyz[:, 1:]
        min_xy = min(xy, axis=0)
        max_xy = max(xy, axis=0)

        x_grid, y_grid = meshgrid(linspace(min_xy[0], max_xy[0], RESOLUTION),
                                  linspace(min_xy[1], max_xy[1], RESOLUTION))

        zi = griddata(xy, values, (x_grid, y_grid), method='linear')

        if v_lim is None:
            vlim = (0, max(values))
        else:
            vlim = v_lim

        img = ax.imshow(zi, vmin=vlim[0], vmax=vlim[1], origin='lower',
                        aspect='equal',
                        extent=[min_xy[0], max_xy[0], min_xy[1], max_xy[1]])

        ax.set_title(subj)
        colorbar(img, ax=ax)

        if mean(xyz[:, 0]) < 0:
            ax.invert_xaxis()


def hist_overlap(all_subj, all_spindles):
    f, subp = subplots(3, 2, figsize=(18, 12))

    x_lim = (0, 60)

    for ax, subj, spindles in zip(subp.flatten(), all_subj, all_spindles):
        s = asarray([x['start_time'] for x in spindles])
        e = asarray([x['end_time'] for x in spindles])
        ov = estimate_overlap(s, e)
        ax.hist(sum(ov, axis=1), bins=(x_lim[1] - x_lim[0]),
                range=x_lim, align='left')
        ax.set_title(subj)
        ax.set_xlim(x_lim)
