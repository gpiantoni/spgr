from numpy import arange, concatenate, histogram, linspace, mean
from vispy.geometry import Rect
from vispy.scene.visuals import Rectangle

from phypno.viz import Viz1

from .constants import (CHAN_TYPE,
                        DATA_OPTIONS,
                        DPI,
                        HIST_FIG_SIZE,
                        HIST_WIDTH,
                        HIST_BAR_COLOR,
                        HIST_BAR_WIDTH,
                        HIST_N_CHAN,
                        SPINDLE_OPTIONS,
                        TICKS_FONT_SIZE)
from .detect_spindles import get_spindles
from .read_data import get_data
from .stats_on_spindles import count_sp_at_any_time


X_MAJOR_TICK = 10
X_MINOR_TICK = 1
Y_MAX = 50
Y_MAJOR_TICK = 10
Y_MINOR_TICK = 1


def make_hist_overlap(subj, reref='avg'):

    nchan = HIST_N_CHAN[reref]

    width=HIST_WIDTH

    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
    data = get_data(subj, period_name='sleep', chan_type=CHAN_TYPE,
                    reref=reref, **DATA_OPTIONS)
    t_range = concatenate(data.axis['time'][:])

    p = count_sp_at_any_time(spindles, t_range)
    p_with_sp = p[p >= 1]
    h_chan, bin_edges = histogram(p_with_sp, bins=arange(.5, nchan, width))

    # normalization (so that are == 1)
    hist_norm = h_chan / sum(h_chan) * 100

    v = Viz1(dpi=DPI, size_mm=HIST_FIG_SIZE)
    plt = v._fig[0, 0]
    plt._configure_2d()

    for left_edge, col in zip(bin_edges[:-1], hist_norm):
        if col > 0:  # vispy doesn't like it when hist has height = 0
            rect = Rectangle(center=(left_edge + width / 2, col / 2),
                             height=col, width=width, border_color='k',
                             color=HIST_BAR_COLOR)
            rect.border.set_data(width=HIST_BAR_WIDTH)  # border_width doesn't work
            plt.view.add(rect)

    plt.view.camera.set_state(rect=Rect(pos=(0, 0), size=(nchan, Y_MAX)))

    def xtick_frac():
        tick_labels = [str(x) for x in arange(0, nchan + X_MAJOR_TICK,
                                              X_MAJOR_TICK)]
        major_tick_fractions = linspace(0, 1, len(tick_labels))
        minor_tick_fractions = linspace(0, 1,  nchan / X_MINOR_TICK + 1)

        return major_tick_fractions, minor_tick_fractions, tick_labels

    plt.xaxis.axis.ticker._get_tick_frac_labels = xtick_frac
    plt.xaxis.axis._text.font_size = TICKS_FONT_SIZE

    def ytick_frac():
        tick_labels = [str(x) for x in arange(0, Y_MAX + Y_MAJOR_TICK,
                                              Y_MAJOR_TICK)]
        major_tick_fractions = linspace(0, 1, len(tick_labels))
        minor_tick_fractions = linspace(0, 1, Y_MAX / Y_MINOR_TICK + 1)

        return major_tick_fractions, minor_tick_fractions, tick_labels

    plt.yaxis.axis.ticker._get_tick_frac_labels = ytick_frac
    plt.yaxis.axis._text.font_size = TICKS_FONT_SIZE

    plt.margin = 25  # otherwise xtick label overlaps with border
    plt.view.border_color = 'w'

    return v, mean(p_with_sp)
