from numpy import arange, ceil, concatenate, histogram, linspace
from vispy.geometry import Rect
from vispy.scene.visuals import Rectangle

from phypno.viz import Viz1

from .constants import DPI, SPINDLE_OPTIONS, TICKS_FONT_SIZE
from .detect_spindles import get_spindles
from .stats_on_spindles import count_sp_at_any_time


from .constants import CHAN_TYPE, DATA_OPTIONS
from .read_data import get_data

X_MAJOR_TICK = 10
X_MINOR_TICK = 1
Y_MAJOR_TICK = .1
Y_MINOR_TICK = .05

FIG_SIZE = 40, 30  # actual size in mm for figure


def make_hist_overlap(subj, reref='avg', width=2, nchan=60):

    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
    data = get_data(subj, period_name='sleep', chan_type=CHAN_TYPE,
                    reref=reref, **DATA_OPTIONS)
    t_range = concatenate(data.axis['time'][:])

    p = count_sp_at_any_time(spindles, t_range)
    p_with_sp = p[p >= 1]
    hist, bin_edges = histogram(p_with_sp, bins=arange(.5, nchan, width))
    h_chan = hist * arange(1, nchan, width)  # TODO: normalization
    hist = h_chan / sum(h_chan)  # TODO: rename hist

    v = Viz1(dpi=DPI, size_mm=FIG_SIZE)
    plt = v._fig[0, 0]
    plt._configure_2d()

    for left_edge, col in zip(bin_edges[:-1], hist):
        if col > 0:  # vispy doesn't like it when hist has height = 0
            rect = Rectangle(center=(left_edge + width / 2, col / 2),
                             height=col, width=width, border_color='k',
                             color=None)
            rect.border.set_data(width=2)  # border_width doesn't work
            plt.view.add(rect)

    max_height = ceil(max(hist) / Y_MAJOR_TICK) * Y_MAJOR_TICK
    plt.view.camera.set_state(rect=Rect(pos=(-0.1, -1),
                                        size=(nchan, max_height)))
    """
    def xtick_frac():
        tick_labels = [str(x) for x in arange(0, nchan + X_MAJOR_TICK, X_MAJOR_TICK)]
        major_tick_fractions = linspace(0, 1, len(tick_labels))
        minor_tick_fractions = linspace(0, 1,  nchan / X_MINOR_TICK + 1)

        return major_tick_fractions, minor_tick_fractions, tick_labels

    plt.xaxis.axis.ticker._get_tick_frac_labels = xtick_frac
    plt.xaxis.axis._text.font_size = TICKS_FONT_SIZE

    def ytick_frac():
        tick_labels = [str(x) for x in arange(0, max_height + Y_MAJOR_TICK, Y_MAJOR_TICK)]
        major_tick_fractions = linspace(0, 1, len(tick_labels))
        minor_tick_fractions = linspace(0, 1,  max_height / Y_MINOR_TICK + 1)

        return major_tick_fractions, minor_tick_fractions, tick_labels

    plt.yaxis.axis.ticker._get_tick_frac_labels = ytick_frac
    plt.yaxis.axis._text.font_size = TICKS_FONT_SIZE
    """
    plt.margin = 25  # otherwise xtick label overlaps with border
    plt.view.border_color = 'w'

    return v
