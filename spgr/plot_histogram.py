from collections import Counter, defaultdict

from numpy import arange, histogram, zeros
from PyQt4.QtCore import Qt
from PyQt4.QtGui import (QBrush,
                         QPen,
                         )
from pyqtgraph import BarGraphItem
from phypno.viz import Viz1

from .constants import SPINDLE_OPTIONS
from .detect_spindles import get_spindles
from .lobes import LOBES, LOBE_COLORS
from .spindle_source import get_chan_with_regions
from .stats_on_spindles import create_spindle_groups


NoPen = QPen()
NoPen.setStyle(Qt.NoPen)


def make_hist_overlap(subj, color='kw', reref='avg', width=2, nchan=60):

    spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
    spindle_group = create_spindle_groups(spindles)

    chan = get_chan_with_regions(subj, reref)

    spindle_size = [len(x) for x in spindle_group]
    spindle_region = [_find_lobe_from_region(_find_region_for_spindle_group(x, chan)) for x in spindle_group]

    group_per_region = defaultdict(list)

    for one_size, one_region in zip(spindle_size, spindle_region):
        group_per_region[one_region].append(one_size)

    group_per_region = dict(group_per_region)

    hist = arange(0, nchan, width)
    y0 = zeros(hist.shape[0] - 1)

    v = Viz1(color='kw')
    p = v._widget.addPlot(title=subj)

    for lobe, x in group_per_region.items():

        h0, h1 = histogram(x, bins=hist)
        bars = BarGraphItem(x0=h1[:-1], y0=y0, height=h0, width=width,
                            pen=NoPen, brush=QBrush(LOBE_COLORS[lobe]))

        p.addItem(bars)
        y0 += h0

    return v


def _find_region_for_spindle_group(gr, chan):
    spindle_regions = []

    for one_chan in gr:
        sp_region = chan(lambda x: x.label == one_chan).return_attr('region')[0]
        if not sp_region.startswith('ctx-'):
            continue

        sp_region = sp_region[len('ctx-Xh-'):]
        spindle_regions.append(sp_region)

    if not spindle_regions:
        return 'Unspecified'
    else:
        c = Counter(spindle_regions)
        return c.most_common(1)[0][0]


def _find_lobe_from_region(region):
    for lobe_name, regions_in_lobe in LOBES.items():
        if region in regions_in_lobe:
            return lobe_name
    else:
        return 'unknown'
