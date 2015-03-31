from numpy import asarray, sum, where, diff, r_, histogram, arange, NaN, nanmean, mean, isnan

from phypno.attr import Freesurfer
from phypno.viz import Viz3

from .constants import DEFAULT_HEMI, FS_AVG
from .stats_on_spindles import estimate_overlap


fs = Freesurfer(FS_AVG)
surf_avg = getattr(fs.read_brain(), DEFAULT_HEMI)


def hist_overlap(spindles, width=2, nchan=70):

    s = asarray([x['start_time'] for x in spindles.spindle])
    e = asarray([x['end_time'] for x in spindles.spindle])
    ov = estimate_overlap(s, e)
    x = sum(ov, axis=1)

    v = diff(r_[1, where(x == 1)[0]])

    hist = arange(0, nchan, width)
    return histogram(v, hist)


def plot_surf(all_values, threshold=(None, None), limits=(0, 2),
              extra_smoothing=True):
    """Plot values onto the surface.

    Parameters
    ----------
    all_values : list of Data
        values for each subject
    threshold : tuple of 2 float or of None
        low and high thresholds to include the values
    limits : tuple of 2 floats
        values used for color scaling

    Returns
    -------
    instance of Viz3
        plot with the surfaces
    """
    for x in all_values:
        if threshold[0]:
            x.data[0][x.data[0] < threshold[0]] = NaN
        if threshold[1]:
            x.data[0][x.data[0] > threshold[1]] = NaN

    values = nanmean(asarray([x.data[0] for x in all_values]), axis=0)

    if extra_smoothing:
        # apply some quick smoothing (but rather strong, useful to avoid the clown fish effect)
        for one_tri in surf_avg.tri:
            if not any(isnan(values[one_tri])):
                values[one_tri] = mean(values[one_tri])

    v = Viz3(color='kw')
    v.add_surf(surf_avg, values=values, limits_c=limits,
               color=(100, 100, 100, 255))
    v._widget.opts.update({'elevation': 15, 'azimuth': 17})

    return v
