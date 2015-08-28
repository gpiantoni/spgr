from numpy import array, asarray, NaN, nanmean, mean, isnan, zeros
import warnings

from phypno.attr import Freesurfer
from phypno.viz import Viz3

from .constants import (AZIMUTH,
                        DEFAULT_HEMI,
                        ELEVATION,
                        FS_AVG,
                        NAN_COLOR,
                        PLOT_COLOR)

fs = Freesurfer(FS_AVG)
surf_avg = getattr(fs.read_brain(), DEFAULT_HEMI)
label_avg = fs.read_label(DEFAULT_HEMI)


def plot_surf(all_values, threshold=(None, None), limits=None,
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # some columns are completely NaN
        values = nanmean(asarray([x.data[0] for x in all_values]), axis=0)

    if extra_smoothing:
        # apply some quick smoothing (but rather strong, useful to avoid the clown fish effect)
        for one_tri in surf_avg.tri:
            if not any(isnan(values[one_tri])):
                values[one_tri] = mean(values[one_tri])

    v = Viz3(color=PLOT_COLOR)
    v.add_surf(surf_avg, values=values, limits_c=limits, color=NAN_COLOR)
    v._canvas.view.camera.elevation = ELEVATION
    v._canvas.view.camera.azimuth = AZIMUTH
    return v


def plot_lmer(coef, pvalues=None, limits=(0, 2), p_threshold=0.05):
    """Plot coefficients computed from lmer.

    Parameters
    ----------
    coef : dict
        dictionary with coefficients (computed with lmer)
    pvalues : dict, optional
        dictionary with pvalues
    limits : tuple of 2 floats
        values used for color scaling
    p_threshold : float, optional
        threshold for the pvalues

    Returns
    -------
    instance of Viz3
        plot with the surfaces
    """
    val = zeros(label_avg[0].shape)
    val.fill(NaN)

    for one_region, one_v in coef.items():
        if pvalues[one_region] <= p_threshold:
            val[array(label_avg[0]) == label_avg[2].index(one_region)] = one_v

    v = Viz3(color=PLOT_COLOR)
    v.add_surf(surf_avg, values=val, limits_c=limits,
               color=(100, 100, 100, 255))
    v._widget.opts.update({'elevation': 15, 'azimuth': 17})
    return v
