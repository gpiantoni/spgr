from numpy import array, asarray, NaN, nanmean, nansum, mean, isnan, zeros
from vispy.color import get_colormap, ColorArray
import warnings

from phypno.viz import Viz3
from phypno.viz.base import normalize

from .constants import (AZIMUTH,
                        COLORMAP,
                        DPI,
                        ELEVATION,
                        NAN_COLOR,
                        avg_surf,
                        avg_vert,
                        avg_regions,
                        )

cm = get_colormap(COLORMAP)


def plot_surf(all_values, size_mm, limits=None, extra_smoothing=True,
              fun='mean'):
    """Plot values onto the surface.

    Parameters
    ----------
    all_values : list of Data
        values for each subject
    size_mm : tuple of 2 int
        size in pixels of the final image
    limits : tuple of 2 floats
        values used for color scaling
    extra_smoothing : bool
        if it should apply some extra smoothing
    fun : str
        'mean' or 'sum'

    Returns
    -------
    instance of Viz3
        plot with the surfaces
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # some columns are completely NaN
        if fun == 'mean':
            values = nanmean(asarray([x.data[0] for x in all_values]), axis=0)
        elif fun == 'sum':
            values = nansum(asarray([x.data[0] for x in all_values]), axis=0)

    if extra_smoothing:
        # apply some quick smoothing (but rather strong, useful to avoid the clown fish effect)
        for one_tri in avg_surf.tri:
            if not any(isnan(values[one_tri])):
                values[one_tri] = mean(values[one_tri])

    v = Viz3(dpi=DPI, show=False, size_mm=size_mm)
    v.add_surf(avg_surf, values=values, limits_c=limits, color=NAN_COLOR,
               colormap=COLORMAP)
    v._plt.view.camera.elevation = ELEVATION
    v._plt.view.camera.azimuth = AZIMUTH
    return v


def plot_lmer(coef, size_mm, pvalues=None, limits=(0, 2), p_threshold=0.05):
    """Plot coefficients computed from lmer.

    Parameters
    ----------
    coef : dict
        dictionary with coefficients (computed with lmer)
    size_mm : tuple of 2 int
        size in pixels of the final image
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
    val = zeros((avg_vert.shape[0], 4))
    val.fill(NaN)

    for one_region, one_v in coef.items():
        idx = array(avg_vert) == avg_regions.index(one_region)
        norm_v = normalize(one_v, limits)
        if pvalues[one_region] < 0.05:
            level = 1
        else:
            level = 0.4
        val[idx, :] = saturate(cm[norm_v], level).rgba

    hasnan = isnan(val).all(axis=1)
    val[hasnan, :] = NAN_COLOR

    for one_region, one_v in coef.items():
        if pvalues[one_region] <= p_threshold:
            val[array(avg_vert) == avg_regions.index(one_region)] = one_v

    v = Viz3(dpi=DPI, show=False, size_mm=size_mm)
    v.add_surf(avg_surf, vertex_colors=val, color=NAN_COLOR)
    v._plt.view.camera.elevation = ELEVATION
    v._plt.view.camera.azimuth = AZIMUTH
    return v


def saturate(c, level=0.5):
    """Only works if colormap has constant saturation"""
    c = c.hsv
    c[0, 1] = level
    return ColorArray(color=c, color_space='hsv')
