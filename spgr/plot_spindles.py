from tempfile import mkstemp, mkdtemp
tmpdir = mkdtemp()

from visvis import subplot, figure, hist, screenshot, title

SUBPLOT_ROW = 3
SUBPLOT_COL = 2

SUBPLOT_HEIGHT = 400
SUBPLOT_WIDTH = 400

RESOLUTION = 200


def plot_inline(fig):
    """Easy work-around to plot figures in

    Parameters
    ----------
    fig : instance of visvis.BaseFigure
        figure to plot inline

    """
    fig.relativeFontSize = 1
    fig.DrawNow()
    plot_file = mkstemp(suffix='.png', dir=tmpdir)[1]
    screenshot(plot_file, sf=1, bg='w')
    fig.Destroy()

    return plot_file


def hist_values(all_subj, all_spindles, get_value, x_lim):
    """Plot an histogram of the values.

    Parameters
    ----------
    all_subj : list of str
        list with the name of the subjects (for the title)
    all_spindles : list
        it can be a list of spindles or of values
    get_value : funct
        function to get values from spindles or values
    x_lim : tuple of two float
        min and max values for the x-axis

    """
    f = figure()
    f.position = (0, 0, SUBPLOT_HEIGHT * SUBPLOT_ROW,
                  SUBPLOT_WIDTH * SUBPLOT_COL)

    for i, spindles in enumerate(all_spindles):
        subplot(SUBPLOT_ROW, SUBPLOT_COL, i + 1)

        if isinstance(spindles, dict):
            value = get_value(spindles)
        else:
            value = [get_value(x) for x in spindles]
            value = [x for x in value if x is not None]

        hist(value, drange=x_lim, bins=100)
        title(all_subj[i])

    return plot_inline(f)
