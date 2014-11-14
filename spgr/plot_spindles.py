from os.path import join, splitext, exists
from tempfile import mkstemp, mkdtemp

from matplotlib.pyplot import bar
from numpy import mean, asarray, sum, where, diff, r_, histogram, arange
from visvis import colorbar, CM_HOT, subplot, figure, hist, screenshot, title

from phypno.attr import Freesurfer
# from phypno.viz.plot_3d import plot_surf, plot_chan, make_movie

from .read_data import REC_DIR, GROUP_DIR, get_chan_used_in_analysis
from .stats_on_spindles import estimate_overlap

from base64 import b64encode
from IPython.display import HTML, Image

tmpdir = mkdtemp()


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


def plot_topography(all_values, limits=None, colormap=CM_HOT):
    """Plot topography as 3D video, from values to html tag.

    Parameters
    ----------
    all_values : dict
        where keys are subject codes, and values is a vector with the data
        for each electrode
    limits : tuple of float, optional
        min and max values for the plots
    colormap : ndarray, optional
        colormap from visvis

    Returns
    -------
    str
        long string to pass to IPython HTML

    """
    all_movies = []
    for subj, values in all_values.items():
        movie_name = create_topography_movie(subj, values, limits=limits,
                                             colormap=colormap)
        all_movies.append(movie_name)

    return show_movies(all_movies)


def create_topography_movie(subj, values, limits, colormap):
    """Read channels and surfaces, and save them as rotating movies

    Parameters
    ----------
    subj : str
        subject code
    values : ndarray
        one value for each channel
    limits : tuple of float, optional
        min and max values for the plots
    colormap : ndarray, optional
        colormap from visvis

    Returns
    -------
    movies : path to file
        file with the movie, it should have an image with the same name but
        ending in _poster.jpg

    """
    chan = get_chan_used_in_analysis(subj)
    if mean(chan.return_xyz()[:, 0]) > 0:
        hemi = 'rh'
    else:
        hemi = 'lh'

    fs_dir = join(REC_DIR, subj, 'mri/proc/freesurfer')
    fs = Freesurfer(fs_dir)
    surf = fs.read_surf(hemi)

    fig = plot_surf(surf)
    plot_chan(chan, fig=fig, values=values, limits=limits, colormap=colormap)
    colorbar()

    movie_name = join(GROUP_DIR, 'powerspectrum', subj + '.mp4')
    make_movie(fig, movie_name, step=3, loop='consistent', rate=10,
               poster=True)
    fig.Destroy()

    return movie_name


def show_movies(movies):
    """Helper function to convert videos to html tag.

    Parameters
    ----------
    movies : path to file
        file with the movie, it should have an image with the same name but
        ending in _poster.jpg

    Returns
    -------
    str
        long string to pass to IPython HTML

    """
    all_tags = []

    for movie_file in movies:

        options = 'controls width="520"'

        video_tag = '<video ' + options
        poster_file = splitext(movie_file)[0] + '_poster.jpg'
        if exists(poster_file):
            with open(poster_file, "rb") as f:
                img = f.read()
            poster_encoded = b64encode(img).decode('ascii')
            video_tag += (' poster="data:image/png;base64,{1}"'
                          ''.format(poster_encoded))
        video_tag += '>'

        with open(movie_file, "rb") as f:
            video = f.read()
        video_encoded = b64encode(video).decode('ascii')

        source_tag = ('<source src="data:video/mp4;base64,{0}" />'
                      ''.format(video_encoded))

        all_tags.append(video_tag + source_tag + '</video>')

    return '\n'.join(all_tags)


def hist_overlap(spindles, width=2, nchan=70):

    s = asarray([x['start_time'] for x in spindles.spindle])
    e = asarray([x['end_time'] for x in spindles.spindle])
    ov = estimate_overlap(s, e)
    x = sum(ov, axis=1)

    v = diff(r_[1, where(x == 1)[0]])

    hist = arange(0, nchan, width)
    return histogram(v, hist)


def hist_overlap2(spindles, width=2, nchan=70):

    s = asarray([x['start_time'] for x in spindles.spindle])
    e = asarray([x['end_time'] for x in spindles.spindle])
    ov = estimate_overlap(s, e)
    x = sum(ov | ov.T, axis=1)

    hist = arange(0, nchan, width)
    h0, h1 = histogram(x, hist)
    bar(h1[:-1], h0, width=width)
