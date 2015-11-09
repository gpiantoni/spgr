from functools import partial
from multiprocessing import Pool
from numpy import (array,
                   c_,
                   exp,
                   isfinite,
                   fill_diagonal,
                   flipud,
                   log,
                   min,
                   nanmean,
                   NaN,
                   r_,
                   seterr,
                   sum,
                   where,
                   zeros)
from vispy.plot import Fig
from vispy.scene.visuals import Image
from vispy.io import write_png

from .constants import (ALL_REREF,
                        COLORMAP,
                        DIR_MAT_RATIO,
                        DIR_SUMMARY_RATIO,
                        DIR_SUMMARY_MINCNT,
                        DIR_SURF_RATIO,
                        HEMI_SUBJ,
                        IMAGE_NAN_COLOR,
                        P_CORRECTION,
                        P_THRESHOLD,
                        SPINDLE_OPTIONS,
                        SURF_PLOT_SIZE)
from .detect_spindles import get_spindles
from .plot_spindles import plot_lmer
from .spindle_source import get_chan_with_regions, get_regions_with_elec

from .log import with_log


from numpy import triu
from numpy.random import binomial, seed
from rpy2.robjects.numpy2ri import activate
from rpy2.robjects import r

activate()  # conversion R <-> numpy
p_adjust = r['p.adjust']


NULL_PROBABILITY = .5
N_RND = 10000

REGIONS = [
 'medialorbitofrontal_1',
 'medialorbitofrontal_2',
 'lateralorbitofrontal_1',
 'lateralorbitofrontal_2',
 'parsorbitalis_1',  # empty
 'rostralmiddlefrontal_1',
 'rostralmiddlefrontal_2',
 'rostralmiddlefrontal_3',
 'caudalmiddlefrontal_1',
 'parsopercularis_1',
 'parstriangularis_1',
 'superiorfrontal_4',
 'precentral_1',
 'precentral_2',
 'precentral_3',
 'precentral_4',
 'postcentral_1',
 'postcentral_2',
 'postcentral_3',
 'superiorparietal_1',
 'superiorparietal_2',
 'supramarginal_1',
 'supramarginal_2',
 'inferiorparietal_1',
 'inferiorparietal_2',
 'inferiorparietal_3',
 'insula_2',  # empty
 'temporalpole_1',  # empty
 'superiortemporal_1',
 'superiortemporal_2',
 'middletemporal_2',
 'middletemporal_1',
 'inferiortemporal_1',
 'inferiortemporal_2',
 'bankssts_1',
 'fusiform_1',
 'lateraloccipital_2',
 'lateraloccipital_3',
 ]


@with_log
def Direction_of_Spindles(lg, images_dir):

    lg.info('## Cooccurrence_of_Spindles')

    for reref in ALL_REREF:

        lg.info('### reref {}'.format(reref))

        if reref == 'avg':
            regions = REGIONS
        else:
            regions = get_regions_with_elec(reref)
        x = zeros((len(regions), len(regions)))
        for subj in HEMI_SUBJ:

            spindles = get_spindles(subj, reref=reref, **SPINDLE_OPTIONS)
            chan = get_chan_with_regions(subj, reref=reref)
            chan_with_regions = dict(zip(spindles.chan_name,
                                         chan.return_attr('region',
                                                          spindles.chan_name)))

            for sp0 in spindles.spindle:  # this is the lead
                for sp1 in spindles.spindle:  # this is the follower

                    if ((sp1['start_time'] > sp0['start_time']) and
                        (sp1['start_time'] < sp0['end_time'])):

                        region0 = chan_with_regions[sp0['chan']]
                        region1 = chan_with_regions[sp1['chan']]
                        try:
                            i0 = regions.index(region0[7:])
                            i1 = regions.index(region1[7:])
                            x[i0, i1] += 1
                        except ValueError:  # if it's not in the regions of interest
                            pass

        old_warnings = seterr(all="ignore")

        _direction_summary(lg, x)

        img = _make_direction_matrix(x)

        png_name = 'direction_image_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        write_png(png_file, img)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))

        d = _calc_dir_summary(x)

        coef = dict(zip(regions, d))

        _partial_null = partial(_compute_null_direction, x)

        with Pool() as p:
            n_d = p.map(_partial_null, range(N_RND))
        n_d = r_[n_d]

        lg.info('\nCorrected at {} {}'.format(P_CORRECTION, P_THRESHOLD))

        uncorr_pv = min(c_[sum(d >= n_d, axis=0) / N_RND,
                           sum(d <= n_d, axis=0) / N_RND], axis=1)
        pv = array(p_adjust(uncorr_pv * 2,  # two-tailed
                            method=P_CORRECTION))

        pvalues = dict(zip(regions, pv))

        for region, _ in sorted(coef.items(), key=lambda x: x[1]):

            if pvalues[region] < P_THRESHOLD:

                if coef[region] > 0:
                    coef_print = '{:5.3f}:1'.format(exp(coef[region]))
                else:
                    coef_print = '1:{:5.3f}'.format(exp(-coef[region]))

                lg.info('{:30} coef = {:7},  p-value = {:.3f}'
                        ''.format(region,
                                  coef_print,
                                  pvalues[region]))
        limits = -log(DIR_SURF_RATIO), log(DIR_SURF_RATIO)
        v = plot_lmer(coef, limits=limits, size_mm=SURF_PLOT_SIZE)
        png_name = 'direction_map_{}.png'.format(reref)
        png_file = str(images_dir.joinpath(png_name))
        v.save(png_file)
        lg.info('![{}]({})'.format('{}'.format(reref),
                png_file))

        seterr(**old_warnings)


def _make_direction_matrix(x):

    from numpy import isnan
    from phypno.viz.base import normalize
    from vispy.color import get_colormap

    cm = get_colormap(COLORMAP)

    c = log(x / x.T)
    c[~isfinite(c)] = NaN
    fill_diagonal(c, NaN)

    # colormap only accepts column vector
    val = flipud(c).reshape((-1, ))
    val = normalize(val, -log(DIR_MAT_RATIO), log(DIR_MAT_RATIO))
    clr = cm[val].rgba

    # convert empty values to NaN color
    hasnan = isnan(clr).all(axis=1)
    clr[hasnan, :] = IMAGE_NAN_COLOR

    # convert back to square matrix + color dimension
    clr = clr.reshape((c.shape[0], c.shape[1], 4))

    f = Fig()
    plt = f[0, 0]
    im = Image(clr)
    plt.view.add(im)

    plt.view.camera = 'panzoom'
    plt.view.camera.aspect = 1
    plt.view.camera.set_range()

    region = (160, 34, 480, 479)
    img = f.render(region)

    """convert img to transparent, if IMAGE_NAN_COLOR has unique values
    img_r = img.reshape(-1, 4)
    nan_vxl = (img_r == (array(IMAGE_NAN_COLOR) * 255).astype('uint8')).all(axis=1)
    img_r[nan_vxl, :] = array([0, 0, 0, 0])
    img = img_r.reshape((img.shape[0], img.shape[1], 4))
    """

    return img


def _compute_null_direction(x, seed_value):
    s_x = triu(x + x.T)

    seed(seed_value)
    x2 = zeros(x.shape)

    for i0 in range(x.shape[0]):
        for i1 in range(x.shape[1]):
            if i0 <= i1:
                b = binomial(s_x[i0, i1], NULL_PROBABILITY)
                x2[i0, i1] = b
                x2[i1, i0] = s_x[i0, i1] - b

    return _calc_dir_summary(x2)


def _calc_dir_summary(x):
    # d = sum(x, axis=1) / (sum(x, axis=0) + sum(x, axis=1)) * 100

    c = log(x / x.T)
    c[~ isfinite(c)] = NaN
    d = nanmean(c, axis=1)

    return d


def _direction_summary(lg, x):
    """Create a table with the summary of the direction results"""

    c = log(x / x.T)
    c[~isfinite(c)] = NaN
    fill_diagonal(c, NaN)
    cnt = x + x.T

    lg.info(' {:<30} {:<30} {:<15} {:<12}'
            ''.format('From', 'To', '# Spindle Pairs', 'Ratio'))
    lg.info('-' * 30 + ' ' + '-' * 30 + ' ' + '-' * 17 + ' ' + '-' * 12 + ' ')

    pairs = []

    for i0, i1 in zip(*where(c > log(DIR_SUMMARY_RATIO))):
        if cnt[i0, i1] >= DIR_SUMMARY_MINCNT:
            one_pair = {'from': _rename_region(REGIONS[i0]),
                        'to': _rename_region(REGIONS[i1]),
                        'cnt': int(cnt[i0, i1]),
                        'ratio': exp(c[i0, i1])}
            pairs.append(one_pair)

    for one_pair in sorted(pairs, key=lambda k: k['ratio'], reverse=True):
        lg.info('{from:<30} {to:<30}{cnt: 17d} {ratio: 10.2f}:1 '
                ''.format(**one_pair))

    lg.info('\n')  # otherwise no figure in html


def _rename_region(s):
    s = ' ('.join(s.split('_')) + ')'
    s = s.replace('medial', 'medial ')
    s = s.replace('lateral', 'lateral ')
    s = s.replace('superior', 'superior ')
    s = s.replace('middle', 'middle ')
    s = s.replace('inferior', 'inferior ')
    s = s.replace('rostral', 'rostral ')
    s = s.replace('caudal', 'caudal ')
    s = s.replace('orbito', 'orbital ')
    s = s.replace('pars', 'pars ')
    s = s.replace('pole', ' pole')
    return s
