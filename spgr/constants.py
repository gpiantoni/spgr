from json import load
from os.path import expanduser
from pathlib import Path

from phypno.attr import Freesurfer

PROJECT = 'spgr'

HOME = Path(expanduser('~'))

REC_PATH = HOME.joinpath('recordings')
XLTEK_FOLDER = Path('eeg/raw/xltek')
ELEC_FOLDER = Path('doc/elec')
FS_FOLDER = Path('mri/proc/freesurfer')

PROJECT_PATH = HOME.joinpath('projects').joinpath(PROJECT)
SCRIPTS_PATH = PROJECT_PATH.joinpath('scripts')
DATA_PATH = PROJECT_PATH.joinpath('subjects')
GROUP_PATH = PROJECT_PATH.joinpath('group')
IMAGES_PATH = GROUP_PATH.joinpath('images')
LOG_PATH = GROUP_PATH.joinpath('log')
LOGSRC_PATH = LOG_PATH.joinpath('src')
SCORES_PATH = GROUP_PATH.joinpath('scores')

PARAMETERS_PATH = SCRIPTS_PATH.joinpath(PROJECT).joinpath('parameters.json')
with open(str(PARAMETERS_PATH), 'r') as f:
    PARAMETERS = load(f)
with open(str(PARAMETERS_PATH), 'r') as f:
    PARAMETERS_TXT = f.read()


if not IMAGES_PATH.exists():
    IMAGES_PATH.mkdir(parents=True)
if not LOGSRC_PATH.exists():
    LOGSRC_PATH.mkdir(parents=True)

# READ DATA-------------------------------------------------------------------#
REC_FOLDER = Path('rec')

PERIOD = 'sleep'
STAGES = ('NREM2', )
MIN_DURATION = 60 * 60

HEMI_SUBJ = {'EM09': 'rh',
             'MG17': 'rh',
             'MG33': 'lh',
             'MG37': 'lh',
             'MG61': 'lh',
             'MG63': 'rh',
             'MG64': 'rh',
             'MG67': 'rh',
             }


# define patients with grids
SUBJECTS = {'grid_scalp': ('EM03',
                           'EM09',
                           'MG07',
                           'MG11',
                           'MG12',
                           'MG15',
                           'MG17',
                           'MG29',
                           'MG31',
                           'MG32',
                           'MG33',
                           'MG37',
                           'MG42',
                           'MG45',
                           'MG48',
                           'MG49',
                           'MG51',
                           'MG56',
                           'MG61',
                           'MG63',
                           'MG65',
                           'MG72',
                           'MG83'),
            'grid_noscalp': ('MG28',
                             'MG64',
                             'MG67',
                             'MG75',
                             'MG81'),
            'many_SZ': ('MG49',
                        'MG74',
                        'MG75',
                        'MG81',
                        'MG83'),
            'not_in_rec': ('MG12', ),
            }


SUBJECTS['all'] = ((set(SUBJECTS['grid_scalp']) |
                    set(SUBJECTS['grid_noscalp'])) -
                   set(SUBJECTS['many_SZ']) -
                   set(SUBJECTS['not_in_rec']))


def SESSION(subj):
    if subj in ('MG64', ):
        return 'B'
    else:
        return 'A'

DATA_OPTIONS = {'resample_freq': 256,
                'hp_filter': 0.5,
                'lp_filter': 50,
                }

# CHAN_TYPE = ('grid', 'strip')
CHAN_TYPE = ('grid', )
ALL_REREF = ('avg', )

# SPINDLE OPTIONS-------------------------------------------------------------#
SPINDLE_FOLDER = Path('spindles')
SPINDLE_OPTIONS = PARAMETERS['SPINDLE_OPTIONS']
SPINDLE_OPTIONS.update(DATA_OPTIONS)

# SURFACE OPTIONS-------------------------------------------------------------#
APARC_FOLDER = Path('aparc')
PROJ_FOLDER = Path('proj')

DEFAULT_HEMI = 'rh'
SMOOTHING_STD = 10
SMOOTHING_THRESHOLD = 20
FS_AVG = GROUP_PATH / 'fsaverage'
MORPH_SMOOTHING = None

fs = Freesurfer(str(FS_AVG))
avg_surf = getattr(fs.read_brain(), DEFAULT_HEMI)
avg_vert, _, avg_regions = fs.read_label(DEFAULT_HEMI,
                                         parc_type=PARAMETERS['PARC_TYPE'])

P_THRESHOLD = 0.05
P_CORRECTION = 'fdr'

# PLOT OPTIONS----------------------------------------------------------------#
DPI = PARAMETERS['DPI']
RAW_LIMITS_Y = (-100, 100)
HIST_WIDTH = 1
HIST_BAR_COLOR = 'w'
HIST_BAR_WIDTH = 4
HIST_N_CHAN = {'avg': 20,
               15: 40}

SINGLE_CHAN_LIMITS = {'density': (0, 1.5),
                      'duration': (1, 1.5),
                      'peak_val': (0, 70),
                      'peak_freq': (11, 13),
                      }

COOCCUR_CHAN_LIMITS = {'avg': (4, 6),
                       15: (5, 10)
                       }

DIR_MAT_RATIO = 2 / 1
DIR_SURF_RATIO = 4 / 3
DIR_SUMMARY_RATIO = 2 / 1
DIR_SUMMARY_MINCNT = 100

COLORMAP = PARAMETERS['COLORMAP']
SATURATION_LEVEL = 0.4
CHAN_COLOR = 0.8, 0.1, 0.1, 1.
SKIN_COLOR = 0.93, 0.82, 0.81, 0.94
NAN_COLOR = 0.2, 0.2, 0.2, 1.
IMAGE_NAN_COLOR = 1., 1., 1., 1.
HIGHLIGHT_COLOR = PARAMETERS['HIGHLIGHT_COLOR']  # color to highlight spindles
HIGHLIGHT_HEIGHT = 200
SURF_PLOT_SIZE = 50, 35
REPR_PLOT_SIZE = 10, 8
SINGLE_SUBJ_SURF = 20, 15
AVERAGE_BW_SURF = 40, 30  # the actual dimensions are 80, 60 but it's too big when saving with read_pixels
HIST_FIG_SIZE = 40, 30  # actual size in mm for figure


TICKS_FONT_SIZE = 8
LABEL_FONT_SIZE = 10

CENTER = (25, -15, 18)
ELEVATION = -6
AZIMUTH = 96
