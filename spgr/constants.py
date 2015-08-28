from json import load
from os.path import expanduser, join
from pathlib import Path
from os import environ

from phypno.attr import Freesurfer

PROJECT = 'spgr'

HOME = Path(expanduser('~'))

REC_PATH = HOME.joinpath('recordings')
XLTEK_FOLDER = Path('eeg/raw/xltek')
SCORE_FOLDER = Path('doc/scores')
ELEC_FOLDER = Path('doc/elec')
FS_FOLDER = Path('mri/proc/freesurfer')
REC_FOLDER = Path('rec')


PROJECT_PATH = HOME.joinpath('projects').joinpath(PROJECT)
SCRIPTS_PATH = PROJECT_PATH.joinpath('scripts')
DATA_PATH = PROJECT_PATH.joinpath('subjects')
GROUP_PATH = PROJECT_PATH.joinpath('group')
IMAGES_PATH = GROUP_PATH.joinpath('images')
LOG_PATH = GROUP_PATH.joinpath('log')
LOGSRC_PATH = LOG_PATH.joinpath('src')


PARAMETERS_PATH = SCRIPTS_PATH.joinpath(PROJECT).joinpath('parameters.json')
with open(str(PARAMETERS_PATH), 'r') as f:
    PARAMETERS = load(f)
with open(str(PARAMETERS_PATH), 'r') as f:
    PARAMETERS_TXT = f.read()


if not IMAGES_PATH.exists():
    IMAGES_PATH.mkdir(parents=True)
if not LOGSRC_PATH.exists():
    LOGSRC_PATH.mkdir(parents=True)


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


HP_FILTER = 1
LP_FILTER = 40
RESAMPLE_FREQ = 256

DATA_OPTIONS = {'resample_freq': RESAMPLE_FREQ,
                'hp_filter': 0.5,
                'lp_filter': 50,
                }

SPINDLE_OPTIONS = PARAMETERS['SPINDLE_OPTIONS']
SPINDLE_OPTIONS.update(DATA_OPTIONS)

# CHAN_TYPE = ('grid', 'strip')
CHAN_TYPE = ('grid', )

DEFAULT_HEMI = 'rh'
SMOOTHING_STD = 10
SMOOTHING_THRESHOLD = 20
FS_AVG = join(environ['FREESURFER_HOME'], 'subjects', 'fsaverage')

PLOT_COLOR = 'kw'
CHAN_COLOR = 0.8, 0.1, 0.1, 1.
SKIN_COLOR = 0.93, 0.82, 0.81, 0.94
NAN_COLOR = 0.4, 0.4, 0.4, 1.
HIGHLIGHT_COLOR = PARAMETERS['HIGHLIGHT_COLOR']  # color to highlight spindles

ELEVATION = 15
AZIMUTH = 15

RAW_LIMITS_Y = (-200, 200)
HISTOGRAM_WIDTH = 1

fs = Freesurfer(FS_AVG)
surf_avg = getattr(fs.read_brain(), DEFAULT_HEMI)
avg_vert, avg_color, avg_regions = fs.read_label(DEFAULT_HEMI,
                                                 parc_type=PARAMETERS['PARC_TYPE'])
