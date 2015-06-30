from os.path import expanduser, join
from pathlib import Path
from os import environ

PROJECT = 'spgr'

HOME = Path(expanduser('~'))

REC_PATH = HOME.joinpath('recordings')
XLTEK_FOLDER = Path('eeg/raw/xltek')
SCORE_FOLDER = Path('doc/scores')
ELEC_FOLDER = Path('doc/elec')
FS_FOLDER = Path('mri/proc/freesurfer')
REC_FOLDER = Path('rec')

DATA_PATH = HOME.joinpath('projects').joinpath(PROJECT).joinpath('subjects')
GROUP_PATH = HOME.joinpath('projects').joinpath(PROJECT).joinpath('group')
IMAGES_PATH = GROUP_PATH.joinpath('images')
LOG_PATH = GROUP_PATH.joinpath('log')
OUTPUT_PATH = LOG_PATH.joinpath('src')


if not IMAGES_PATH.exists():
    IMAGES_PATH.mkdir(parents=True)
if not OUTPUT_PATH.exists():
    OUTPUT_PATH.mkdir(parents=True)


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

SPINDLE_OPTIONS = {'method': 'Nir2011',
                   'frequency': (9, 16),
                   'duration': (0.5, 2),
                   }
SPINDLE_OPTIONS.update(DATA_OPTIONS)

# CHAN_TYPE = ('grid', 'strip')
CHAN_TYPE = ('grid', )

DEFAULT_HEMI = 'rh'
SMOOTHING_STD = 10
SMOOTHING_THRESHOLD = 20
FS_AVG = join(environ['FREESURFER_HOME'], 'subjects', 'fsaverage')
