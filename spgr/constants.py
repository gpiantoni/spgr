from logging import getLogger
from os.path import expanduser
from pathlib import Path

PROJECT = 'spgr'
lg = getLogger(PROJECT)

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


HEMI_SUBJ = {'EM09': 'rh',
             'MG17': 'rh',
             'MG33': 'lh',
             'MG37': 'lh',
             'MG61': 'lh',
             'MG63': 'rh',
             }

HP_FILTER = 1
LP_FILTER = 40
RESAMPLE_FREQ = 256

DATA_OPTIONS = {'resample_freq': RESAMPLE_FREQ,
                'hp_filter': 0.5,
                'lp_filter': 50,
                }

CHAN_TYPE = ('grid', 'strip')

DEFAULT_HEMI = 'rh'
SMOOTHING_STD = 10
SMOOTHING_THRESHOLD = 20
FS_AVG = '/home/gio/recordings/EM09/mri/proc/fsaverage'
