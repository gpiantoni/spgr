from glob import glob
from os.path import join
from pickle import load

from phypno.detect import DetectSpindle

# from lsf import map_lsf

from .read_data import DATA_DIR, REC_FOLDER, STAGES

STAGE = 'sleep'
thresh = 300


def calc_spindle_values(subj=None, detection_options=None, ref_to_avg=None):

    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    if ref_to_avg:
        data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]
    else:
        data_file = glob(join(subj_dir, '*_' + STAGE + '_avg.pkl'))[0]

    with open(data_file, 'rb') as f:
        data = load(f)

    detsp = DetectSpindle(**detection_options)

    spindles = detsp(data)
    return spindles



