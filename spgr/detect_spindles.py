from glob import glob
from os import mkdir
from os.path import join, exists
from pickle import load, dump

from lsf import map_lsf

from phypno.detect import DetectSpindle
from phypno.trans import Select, Math

from .read_data import DATA_DIR, REC_FOLDER, STAGES


STAGE = 'sleep'


def det_sp_in_one_epoch(one_trial_data=None, detection_options=None):

    calc_std = Math(operator_name='std', axis='time')
    std_per_chan = calc_std(one_trial_data)
    good_chan = where((std_per_chan(trial=0) > .001) & (std_per_chan(trial=0) < thresh))[0]
    normal_chan = Select(chan=one_trial_data.axis['chan'][0][good_chan])

    detsp = DetectSpindle(**detection_options)

    spindles = detsp(one_trial_data)
    return spindles.spindle


def calc_spindle_values(subj):

    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]  # possible multiple files

    with open(data_file, 'rb') as f:
        data = load(f)

    all_sp = map_lsf(det_sp_in_one_epoch, iter(data),
                     queue='short',
                     variables={'detsp': detsp, 'thresh': 2},
                     imports={'phypno.trans': ('Math', 'Select'), 'numpy': 'where'})

    spindles = [x for sp in all_sp for x in sp]
    return spindles


"""
multiprocessing works sometimes, not clear when it does and when it doesn't.
Maybe in a new kernel, you need to be really clean about starting subprocesses
and closing them.

all_spindles = list(map(calc_spindle_values, all_subj))
"""
