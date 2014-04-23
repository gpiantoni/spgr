from glob import glob
from os.path import join
from pickle import load, dump

from phypno.detect import DetectSpindle
from lsf import map_lsf

from .read_data import DATA_DIR, REC_FOLDER, STAGES


STAGE = 'sleep'


def det_sp_in_one_epoch(one_trial_data=None):
    from phypno.trans import Select, Math, Montage

    calc_std = Math(operator_name='std', axis='time')
    std_per_chan = calc_std(one_trial_data)
    good_chan = where((std_per_chan(trial=0) > .001) & (std_per_chan(trial=0) < thresh))[0]
    normal_chan = Select(chan=one_trial_data.axis['chan'][0][good_chan])

    if ref_to_avg:
        reref = Montage(ref_to_avg=True)
        one_trial_data = reref(one_trial_data)

    one_trial_data = normal_chan(one_trial_data)

    spindles = detsp(one_trial_data)
    return spindles.spindle


def calc_spindle_values(subj=None, detection_options=None, ref_to_avg=None):

    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]  # multiple files

    with open(data_file, 'rb') as f:
        data = load(f)

    detsp = DetectSpindle(**detection_options)

    all_sp = map_lsf(det_sp_in_one_epoch, iter(data),
                     queue='short',
                     variables={'detsp': detsp,
                                'thresh': 2,
                                'ref_to_avg': ref_to_avg},
                     imports={'numpy': 'where'})

    spindles = [x for sp in all_sp for x in sp]
    return spindles


"""
multiprocessing works sometimes, not clear when it does and when it doesn't.
Maybe in a new kernel, you need to be really clean about starting subprocesses
and closing them.

all_spindles = list(map(calc_spindle_values, all_subj))
"""
