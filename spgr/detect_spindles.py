from logging import getLogger, INFO, DEBUG
lg = getLogger('spgr')

from copy import deepcopy
from glob import glob
from os.path import join
from pickle import load

from numpy import asarray

from phypno.detect import DetectSpindle

from lsf import map_lsf

from .read_data import DATA_DIR, REC_FOLDER, STAGES

STAGE = 'sleep'

def get_one_chan(data):
    """This could also work but it's very slow:

    def get_one_chan(data):
        for i, chan in enumerate(data.axis['chan'][0]):

            sel_one_chan = Select(chan=(chan, ))
            yield sel_one_chan(data)
    """
    for i_chan in range(data.number_of('chan')[0]):
        one_chan = deepcopy(data)
        for i_trl in range(data.number_of('trial')):
            one_chan.axis['chan'][i_trl] = asarray((data.axis['chan'][i_trl][i_chan], ))
            one_chan.data[i_trl] = data.data[i_trl][(i_chan, ), :]

        yield one_chan


def det_sp_in_one_chan(data):
    spindles = detsp(data)
    return spindles.spindle


def calc_spindle_values(subj=None, detsp=None, ref_to_avg=None):

    assert STAGE in STAGES.keys()

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    if ref_to_avg:
        data_file = glob(join(subj_dir, '*_' + STAGE + '_avg.pkl'))[0]
    else:
        data_file = glob(join(subj_dir, '*_' + STAGE + '.pkl'))[0]

    lg.info('Subj %s, reading data: %s', subj, data_file)
    with open(data_file, 'rb') as f:
        data = load(f)

    lg.info('Subj %s, ready to submit %d jobs', subj,
            data.number_of('chan')[0])
    all_sp = map_lsf(det_sp_in_one_chan, get_one_chan(data),
                     queue='short',
                     variables={'detsp': detsp})
    return [item for sublist in all_sp for item in sublist]
