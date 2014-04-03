#!/usr/bin/env python3

from getpass import getuser
from json import dump
from numpy import asarray, min, max, mean, median
from os import listdir
from os.path import join
from pickle import load
from sys import argv, path

if getuser() == 'gio':
    PHYPNO_PATH = '/home/gio/tools/phypno'
    log_dir = '/home/gio/projects/spgr/group/log_detect_spindle'
    data_dir = '/home/gio/projects/spgr/subjects'
else:
    PHYPNO_PATH = '/PHShome/gp902/toolbox/python/usr/phypno'
    log_dir = '/PHShome/gp902/projects/spgr/group/log_detect_spindle'
    data_dir = '/PHShome/gp902/projects/spgr/subjects'


log_file = join(log_dir, 'log_' + '_'.join(argv[1:]) + '.json')
path.append(PHYPNO_PATH)
from phypno.detect import DetectSpindle


def describe_spindles(spindle):
    """function to summarize spindles

    Parameters
    ----------
    spindle : list of dict
        output of DetectSpindle

    Returns
    -------
    dict
        dictionary to export to json.

    """
    output = {}
    output['number'] = len(spindle)
    all_dur = asarray([sp['end_time'] - sp['start_time'] for sp in spindle])

    output['duration'] = {}
    if len(spindle) > 0:
        output['duration']['min'] = min(all_dur)
        output['duration']['mean'] = mean(all_dur)
        output['duration']['median'] = median(all_dur)
        output['duration']['max'] = max(all_dur)

    output['chan'] = list(set([sp['chan'] for sp in spindle]))

    return output


data_types = ('wake_data.pkl', 'sleep_data.pkl')

OPTIONS = {'frequency': (float(argv[2]), float(argv[3])),
           'method': argv[1],

           'threshold': argv[4],
           'threshold_options': {'peak_width': 1,
                                 'select_width': 1,
                                 },
           'criteria': {'duration': (float(argv[7]), float(argv[8])),
                        },
           }


if argv[1] == 'hilbert':
    OPTIONS['method_options'] = {}

elif argv[1] == 'wavelet':
    OPTIONS['method_options'] = {'detection_wavelet': {'M_in_s': 1,
                                                       'w': 7,
                                                       },
                                 'detection_smoothing': {'window': 'boxcar',
                                                         'length': 1,
                                                         },
                                 }

if argv[4] in ('relative', 'absolute'):
    OPTIONS['threshold_options'] = {'detection_value': float(argv[5]),
                                    'selection_value': float(argv[6]),
                                     }
elif argv[4] in ('maxima', ):

    OPTIONS['threshold_options'] = {'peak_width': float(argv[5]),
                                    'select_width': float(argv[6]),
                                    }

if argv[9] == '1':
    OPTIONS['criteria']['peak_in_fft'] = {'length': 1}

ALL_SUBJ = []

for subj in listdir(data_dir):

    subj_dir = join(data_dir, subj)

    sp_one_subj = {}
    for data_type in data_types:
        with open(join(subj_dir, data_type), 'rb') as data_file:
            data = load(data_file)

        det_sp = DetectSpindle(**OPTIONS)
        spindles = det_sp(data)

        cond_name = data_type.split('_')[0]
        sp_one_subj[cond_name] = describe_spindles(spindles.spindle)

    sp_one_subj['subj'] = subj
    ALL_SUBJ.append(sp_one_subj)


to_json = {'parameters': OPTIONS,
           'spindles': ALL_SUBJ}

with open(log_file, 'w') as j:
    dump(to_json, j, indent=2)
