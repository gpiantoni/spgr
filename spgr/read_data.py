from glob import glob
from logging import getLogger
from os import listdir, makedirs
from os.path import isdir, join, expanduser, basename, splitext, exists
from pickle import dump, load
from re import match

from numpy import where, mean, std, min, max

from phypno import Dataset
from phypno.attr import Scores, Channels
from phypno.trans import Filter, Select, Math, Montage, Resample

lg = getLogger('spgr')

HOME = expanduser('~')

REC_DIR = join(HOME, 'recordings')
xltek_path = 'eeg/raw/xltek'
score_path = 'doc/scores'

STAGES = {'sleep': ('NREM2', 'NREM3'),
          }
DATA_DIR = join(HOME, 'projects/spgr/subjects')

REC_FOLDER = 'rec'
MIN_EPOCHS = 60

HP_FILTER = 1
LP_FILTER = 40

RESAMPLE_FREQ = 100

thresh = 300


def is_grid(label):
    """Determine if a channel is a grid channel based on simple heuristics."""
    G1 = match('.*G[0-9]{1,2}$', label)
    CING1 = match('.*CING[0-9]$', label)
    EMG1 = match('.*EMG[0-9]$', label)
    GR1 = match('.*GR[0-9]{1,2}$', label.upper())
    return (G1 and not CING1 and not EMG1) or GR1


def save_wake_sleep_data(xltek_file, subj, epochs):

    d = Dataset(xltek_file)
    gr_chan = [x for x in d.header['chan_name'] if is_grid(x)]

    subj_dir = join(DATA_DIR, subj, REC_FOLDER)
    if not isdir(subj_dir):
        makedirs(subj_dir)

    for stage, epochs_in_stage in epochs.items():
        start_time = [x['start_time'] for x in epochs_in_stage]
        end_time = [x['end_time'] for x in epochs_in_stage]
        data = d.read_data(begtime=start_time, endtime=end_time, chan=gr_chan)

        hp_filt = Filter(low_cut=HP_FILTER, s_freq=data.s_freq)
        lp_filt = Filter(high_cut=LP_FILTER, s_freq=data.s_freq)
        data = lp_filt(hp_filt(data))

        res = Resample(s_freq=RESAMPLE_FREQ)
        data = res(data)

        # remove bad channels
        calc_std = Math(operator_name='std', axis='time')
        std_per_chan = calc_std(data)

        chan_val = std_per_chan(trial=0)
        lg.debug('Channels values: m %f, std %f, min %f, max %f',
                 mean(chan_val), std(chan_val), min(chan_val), max(chan_val))

        good_chan_idx = where((chan_val > .001) & (chan_val < thresh))[0]
        good_chan = data.axis['chan'][0][good_chan_idx]
        normal_chan = Select(chan=good_chan)

        # save channels used in the analysis
        pkl_file = join(subj_dir, splitext(basename(xltek_file))[0] + '_' +
                        stage + '_chan' + '.pkl')
        with open(pkl_file, 'wb') as f:
            dump(good_chan, f)

        data = normal_chan(data)

        for reref in ('', '_avg'):
            if reref == '_avg':
                to_avg = Montage(ref_to_avg=True)
                data = to_avg(data)

            pkl_file = join(subj_dir, splitext(basename(xltek_file))[0] + '_' +
                            stage + reref + '.pkl')
            with open(pkl_file, 'wb') as f:
                dump(data, f)


def read_score_per_subj(subj, save_data=False):

    score_dir = join(REC_DIR, subj, score_path)
    all_xml = listdir(score_dir)

    good_xltek = None

    if all_xml:
        first_dataset = True  # save only one dataset per subject

        for one_xml in all_xml:
            score = Scores(join(score_dir, one_xml))

            lg.info(one_xml)
            epochs = {}
            enough_epochs = True
            for stage_name, stages in STAGES.items():
                epochs_in_stage = score.get_epochs(stages)
                lg.info('    %s has % 5.1f min', stage_name,
                        len(epochs_in_stage) / 2.)

                epochs[stage_name] = epochs_in_stage[:MIN_EPOCHS]  # REMOVE THIS
                if len(epochs[stage_name]) < MIN_EPOCHS:
                    enough_epochs = False

            if enough_epochs and first_dataset:
                xltek_file = join(REC_DIR, subj, xltek_path, one_xml[:-11])
                if save_data:
                    lg.info('Enough epochs, saving data')
                    save_wake_sleep_data(xltek_file, subj, epochs)
                first_dataset = False
                good_xltek = xltek_file

    return good_xltek


def get_chan_used_in_analysis(all_subj):
    """If we use multiple datasets, we need to take into account, also above
    at: dump(good_chan, f)

    """
    all_chan = []

    for subj in all_subj:
        chan_file = join(REC_DIR, subj, 'doc/elec/elec_pos_adjusted_renamed.csv')
        if not exists(chan_file):
            chan_file = join(REC_DIR, subj, 'doc/elec/elec_pos_adjusted.csv')
        chan = Channels(chan_file)

        subj_dir = join(DATA_DIR, subj, REC_FOLDER)
        gr_chan_file = glob(join(subj_dir, '*_' + 'chan' + '.pkl'))[0]

        with open(gr_chan_file, 'rb') as f:
            gr_chan = load(f)

        chosen_chan = chan(lambda x: x.label in gr_chan)
        lg.info('%s analysis chan %d, with location %d',
                subj, len(gr_chan), chosen_chan.n_chan)
        all_chan.append(chosen_chan)

    return all_chan
