from spgr.constants import (CHAN_TYPE,
                            DATA_OPTIONS,
                            HEMI_SUBJ,
                            MIN_DURATION,
                            PERIOD,
                            STAGES,
                            SUBJECTS,
                            )
from spgr.read_data import save_data, select_scores

from .log import with_log


@with_log
def Read_ECoG_Recordings(lg, img_dir):

    lg.info('## Selection of Patients')

    scores = select_scores(STAGES, MIN_DURATION,
                           sorted(SUBJECTS['all'], reverse=True))

    lg.info('## Patients with sufficient recordings:')
    lg.info('\t'.join(sorted(scores.keys())))

    good_subj = sorted(list(scores.keys()))

    lg.info('Grid patients with scalp : ' +
            ', '.join(sorted(set(good_subj) & set(SUBJECTS['grid_scalp']))))
    lg.info('Grid patients without scalp : ' +
            ', '.join(sorted(set(good_subj) & set(SUBJECTS['grid_noscalp']))))
    lg.info('Depth patients without scalp : ' +
            ', '.join(sorted(set(good_subj) -
                             set(SUBJECTS['grid_scalp']) -
                             set(SUBJECTS['grid_noscalp']))))

    lg.info('## Read Recordings')
    for subj in HEMI_SUBJ:
        save_data(subj, scores[subj], PERIOD, STAGES, chan_type=CHAN_TYPE,
                  **DATA_OPTIONS)
