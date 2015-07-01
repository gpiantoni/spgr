from spgr.constants import PERIOD, STAGES, MIN_DURATION, HEMI_SUBJ, DATA_OPTIONS, CHAN_TYPE
from spgr.read_data import save_data, select_scores

from .log import with_log


@with_log
def Read_ECoG_Recordings(lg):

    lg.info('## Selection of Patients')

    # define patients with grids
    subj_with_grid_scalp = ('EM03', 'EM09', 'MG07', 'MG11', 'MG12', 'MG15', 'MG17',
                            'MG29', 'MG31', 'MG32', 'MG33', 'MG37', 'MG42', 'MG45',
                            'MG48', 'MG49', 'MG51', 'MG56', 'MG61', 'MG63', 'MG65', 'MG72', 'MG83')
    subj_with_grid_only = ('MG28', 'MG64', 'MG67', 'MG75', 'MG81')

    # define patients to exclude
    subj_with_many_SZ = ('MG49', 'MG74', 'MG75', 'MG81', 'MG83')
    not_in_rec = ('MG12', )

    all_subj = ((set(subj_with_grid_scalp) | set(subj_with_grid_only)) -
                set(subj_with_many_SZ) - set(not_in_rec))

    scores = select_scores(STAGES, MIN_DURATION,
                           sorted(all_subj, reverse=True))

    lg.info('## Patients with sufficient recordings:')
    lg.info('\t'.join(sorted(scores.keys())))

    good_subj = sorted(list(scores.keys()))

    lg.info('Grid patients with scalp : ' +
            ', '.join(sorted(set(good_subj) & set(subj_with_grid_scalp))))
    lg.info('Grid patients without scalp : ' +
            ', '.join(sorted(set(good_subj) & set(subj_with_grid_only))))
    lg.info('Depth patients without scalp : ' +
            ', '.join(sorted(set(good_subj) - set(subj_with_grid_scalp) -
                             set(subj_with_grid_only))))

    lg.info('## Read Recordings')

    for subj in HEMI_SUBJ:
        save_data(subj, scores[subj], PERIOD, STAGES, chan_type=CHAN_TYPE,
                  **DATA_OPTIONS)
