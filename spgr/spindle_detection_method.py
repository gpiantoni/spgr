from .constants import (CHAN_TYPE,
                        HEMI_SUBJ,
                        SPINDLE_OPTIONS,
                        )
from .detect_spindles import get_spindles

from .log import with_log


@with_log
def Spindle_Detection_Method(lg, images_dir):

    lg.info('## Method description')
    lg.info('Spindles were detected on each electrode independently, using '
            'previously reported methods with stringent criteria ({method}). '
            'The raw signal was filtered between {freq0} and {freq1} Hz. '
            'Spindle duration had to be between {dur0} and {dur1} s.'
            ''.format(method=SPINDLE_OPTIONS['method'],
                      freq0=SPINDLE_OPTIONS['frequency'][0],
                      freq1=SPINDLE_OPTIONS['frequency'][1],
                      dur0=SPINDLE_OPTIONS['duration'][0],
                      dur1=SPINDLE_OPTIONS['duration'][1]))

    for ref in ('avg', 15):
        for subj in HEMI_SUBJ:
            print(subj)
            sp = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref,
                              **SPINDLE_OPTIONS)
