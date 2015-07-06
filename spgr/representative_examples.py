from spgr.constants import (SPINDLE_OPTIONS,
                            )
from spgr.detect_spindles import get_spindles

from .log import with_log


@with_log
def Representative_Examples(lg, images_dir):

    lg.info('## Representative Spindles')

    lg.info('### Spindles in each brain region')

    spindles = get_spindles(subj, reref=REREF, **SPINDLE_OPTIONS)
