#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import OrderedDict
from subprocess import check_call

from spgr import read_ecog_recordings
from spgr.constants import LOGSRC_PATH

## ALWAYS GIT COMMIT


all_func = OrderedDict([('-r', 'read_ecog_recordings'),
                         ])

parser = ArgumentParser(prog='spgr',
                        description='SPGR analysis')
for abbr, func_name in all_func.items():
    parser.add_argument(abbr, help=func_name.replace('_', ' '),
                        action='store_true')
args = parser.parse_args()


if __name__ == '__main__':
    for abbr, func_name in all_func.items():
        if getattr(args, abbr[1:]):
            eval(func_name + '()')


## PREPARE PANDOC FILE
check_call(['pandoc', '-s', '-S', '--toc',
            str(LOGSRC_PATH.joinpath('read_ecog_recordings.md')),
            '-o', '/home/gio/projects/spgr/group/log/test.html'])
