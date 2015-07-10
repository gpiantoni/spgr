#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import OrderedDict
from datetime import datetime
from subprocess import check_call

from PyQt4.QtGui import QApplication

from spgr import (Read_ECoG_Recordings,
                  Representative_Examples,
                  Electrode_Locations,
                  Single_Channel_Statistics,
                  Cooccurrence_of_Spindles,
                  )
#                  Spindle_Detection_Method,

from spgr.constants import LOG_PATH, LOGSRC_PATH, PROJECT


# ALWAYS GIT COMMIT
all_func = OrderedDict([('-r', 'Read_ECoG_Recordings'),
                        ('-x', 'Representative_Examples'),
                        ('-e', 'Electrode_Locations'),
                        ('-d', 'Spindle_Detection_Method'),  # does not work
                        ('-s', 'Single_Channel_Statistics'),
                        ('-c', 'Cooccurrence_of_Spindles'),
                        ])

parser = ArgumentParser(prog='spgr',
                        description='SPGR analysis')
parser.add_argument('--to', default='html',
                    help='format to export to: html (default), pdf')

for abbr, func_name in all_func.items():
    parser.add_argument(abbr, help=func_name.replace('_', ' '),
                        action='store_true')
args = parser.parse_args()


if __name__ == '__main__':

    app = QApplication([])

    for abbr, func_name in all_func.items():
        if getattr(args, abbr[1:]):
            eval(func_name + '()')

    # PREPARE PANDOC FILE
    md_files = []
    for func_name in all_func.values():
        one_md_file = LOGSRC_PATH.joinpath(func_name + '.md')
        if one_md_file.exists():
            md_files.append(str(one_md_file))

    LOGOUTPUT_PATH = LOG_PATH.joinpath(args.to)
    if not LOGOUTPUT_PATH.exists():
        LOGOUTPUT_PATH.mkdir()
    t = datetime.now()
    output_name = t.strftime(PROJECT + '_%y%m%d_%H%M%S') + '.' + args.to

    cmd = ['pandoc', '-s', '-S', '--toc']
    cmd.extend(md_files)
    cmd.extend(['-o', str(LOGOUTPUT_PATH.joinpath(output_name))])
    check_call(cmd)
