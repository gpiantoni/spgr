#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import OrderedDict
from subprocess import check_call

from PyQt4.QtGui import QApplication

from spgr import (Read_ECoG_Recordings,
                  Electrode_Locations,
                  Single_Channel_Statistics,
                  Cooccurrence_of_Spindles,
                  )
#                  Spindle_Detection_Method,

from spgr.constants import LOG_PATH, LOGSRC_PATH


# ALWAYS GIT COMMIT
all_func = OrderedDict([('-r', 'Read_ECoG_Recordings'),
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


app = QApplication([])


md_files = []
if __name__ == '__main__':
    for abbr, func_name in all_func.items():
        if getattr(args, abbr[1:]):
            eval(func_name + '()')
            md_files.append(str(LOGSRC_PATH.joinpath(func_name + '.md')))

# app.exec_()

# PREPARE PANDOC FILE
cmd = ['pandoc', '-s', '-S', '--toc']
cmd.extend(md_files)
cmd.extend(['-o', str(LOG_PATH.joinpath('main.' + args.to))])
check_call(cmd)
