#!/usr/bin/env python3

from argparse import ArgumentParser
from spgr import read_ecog_recordings
from spgr.constants import OUTPUT_PATH

## ALWAYS GIT COMMIT

parser = ArgumentParser(prog='spgr',
                        description='Full SPGR analysis')
args = parser.parse_args()

if __name__ == '__main__':
    read_ecog_recordings()



## PREPARE PANDOC FILE
# change ### into \n###

# pandoc -s -S --toc test.md -o test.html