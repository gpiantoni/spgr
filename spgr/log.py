from datetime import datetime
from logging import getLogger, INFO, FileHandler, StreamHandler, Formatter

from spgr.constants import OUTPUT_PATH, PROJECT

from socket import gethostname
from subprocess import check_output
import phypno
import spgr


def git_hash(package):
    package_path = package.__path__[0] + '/../.git'
    version = check_output('git --git-dir ' + package_path +
                           ' rev-parse HEAD', shell=True).decode('utf-8')
    return version

def git_info(lg):

    lg.info('{} {} '.format(PROJECT, git_hash(spgr)))
    lg.info('{} {} '.format('phypno', git_hash(phypno)))
    lg.info(check_output('pip freeze', shell=True).decode('utf-8'))


def with_log(function):

    def add_log():
        log_file = OUTPUT_PATH.joinpath(function.__name__ + '.md')
        if log_file.exists():
            log_file.unlink()

        lg = getLogger(PROJECT)
        lg.setLevel(INFO)

        formatter = Formatter('%(message)s  ')

        lg.handlers = []

        fh = FileHandler(str(log_file))
        fh.setFormatter(formatter)
        lg.addHandler(fh)

        ch = StreamHandler()
        ch.setFormatter(formatter)
        lg.addHandler(ch)

        lg.info('# {}'.format(function.__name__.replace('_', ' ')))

        git_info(lg)

        lg.info('## Started')
        t0 = datetime.now()
        lg.info('{} on {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S'),
                                  gethostname()))

        function(lg)

        lg.info('## Finished')
        t1 = datetime.now()
        lg.info('{} after {}'.format(t1.strftime('%Y-%m-%d %H:%M:%S'),
                                     str(t1 - t0)[:-7]))

        with open(str(log_file), 'r') as f:
            s = f.read()
        s = s.replace('\n#', '\n\n#')
        with open(str(log_file), 'w') as f:
            f.write(s)

    return add_log





