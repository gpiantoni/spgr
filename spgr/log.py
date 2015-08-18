from datetime import datetime
from logging import (DEBUG,
                     INFO,
                     FileHandler,
                     Formatter,
                     getLogger,
                     StreamHandler,
                     )

from spgr.constants import LOGSRC_PATH, PROJECT, IMAGES_PATH

from base64 import b64encode
from re import sub
from socket import gethostname
from subprocess import check_output
from shutil import rmtree
import phypno
import spgr


def embed_images_in_html(html_file):
    """read images from png file and embed them into the html.

    Parameters
    ----------
    html_file : path to file
        path to html file

    """
    with open(html_file, 'r') as f:
        s = f.read()

    s1 = sub('<img src="([a-zA-Z0-9_/\.]*)" ', _embed_png, s)

    with open(html_file, 'w') as f:
        f.write(s1)


def _embed_png(matched):
    """Take a regex object with img tag and convert the png path to base64 data.

    Parameters
    ----------
    matched : regex match object
        matched regex of the img tag

    Returns
    -------
    str
        string to replace the whole img tag.
    """
    string = matched.group(0)
    image_path = matched.group(1)

    with open(image_path, 'rb') as f:
        image_data = b64encode(f.read()).decode()

    return string.replace(image_path, 'data:image/png;base64,' + image_data)


def git_hash(package):
    package_path = package.__path__[0] + '/../.git'
    version = check_output('git --git-dir ' + package_path +
                           ' rev-parse HEAD', shell=True).decode('utf-8')
    return version


def git_branch(package):
    package_path = package.__path__[0] + '/../.git'
    branch = check_output('git --git-dir ' + package_path +
                          ' rev-parse --abbrev-ref HEAD',
                          shell=True).decode('utf-8')
    return branch[:-1]  # remove new line


def git_info(lg):

    lg.info('{}: {} {} '.format(PROJECT, git_branch(spgr), git_hash(spgr)))
    lg.info('{}: {} {} '.format('phypno', git_branch(phypno), git_hash(phypno)))
    lg.info(check_output('pip freeze', shell=True).decode('utf-8'))


def with_log(function):

    def add_log():
        log_file = LOGSRC_PATH.joinpath(function.__name__ + '.md')
        if log_file.exists():
            log_file.unlink()

        images_dir = IMAGES_PATH.joinpath(function.__name__)
        try:
            rmtree(str(images_dir))
        except FileNotFoundError:
            pass
        images_dir.mkdir(parents=True)

        lg = getLogger(PROJECT)
        lg.setLevel(DEBUG)

        formatter = Formatter('%(message)s  ')

        lg.handlers = []

        fh = FileHandler(str(log_file))
        fh.setLevel(INFO)
        fh.setFormatter(formatter)
        lg.addHandler(fh)

        ch = StreamHandler()
        ch.setLevel(DEBUG)
        ch.setFormatter(formatter)
        lg.addHandler(ch)

        lg.info('# {}'.format(function.__name__.replace('_', ' ')))

        git_info(lg)

        lg.info('## Started')
        t0 = datetime.now()
        lg.info('{} on {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S'),
                                  gethostname()))

        function(lg, images_dir)

        lg.info('## Finished')
        t1 = datetime.now()
        lg.info('{} after {}'.format(t1.strftime('%Y-%m-%d %H:%M:%S'),
                                     str(t1 - t0)[:-7]))
        fh.close()

        with open(str(log_file), 'r') as f:
            s = f.read()
        s = s.replace('\n#', '\n\n#')
        with open(str(log_file), 'w') as f:
            f.write(s)

    return add_log
