
from multiprocessing import set_start_method
set_start_method('forkserver')


# In[2]:

import version_control


# In[3]:

from spgr import lg
from logging import INFO, DEBUG
lg.setLevel(INFO)


# ## Electrode Locations

# ### Locations for all the subjects

# In[4]:

from phypno.attr import Freesurfer, Channels

from spgr.constants import (IMAGES_PATH,
                            REC_PATH,
                            FS_FOLDER,
                            HEMI_SUBJ,
                            CHAN_TYPE,
                            DATA_OPTIONS,
                            )
from spgr.read_data import get_chan_used_in_analysis


# In[6]:

images_dir = IMAGES_PATH.joinpath('electrode_locations')
if not images_dir.exists():
    images_dir.mkdir(parents=True)

PLOT_COLOR = 'kw'


# In[7]:

from phypno.viz import Viz3
from IPython.display import display
SKIN_COLOR = (239, 208, 207, 240)


# In[8]:

for subj, hemi in HEMI_SUBJ.items():
    chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, **DATA_OPTIONS)
    fs = Freesurfer(str(REC_PATH.joinpath(subj).joinpath(FS_FOLDER)))
    surf = getattr(fs.read_brain(), hemi)

    v = Viz3(color='kw')
    v.add_chan(chan, color=(200, 20, 20, 255))
    v.add_surf(surf, color=SKIN_COLOR)
    v.save(str(images_dir.joinpath(subj + '.png')))
    display(v)


# ## Coverage

# In[9]:

from numpy import ones
from spgr.constants import HEMI_SUBJ, CHAN_TYPE, DATA_OPTIONS
from spgr.read_data import get_chan_used_in_analysis
from spgr.spindle_source import get_morph_linear, lg
from spgr.plot_spindles import plot_surf


# In[10]:

morphed = []

for subj in HEMI_SUBJ:

    # set channel values to 1
    chan = get_chan_used_in_analysis(subj, 'sleep', CHAN_TYPE, **DATA_OPTIONS)
    values = ones(chan.n_chan)

    morphed.append(get_morph_linear(subj, values, reref='avg'))

v = plot_surf(morphed, limits=(0, 1), extra_smoothing=False)
v.save(str(images_dir.joinpath('coverage_average.png')))
v

