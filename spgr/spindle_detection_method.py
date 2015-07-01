from phypno.detect import DetectSpindle
from phypno.detect.spindle import transform_signal
from phypno.trans import Select
from copy import deepcopy
from numpy import where, NaN
from pathlib import Path

from spgr.read_data import DATA_DIR, REC_DIR, FS_PATH, ELEC_PATH, get_data, GROUP_DIR

images_dir = IMAGES_PATH.joinpath('spindle_detection_method')
try:
    rmtree(str(images_dir))
except FileNotFoundError:
    pass
images_dir.mkdir(parents=True)

PLOT_COLOR = 'kw'

options = {'reref': 'avg',
           'resample_freq': 256,
           'hp_filter': 0.5,
           'lp_filter': 50,
           }
chan_type = ('grid', 'strip')


from .log import with_log

@with_log
def Spindle_Detection_Method(lg):


    lg.info('## Method description')
    lg.info('Spindles were detected on each electrode independently, using previously reported methods with stringent criteria [@Nir2011_regional].')
    lg.info('The raw signal was filtered between 9 and 16 Hz and the instantaneous amplitude was calculated from the analytic signal via Hilbert transform.')
    lg.info('The amplitude had to go above a threshold for detection (set at 3 times the S.D.).')
    lg.info('The beginning and end of the spindle were defined by a threshold for selection (set at 1 time the S.D.).')
    lg.info('Spindle duration had to be between 0.5 and 2s.')
    lg.info('We computed the power spectrum over the spindle interval and the peak in the power spectrum had to lie between 9 and 16 Hz (Fig. +[spgr_detect_method]).')

# In[7]:

LIMITS = (9, 16)

detsp = DetectSpindle(method='Nir2011', frequency=LIMITS, duration=(0.5, 2))


# In[8]:

PLOT_SIZE = 960, 240


# In[9]:

subj = 'EM09'
chan_name = 'GR28'
data = get_data(subj, 'sleep', chan_type, **options)
sel = Select(chan=(chan_name, ))
data = sel(data)


# In[10]:

good_trial = [i for i, trl in enumerate(data) if 42710 // 30 * 30 == trl.axis['time'][0][0]][0]
good_trial = 73


# In[11]:

sel = Select(trial=(good_trial, ))
sel_data = sel(data)


# In[12]:

from phypno.viz import Viz1


# In[13]:

select_beginning = Select(time=(None, 44852))
one_spindle_data = select_beginning(sel_data)


# In[20]:

v = Viz1(color=PLOT_COLOR)
v.size = PLOT_SIZE
v.add_data(one_spindle_data, limits_x=(44850, 44880), limits_y=(-200, 200))
v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')
v.save(str(images_dir.joinpath('raw_spindle.png' )))
v


# In[19]:

v = Viz1(color=PLOT_COLOR)
v.size = PLOT_SIZE
v.add_data(sel_data, limits_x=(44850, 44880), limits_y=(-200, 200))
v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')
v.save(str(images_dir.joinpath('raw.png' )))
v


# In[23]:

sp = detsp(sel_data)


# In[25]:

from phypno.trans import Filter, Math

filt = Filter(low_cut=detsp.det_butter['freq'][0], high_cut=detsp.det_butter['freq'][1],
              order=detsp.det_butter['order'], s_freq=data.s_freq)
x_filt = filt(sel_data)


# In[26]:

from numpy import ones, take, swapaxes, arange, atleast_2d
from scipy.signal import fftconvolve

def move_avg(dat, axis):
    """Function to compute moving average."""
    flat = ones(detsp.moving_avg['dur'] * data.s_freq)

    dat = swapaxes(dat, axis, -1)
    for i, x in enumerate(dat):
        dat[i, :] = fftconvolve(x, flat / sum(flat), mode='same')
    return swapaxes(dat, -1, axis)

hilbert_abs = Math(operator_name=('hilbert', 'abs'), axis='time')
moving_avg = Math(operator=move_avg, axis='time')

dat_det = moving_avg(hilbert_abs(x_filt))


# In[30]:

from pyqtgraph import mkPen

v = Viz1(color=PLOT_COLOR)
v.size = PLOT_SIZE
v.add_data(x_filt, color=mkPen('b', width=1), limits_x=(44850, 44880), limits_y=(-100, 100))
v.add_data(dat_det, color=mkPen('r', width=4))
v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')
v.save(str(images_dir.joinpath('envelope.png' )))
v


# In[33]:

t = dat_det.axis['time'][0]
good_sp = [x for x in sp.spindle if x['start_time'] < t[-1] and x['end_time'] > t[0]]

detected_spindles = deepcopy(dat_det)

for one_sp in good_sp:
    i0 = where(t == one_sp['start_time'])[0]
    i1 = where(t == one_sp['end_time'])[0]
    detected_spindles.data[0][0, i0:i1] = NaN


# In[35]:

v = Viz1(color=PLOT_COLOR)
v.size = PLOT_SIZE
v.add_data(sel_data, limits_x=(44850, 44880), limits_y=(-200, 200))
v._plots[chan_name].setLabels(left='amplitude (μV)', bottom='time (s)')
v.add_graphoelement(sp)
v.save(str(images_dir.joinpath('detected.png' )))
v


# In[36]:

from phypno.trans import Math, Filter
from numpy import pi


# In[38]:

sel = Select(trial=(73, ))
sel_data = sel(data)
sp = detsp(sel_data)

hp = Filter(low_cut=10, high_cut=15, s_freq=sel_data.s_freq)
h = Math(operator_name=('hilbert', 'angle', 'unwrap', 'diff'), axis='time')
ddata = h(hp(sel_data))
ddata.data[0] *= (ddata.s_freq / 2 / pi)

v = Viz1(color=PLOT_COLOR)
v.size = PLOT_SIZE
v.add_data(ddata, limits_x=(44850, 44880), limits_y=(5, 20))
v.add_graphoelement(sp)
v._plots[chan_name].setLabels(left='frequency (Hz)', bottom='time (s)')
v.save(str(images_dir.joinpath('inst_freq.png' )))
v


# ### All Subjects

# In[3]:

from spgr.constants import HEMI_SUBJ, SPINDLE_OPTIONS, CHAN_TYPE
from spgr.detect_spindles import get_spindles


# In[4]:

for ref in ('avg', 15):
    for subj in HEMI_SUBJ:
        print(subj)
        sp = get_spindles(subj, chan_type=CHAN_TYPE, reref=ref, **SPINDLE_OPTIONS)

