
# coding: utf-8

# # Write GNOME ugrid output file from Vdatum ADCIRC tidal constituent files

# The code below using functionality from both `iris` and `libgoods`. You can install `iris` with conda (it's on the `conda-forge` channel), but unfortunately the `libgoods` library is not installable via `conda` or `pypi`, so this is how I installed it:
# ```
# git clone https://github.com/NOAA-ORR-ERD/GnomeTools
# pip install ./GnomeTools/libgoods
# ```

# In[1]:

from __future__ import print_function

import os
import warnings
import numpy as np

import iris

from libgoods import tri_grid
from libgoods import data_files_dir


# In[2]:

#ncfile = os.path.join(data_files_dir,'vdatum','vdatum_fl_sab_grid.nc')


# In[3]:

ncfile = ('http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/'
              'NYsndbght02_adcirc54.nc')
print(ncfile)
sl = 40.5457896; wl = -73.664; nl = 40.6990759; el = -73.3376574
ofile = 'fire_island.nc'


# In[4]:

ncfile = ('http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/'
         'vdatum_gulf_of_maine/adcirc54_38_orig.nc')
print(ncfile)
wl = -70.7234; el = -70.4532; sl = 41.4258; nl = 41.5643  # Vineyard sound 2.
ofile = 'gulf_of_maine.nc'


# In[5]:

ncfile = ('http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/'
              'vdatum_fl_sab/adcirc54.nc')
print(ncfile)
wl = -82; el = -80.6; sl = 29; nl = 30.6;
ofile = 'fl_sab.nc'


# In[6]:

ncfile = ('http://geoport.whoi.edu/thredds/dodsC/usgs/vault0/models/tides/'
              'DEdelches01_adcirc54.nc')
print(ncfile)
wl = -74.537378; el = -74.0315462; sl = 39.354624; nl = 39.704567 # Great Bay, NJ
ofile = 'great_bay.nc'


# In[7]:

var_map = {'latitude':'lat','longitude':'lon','nodes_surrounding_ele':'ele','eles_surrounding_ele':''}
adcirc = tri_grid.ugrid(ncfile)
adcirc.get_dimensions(var_map, get_time=False)
adcirc.get_grid_topo(var_map)

adcirc.find_nodes_eles_in_ss(nl, sl, wl, el)


# In[8]:

print(sl)


# In[9]:

bnd = adcirc.find_bndry_segs(subset=True)
print('Size of boundary: ', len(bnd))


# In[10]:

seg_types = [0] * len(bnd)
adcirc.order_boundary(bnd,seg_types)


# In[11]:

def parse_string(name):
    lista = [e.decode().strip() for e in name.tolist()]
    return ''.join(lista)


# In[12]:

#adcirc.update(ncfile)

names = []
const = adcirc.Dataset.variables['tidenames'][:]
for name in const:
    names.append(parse_string(name))
    


# In[13]:

from utide import _ut_constants_fname
from utide.utilities import loadbunch

con_info = loadbunch(_ut_constants_fname)['const']


# In[14]:

k = 0
ind_nc, ind_ttide = [], []

const_name = [e.strip() for e in con_info['name'].tolist()]

consts = ['STEADY', 'M2', 'S2', 'N2', 'K1', 'O1', 'P1', 'M4', 'M6']
for name in consts:
    try:
        if name == 'STEADY':
            indx = const_name.index('Z0')
        else:
            indx = const_name.index(name)
        k += 1
        ind_ttide.append(indx)
        ind_nc.append(names.index(name))
    except ValueError:
        pass  # `const` not found.
 


# In[15]:

import pytz
from datetime import datetime
from pandas import date_range

start = datetime.strptime('18-Sep-2015 05:00',
                          '%d-%b-%Y %H:%M').replace(tzinfo=pytz.utc)
stop = datetime.strptime('19-Sep-2015 05:00',  # '18-Sep-2015 18:00'
                         '%d-%b-%Y %H:%M').replace(tzinfo=pytz.utc)
dt = 1.0  # Hours.
glocals = date_range(start, stop, freq='1H').to_pydatetime()
ntimes = len(glocals)


# In[16]:

lon = adcirc.data['lon'][:]
lat = adcirc.data['lat'][:]


# In[17]:

inbox = np.logical_and(np.logical_and(lon >= wl,
                                      lon <= el),
                       np.logical_and(lat >= sl,
                                      lat <= nl))


# In[18]:

lon = lon[inbox]
lat = lat[inbox]


# In[19]:

print(inbox.shape)
print(inbox)
print(adcirc.nodes_in_ss.shape)
print(adcirc.nodes_in_ss)


# In[20]:

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cubes = iris.load_raw(ncfile)


# In[21]:

ua = cubes.extract_strict('Eastward Water Velocity Amplitude')
up = cubes.extract_strict('Eastward Water Velocity Phase')
va = cubes.extract_strict('Northward Water Velocity Amplitude')
vp = cubes.extract_strict('Northward Water Velocity Phase')


# In[22]:

uamp = ua.data[0, inbox, :][:, ind_nc]
vamp = va.data[0, inbox, :][:, ind_nc]
upha = up.data[0, inbox, :][:, ind_nc]
vpha = vp.data[0, inbox, :][:, ind_nc]


# In[23]:

# this takes 12 minutes instead of Iris 14 seconds
#uamp = adcirc.Dataset.variables['u_amp'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#vamp = adcirc.Dataset.variables['v_amp'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#upha = adcirc.Dataset.variables['u_phase'][0,adcirc.nodes_in_ss,:][:,ind_nc]
#vpha = adcirc.Dataset.variables['v_phase'][0,adcirc.nodes_in_ss,:][:,ind_nc]


# In[24]:

freq_nc = adcirc.Dataset.variables['tidefreqs'][:][ind_nc] 
freq_ttide = con_info['freq'][ind_ttide]
t_tide_names = np.array(const_name)[ind_ttide]
omega_ttide = 2*np.pi * freq_ttide  # Convert from radians/s to radians/hour.
omega = freq_nc * 3600


# In[25]:

from matplotlib.dates import date2num
from utide.harmonics import FUV
v, u, f = FUV(t=np.array([date2num(start)]), tref=np.array([0]),
              lind=np.array([ind_ttide]),
              lat=55, ngflgs=[0, 0, 0, 0])
              
# Convert phase in radians.
v, u, f = (np.squeeze(i) for i in (v, u, f))
v = v * 2 * np.pi
u = u * 2 * np.pi

thours = np.array([d.total_seconds() for d in
                   (glocals - glocals[0])]) / 60 / 60.
   
adcirc.data['u'] = np.ones((len(thours),len(adcirc.nodes_in_ss)),)   
adcirc.data['v'] = np.ones((len(thours),len(adcirc.nodes_in_ss)),)


# In[26]:

k = 0
for k in range(len(thours)):          
    U = (f * uamp * np.cos(v + thours[k] * omega + u - upha * np.pi/180)).sum(axis=1)
    V = (f * vamp * np.cos(v + thours[k] * omega + u - vpha * np.pi/180)).sum(axis=1) 
    adcirc.data['u'][k,:] = U
    adcirc.data['v'][k,:] = V


# In[27]:

adcirc.data['time'] = thours
adcirc.atts['time'] = {'units':'hours since 2015-09-18 05:00'}
adcirc.atts['u'] = {'units':'m/s','_FillValue':999}
adcirc.atts['v'] = {'units':'m/s','_FillValue':999}
adcirc.atts['nbe'] = {'order':'ccw'}
adcirc.write_unstruc_grid(ofile)                           

