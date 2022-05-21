# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: ewglobal
#     language: python
#     name: ewglobal
# ---

# %% [markdown]
# # Segment the drifter data
#
# The drifter data consist of many tracks from individual drifters. For each drifter, we want to split its track into several segments to be analysed separately. Each segment should be sufficiently long to estimate some quantities of both the near-inertial wave field and the eddy field. We might want to capture 4 inertial periods, for example.
#
# A problem is that we don't know how many segments we need to create and that segment length will depend on latitude.

# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import gsw
import utils
from scipy import stats

# %load_ext autoreload
# %autoreload 2

# %%
ds = xr.open_dataset("../data/internal/hourly_GPS_1.04.nc")
glorys = xr.open_dataset("../data/internal/hourly_GPS_1.04_GLORYS_variables.nc")

# %% [markdown]
# Are all of the data listed by float order?

# %%
fig, ax = plt.subplots()
ax.plot(np.abs(np.diff(ds.ID[:100000])) > 0)

fig, ax = plt.subplots()
ax.plot(ds.ID)

# %%
# always increasing ID number?
(np.diff(ds.ID) >= 0).all()

# %% [markdown]
# Find the indices for each unique float.

# %%
# idxs = np.hstack((0, np.nonzero(np.diff(ds.ID.data))[0] + 1, ds.ID.size))
# idxs.shape

# %% [markdown]
# ## Even segments
#
# How big to make the segments? Well they need to capture at least a few inertial periods. The longest inertial period is at our low latitude cut off, 10 S.

# %%
Tmax = np.pi*2/gsw.f(10)
print(Tmax)

# %% [markdown]
# Sample period is

# %%
Ts = 60*60

# %% [markdown]
# How many inertial periods do we want?

# %%
IP = 5

# %% [markdown]
# How many samples do we need?

# %%
n = int(np.ceil(Tmax*IP/Ts))
print(f"Number of days = {n*Ts/86400:.1f}")
print(f"Number of samples = {n}\n")

n2powx = utils.ceil_to_2_power_x(n)
print(f"Nearest 2**x = {n2powx}")
print(f"Total time {n2powx*Ts/86400:.1f} days")
print(f"Min number of inertial periods = {n2powx*Ts/Tmax:.1f}")


# %% [markdown]
# Create segments...

# %%
n = n2powx  # Number of samples
latamin = 10  # Absolute minimum latitude
overlap = n//2  # Segment overlap
skip_low_lat = False  # Flat to apply absolute minimum latitude limit

# Unique float start indices:
idxs = np.hstack((0, np.nonzero(np.diff(ds.ID.data))[0] + 1, ds.ID.size))

segment_start_idxs = []
segment_ID = np.full(ds.ID.shape, -1, "int32")

for i in tqdm(range(ds.ID_unique.size)):
    ID = ds.ID_unique[i].data
    # selecting with a slice into isel is by far the fasted way of getting data
    dsi = ds.isel(ID=slice(idxs[i], idxs[i+1]))
    
    # print(f"Float {ID} has {dsi.ID.size} samples and starts at {dsi.lon[0].data:1.1f} E, {dsi.lat[0].data:1.1f} N")
    
    try:
        if skip_low_lat:
            ic = utils.find_start_idx(dsi.lat.data, latamin)
        else:
            ic = 0    
    except RuntimeError:
        # print("  All data is equatorial, skipping")
        continue
        
    for j in range(50000):  # This range sets an upper limit, we should break out first!
        enough_data_remains = dsi.lat[ic:].data.size > n
        
        if enough_data_remains:
            # print(f"  Segment start index is {ic}")
            segment_start_idxs.append(idxs[i] + ic)
            ic += n - overlap
            # print(f"  Adding {n - overlap} to ic, ic is now {ic}")
            try:
                if skip_low_lat:
                    add = utils.find_start_idx(dsi.lat[ic:].data, latamin)
                    ic += add
                    # print(f"  Adding {add} to ic for latitude skip")
                else:
                    pass
                
            except (RuntimeError, IndexError):
                # print(f"  Found {j} segments")
                # print("  Remaining data are equatorial, moving on")
                break
            
        else:
            # print(f"  Found {j} segments")
            # print("  Insufficient data remaining, moving on.")
            break
            
            
segs_e = xr.Dataset(dict(start_idx=("segment", np.asarray(segment_start_idxs, "int32")), length=n, overlap=overlap), {}, dict(absolute_min_latitude=latamin, min_inertial_periods=n*Ts/Tmax))

# ds["segment_even"] = ("ID", segment_ID, dict(long_name="segment number", info=f"same size, segment size = {n}, absolute minimum latitude = {latamin}, min inertial periods = {n*Ts/Tmax:.1f}", missing_value=-1))

# %% [markdown]
# Save

# %%
segs_e.to_netcdf("../data/internal/segments_even_hourly_GPS_1.04.nc")

# %% [markdown]
# Do any segments cover more than one float?

# %%
n = segs_e.length.data
for i0 in tqdm(segs_e.start_idx.data):
    if np.unique(ds.ID[i0:i0+n]).size > 1:
        raise RuntimeError(f"Oh no! Segment with start index = {i0} covers more than one float ID")

print("All is well...")

# %% [markdown]
# # Split dataset into even segments

# %%
# Using numpy arrays makes this MUCH faster
u = ds.u.data
v = ds.v.data
u_err = ds.u_err.data
v_err = ds.v_err.data
lon = ds.lon.data
lat = ds.lat.data
lon_err = ds.lon_err.data
lat_err = ds.lat_err.data
time = ds.time.data
DROGUE = ds.DROGUE.data
ID = ds.ID.data
gap = ds.gap.data
n = segs_e.length.data
idxs = segs_e.start_idx.data

# GLORYS
ug = glorys.u.data
vg = glorys.v.data
dudx = glorys.dudx.data
dudy = glorys.dudy.data
dvdx = glorys.dvdx.data
dvdy = glorys.dvdy.data
interp_flag = glorys.interp_flag.data

data_vars = {
    "u": (["segment", "sample"], [u[i:i+n] for i in idxs]),
    "v": (["segment", "sample"], [v[i:i+n] for i in idxs]),
    "lon": (["segment", "sample"], [lon[i:i+n] for i in idxs]),
    "lat": (["segment", "sample"], [lat[i:i+n] for i in idxs]),
    "u_err": (["segment", "sample"], [u_err[i:i+n] for i in idxs]),
    "v_err": (["segment", "sample"], [v_err[i:i+n] for i in idxs]),
    "lon_err": (["segment", "sample"], [lon_err[i:i+n] for i in idxs]),
    "lat_err": (["segment", "sample"], [lat_err[i:i+n] for i in idxs]),
    "time": (["segment", "sample"], [time[i:i+n] for i in idxs]),
    "DROGUE": (["segment", "sample"], [DROGUE[i:i+n] for i in idxs]),
    # "gap": (["segment", "sample"], [gap[i:i+n] for i in idxs]),
    "ID": (["segment"], [ID[i] for i in idxs]),
    "ug": (["segment", "sample"], [ug[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "vg": (["segment", "sample"], [vg[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "dudx": (["segment", "sample"], [dudx[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "dudy": (["segment", "sample"], [dudy[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "dvdx": (["segment", "sample"], [dvdx[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "dvdy": (["segment", "sample"], [dvdy[i:i+n] for i in idxs], dict(info="GLORYS derived variable")),
    "interp_flag": (["segment", "sample"], [interp_flag[i:i+n] for i in idxs]),
}

dss = xr.Dataset(data_vars)
dss["ID_unique"] = np.unique(dss.ID.data)
dss

# %% [markdown]
# Calculate some more useful quantities.
#
# Shear strain
#
# \begin{equation}
# S_s = \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}
# \end{equation}
#
# Normal strain
#
# \begin{equation}
# S_n = \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}
# \end{equation}
#
# Divergence
#
# \begin{equation}
# D = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}
# \end{equation}
#
# Vorticity
#
# \begin{equation}
# \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
# \end{equation}

# %%
dss["sstrain"] = dss.dvdx + dss.dudy
dss["nstrain"] = dss.dudx - dss.dvdy
dss["div"] = dss.dudx + dss.dvdy
dss["vort"] = dss.dvdx - dss.dudy

dss["lon_mean"] = ("segment", stats.circmean(dss.lon, 360, 0, axis=1))
dss["lat_mean"] = ("segment", stats.circmean(dss.lat, 90, -90, axis=1))
dss["time_mean"] = dss.time.mean("sample")

dss["time_flag"] = ("segment", np.any(np.diff(dss.time, axis=1).astype("timedelta64[h]").astype(float) > 1, axis=1), dict(info="a gap in the time record was detected"))

dss["f_mean"] = gsw.f(dss.lat_mean)

# %% [markdown]
# Save

# %%
dss.to_netcdf("../data/internal/hourly_GPS_1.04_evenly_segmented.nc")

# %% [markdown]
# How many segments have gaps in time?

# %%
tgap = dss["time_flag"].data
print(f"{tgap.sum()/tgap.size:.1%} segments have a time gap")

# %% [markdown]
# ## Uneven segments (varying with $f$)
#
# This code is commented out because I don't use uneven segments. 

# %%
# IP = 4  # Number of inertial periods
# latamin = 10  # Absolute minimum latitude
# ceil_to_2powx = False

# # Unique float start indices:
# idxs = np.hstack((0, np.nonzero(np.diff(ds.ID.data))[0] + 1, ds.ID.size))

# segment_start_idxs = []
# segment_lengths = []
# segment_ID = np.full(ds.ID.shape, -1, "int32")

# for i in tqdm(range(ds.ID_unique.size)):
#     ID = ds.ID_unique[i].data
#     # selecting with a slice into isel is by far the fasted way of getting data
#     dsi = ds.isel(ID=slice(idxs[i], idxs[i+1]))
    
#     # print(f"Float {ID} has {dsi.ID.size} samples and starts at {dsi.lon[0].data:1.1f} E, {dsi.lat[0].data:1.1f} N")
    
#     try:
#         ic = utils.find_start_idx(dsi.lat.data, latamin)
#     except RuntimeError:
#         # print("  All data is equatorial, skipping")
#         continue
        
#     for j in range(50000):  # This range sets an upper limit
#         nwin, enough_data_remains = utils.check_data_remaining(dsi.lat[ic:].data, IP=IP, ceil_to_2powx=ceil_to_2powx)
        
#         if enough_data_remains:
#             segment_start_idxs.append(idxs[i] + ic)
#             segment_lengths.append(nwin)
#             segment_ID[idxs[i] + ic:idxs[i] + ic + nwin] = len(segment_lengths) - 1
#             ic += nwin
            
#             try:
#                 ic += utils.find_start_idx(dsi.lat[ic:].data, latamin)
#             except RuntimeError:
#                 # print(f"  Found {j} segments")
#                 # print("  Remaining data are equatorial, moving on")
#                 break
            
#         else:
#             # print(f"  Found {j} segments")
#             # print("  Insufficient data remaining, moving on.")
#             break
            
            
# segs_u = xr.Dataset(dict(start_idx=("segment", np.asarray(segment_start_idxs, "int32")), length=("segment", np.asarray(segment_lengths, "int32"))))

## ds["segment_uneven"] = ("ID", segment_ID, dict(long_name="segment number", info=f"uneven size, absolute minimum latitude = {latamin}, inertial periods = {IP}, ceiling to 2**x = {ceil_to_2powx}", missing_value=-1))

# %% [markdown]
# Save outputs.

# %%
# segs_u.to_netcdf("../data/internal/segments_uneven_hourly_GPS_1.04.nc")

# %% [markdown]
# How many segments did we find?

# %%
# segs_u.segment.size

# %% [markdown]
# Do the lengths and start indexes make sense?

# %%
# fig, ax = plt.subplots()
# ax.plot(segs_u.start_idx, '.')

# fig, ax = plt.subplots()
# ax.plot(segs_u.length, '.')

# %% [markdown]
# Do any segments cover more than one float?

# %%
# for i0, n in tqdm(zip(segs_u.start_idx.data, segs_u.length.data)):
#     if np.unique(ds.ID[i0:i0+n]).size > 1:
#         raise RuntimeError(f"Oh no! Segment with start index = {i0} covers more than one float ID")

# print("All is well...")

# %% [markdown]
# Function for checking data for a particular segment.

# %%
# def check_segment(i, ds, start_idxs, lengths):
#     i0 = start_idxs[i]
#     n = lengths[i]
#     dsi = ds.isel(ID=slice(i0, i0+n))
#     dsi = dsi.assign_coords(dict(time=dsi.time))
    
#     print(f"np.unique(ID) = {np.unique(dsi.ID)}")
    
#     fig, axs = plt.subplots(3, 1, figsize=(16, 20))
#     axs[0].plot(dsi.lon, dsi.lat, "-o", ms=3)
#     axs[0].plot(dsi.lon[0], dsi.lat[0], "go", ms=10)
#     axs[0].plot(dsi.lon[-1], dsi.lat[-1], "ro", ms=10)
#     dsi.u.plot(ax=axs[1], x="time", label="u")
#     dsi.v.plot(ax=axs[1], x="time", label="v")
#     axs[1].legend()
#     dsi.u_err.plot(ax=axs[2], x="time", label="u_err")
#     dsi.v_err.plot(ax=axs[2], x="time", label="v_err")

# %%
# check_segment(5678, ds, segs_u.start_idx.data, segs_u.length.data)
