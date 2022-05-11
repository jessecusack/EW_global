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

# %%
ds = xr.open_dataset("../data/internal/hourly_GPS_1.04.nc")
ds

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
# Great, the IDs are all listed in increasing order. Now find the segments

# %%
latamin = 10  # Absolute minimum latitude
ceil_to_2powx = False

def win_length(f, IP=4, dt=3600., ceil_to_2powx=True):
    """
    Parameters
    ----------
        f : float
            Coriolis frequency [rad s-1]
        NI : int
            Number of inertial periods needed in window.
        dt : float, optional
            Sampling period [s]
        ceil_to_2powx : float, option
            Do ceiling of result to nearest 2**x where x is an integer (e.g. for use with fft).
            
    returns
        n : int
            Number of samples in window.
            
    """
    
    # Inertial period [s]
    Tf = np.pi*2/np.abs(f)
    # Time needed [s]
    T = IP*Tf
    
    if ceil_to_2powx:
        n = np.ceil(T/dt)
        x = int(np.ceil(np.log2(n)))
        n = 2**x
    else:
        n = int(np.ceil(T/dt))
    
    return n


def find_start_idx(lats, latamin):
    """
    Parameters
    ----------
        lats : ndarray
            Array of latitudes [degrees_north]
        latamin : float
            Absolute minimum latitude [degrees]
            
    returns
        idx : int
            Number of samples in window.
            
    """
    
    lowlat = np.abs(lats[0]) < latamin
    
    if not lowlat:
        return 0
    
    idx = np.argmax(np.abs(lats) > latamin)
    if idx == 0:
        raise RuntimeError("No latitude is above the absolute minimum")

    return idx
    
def check_data_remaining(lats, **win_length_kwargs):
    """
    Parameters
    ----------
        lats : ndarray
            Array of latitudes [degrees_north]
        **win_length_kwargs : dict, optional
            Arguments to win_length.
            
    returns
        n : int
            Number of samples in window from function win_length.
        data_remains : bool
            True if enough data remains, otherwise False.
            
    """
    
    f = gsw.f(lats[0])
    n = win_length(f, **win_length_kwargs)
    return n, lats.size > n


# Unique float start indices:
idxs = np.hstack((0, np.nonzero(np.diff(ds.ID.data))[0] + 1, ds.ID.size))

segment_start_idxs = []
segment_lengths = []
segment_ID = np.full(ds.ID.shape, -1, "int32")

for i in tqdm(range(ds.ID_unique.size)):
    ID = ds.ID_unique[i].data
    # selecting with a slice into isel is by far the fasted way of getting data
    dsi = ds.isel(ID=slice(idxs[i], idxs[i+1]))
    
    # print(f"Float {ID} has {dsi.ID.size} samples and starts at {dsi.lon[0].data:1.1f} E, {dsi.lat[0].data:1.1f} N")
    
    try:
        ic = find_start_idx(dsi.lat.data, latamin)
    except RuntimeError:
        # print("  All data is equatorial, skipping")
        continue
        
    for j in range(50000):  # This range sets an upper limit
        nwin, enough_data_remains = check_data_remaining(dsi.lat[ic:].data, ceil_to_2powx=ceil_to_2powx)
        
        if enough_data_remains:
            segment_start_idxs.append(idxs[i] + ic)
            segment_lengths.append(nwin)
            segment_ID[idxs[i] + ic:idxs[i] + ic + nwin] = len(segment_lengths) - 1
            ic += nwin
            
            try:
                ic += find_start_idx(dsi.lat[ic:].data, latamin)
            except RuntimeError:
                # print(f"  Found {j} segments")
                # print("  Remaining data are equatorial, moving on")
                break
            
        else:
            # print(f"  Found {j} segments")
            # print("  Insufficient data remaining, moving on.")
            break
            
            
segs = xr.Dataset(dict(start_idx=("segment", np.asarray(segment_start_idxs, "int32")), length=("segment", np.asarray(segment_lengths, "int32"))))

ds["segment"] = ("ID", segment_ID, dict(long_name="segment number"))

# %% [markdown]
# Save outputs.

# %%
segs.to_netcdf("../data/internal/segments_hourly_GPS_1.04.nc")
# To overwrite existing file, have to load data into ram...
ds = ds.load()
ds.close()
ds.to_netcdf("../data/internal/hourly_GPS_1.04.nc")  

# %% [markdown]
# How many segments did we find?

# %%
segs.segment.size

# %% [markdown]
# Do the lengths and start indexes make sense?

# %%
fig, ax = plt.subplots()
ax.plot(segs.start_idx, '.')

fig, ax = plt.subplots()
ax.plot(segs.length, '.')

# %% [markdown]
# Do any segments cover more than one float?

# %%
for i0, n in tqdm(zip(segs.start_idx.data, segs.length.data)):
    if np.unique(ds.ID[i0:i0+n]).size > 1:
        raise RuntimeError(f"Oh no! Segment with start index = {i0} covers more than one float ID")

print("All is well...")
    

# %% [markdown]
# Function for checking data for a particular segment.

# %%
def check_segment(i, ds, start_idxs, lengths):
    i0 = start_idxs[i]
    n = lengths[i]
    dsi = ds.isel(ID=slice(i0, i0+n))
    dsi = dsi.assign_coords(dict(time=dsi.time))
    
    print(f"np.unique(ID) = {np.unique(dsi.ID)}")
    
    fig, axs = plt.subplots(3, 1, figsize=(16, 20))
    axs[0].plot(dsi.lon, dsi.lat, "-o", ms=3)
    axs[0].plot(dsi.lon[0], dsi.lat[0], "go", ms=10)
    axs[0].plot(dsi.lon[-1], dsi.lat[-1], "ro", ms=10)
    dsi.u.plot(ax=axs[1], x="time", label="u")
    dsi.v.plot(ax=axs[1], x="time", label="v")
    axs[1].legend()
    dsi.u_err.plot(ax=axs[2], x="time", label="u_err")
    dsi.v_err.plot(ax=axs[2], x="time", label="v_err")


# %%
check_segment(5678, ds, segs.start_idx.data, segs.length.data)
