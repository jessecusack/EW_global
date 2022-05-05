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

# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import getpass
import copernicusmarine as cm
import numpy as np
import scipy.signal as sig
import gsw

# %%
ds = xr.open_dataset("test.nc")
ds["time"] = np.datetime64("1970-01-01T00:00:00") + ds.time.data.astype(np.timedelta64)
dt = (ds.time[-1] - ds.time[0]).data.astype("timedelta64[s]").astype(float)

print(f"Record length = {dt/86400:1.2f} days")

ds

# %%
ds.u.plot()

# %%
ds.v.isel(depth=0).plot()

# %%
nperseg = 2**12

dt = (ds.time[1] - ds.time[0]).data.astype("timedelta64[s]").astype(float)
fs = 1/dt
freq, Puu = sig.welch(ds.u.isel(depth=0), fs, window="hann", nperseg=nperseg)
freq, Pvv = sig.welch(ds.v.isel(depth=0), fs, window="hann", nperseg=nperseg)

fig, ax = plt.subplots()
ax.loglog(freq, Puu)
ax.loglog(freq, Pvv)

freq_cor = gsw.f(ds.lat.data)/(np.pi*2)  # Hz
ax.axvline(freq_cor)

ax.set_xlabel("Frequency [Hz]")

print(f"Coriolis period = {1/(freq_cor*3600)} hr")

# %% [markdown]
# ## Extract altimetry

# %%
USERNAME = 'jcusack1'
PASSWORD = getpass.getpass('Enter your password: ')

# %%
DATASET_ID = "dataset-uv-rep-hourly"
ds1 = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD))
ds1

# %%
lon0 = ds.lon.data
lat0 = ds.lat.data
lon_slice = slice(lon0 - 2, lon0 + 2)
lat_slice = slice(lat0 - 2, lat0 + 2)
depth_level = 1
lon_step = 1
lat_step = 1
time_slice = slice(ds.time.data[0], ds.time.data[-1])
itime_plot = -1

ds_ = ds1.isel(depth=depth_level, longitude=slice(0, None, lon_step), latitude=slice(0, None, lat_step))
ds_ = ds_.sel(longitude=lon_slice, latitude=lat_slice, time=time_slice)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))

# make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_extent((lon_slice.start, lon_slice.stop, lat_slice.start, lat_slice.stop))
ax.coastlines()
ax.gridlines(draw_labels=True)

Q = ax.quiver(ds_.longitude, ds_.latitude, ds_.uo.isel(time=itime_plot), ds_.vo.isel(time=itime_plot), transform=ccrs.PlateCarree())
plt.quiverkey(Q, 0.75, 1, 0.2, label="0.2 m s$^{-1}$", coordinates="figure", color="red")

ax.plot(ds.lon, ds.lat, 'ro')

# %%
dl = 0.25
lon0 = ds.lon.data
lat0 = ds.lat.data
lon_slice = slice(lon0 - dl, lon0 + dl)
lat_slice = slice(lat0 - dl, lat0 + dl)
depth_level = 1
lon_step = 1
lat_step = 1
time_slice = slice(ds.time.data[0], ds.time.data[-1])

ds_ = ds1.isel(depth=depth_level).sel(longitude=lon_slice, latitude=lat_slice).sel(time=time_slice)
uo = ds_.uo.interp(longitude=lon0, latitude=lat0)
vo = ds_.vo.interp(longitude=lon0, latitude=lat0)

# %%
idepth_obs = 0
nroll = 48

fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10, 5))
axs[0].plot(ds.time, ds.u.isel(depth=idepth_obs).rolling(time=nroll, center=True).mean(), label="obs")
axs[0].plot(ds_.time, uo, label="alt")

axs[1].plot(ds.time, ds.v.isel(depth=idepth_obs).rolling(time=nroll, center=True).mean())
axs[1].plot(ds_.time, vo)

axs[0].legend()

# %% [markdown]
# Glorys Reanalysis

# %%
DATASET_ID = "cmems_mod_glo_phy_my_0.083_P1D-m"

# cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M
ds3 = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD))
ds3

# %%
lon0 = ds.lon.data
lat0 = ds.lat.data
lon_slice = slice(lon0 - 2, lon0 + 2)
lat_slice = slice(lat0 - 2, lat0 + 2)
depth_instrument = ds.depth[0]
lon_step = 1
lat_step = 1
time_slice = slice(ds.time.data[0], ds.time.data[-1])
itime_plot = -1

ds_ = ds3.isel(longitude=slice(0, None, lon_step), latitude=slice(0, None, lat_step)).sel(depth=depth_instrument, method="nearest")
ds_ = ds_.sel(longitude=lon_slice, latitude=lat_slice, time=time_slice)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))

# make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_extent((lon_slice.start, lon_slice.stop, lat_slice.start, lat_slice.stop))
ax.coastlines()
ax.gridlines(draw_labels=True)

Q = ax.quiver(ds_.longitude, ds_.latitude, ds_.uo.isel(time=itime_plot), ds_.vo.isel(time=itime_plot), transform=ccrs.PlateCarree())
plt.quiverkey(Q, 0.75, 1, 0.2, label="0.2 m s$^{-1}$", coordinates="figure", color="red")

ax.plot(ds.lon, ds.lat, 'ro')

# %%
dl = 0.25
lon0 = ds.lon.data
lat0 = ds.lat.data
lon_slice = slice(lon0 - dl, lon0 + dl)
lat_slice = slice(lat0 - dl, lat0 + dl)
depth_level = 1
lon_step = 1
lat_step = 1
time_slice = slice(ds.time.data[0], ds.time.data[-1])
depth_instrument = ds.depth[0]

ds_ = ds3.sel(longitude=lon_slice, latitude=lat_slice, time=time_slice).sel(depth=depth_instrument, method="nearest")
uo = ds_.uo.interp(longitude=lon0, latitude=lat0)
vo = ds_.vo.interp(longitude=lon0, latitude=lat0)

# %%
idepth_obs = 0
nroll = 48

fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10, 5))
axs[0].plot(ds.time, ds.u.isel(depth=idepth_obs).rolling(time=nroll, center=True).mean(), label="obs")
axs[0].plot(ds_.time, uo, label="reanalysis")

axs[1].plot(ds.time, ds.v.isel(depth=idepth_obs).rolling(time=nroll, center=True).mean())
axs[1].plot(ds_.time, vo)

axs[0].legend()
