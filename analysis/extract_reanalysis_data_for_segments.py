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
# # Extract reanalysis data around drifter tracks

# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import getpass
import copernicusmarine as cm
import numpy as np

# %%
USERNAME = 'jcusack1'
PASSWORD = getpass.getpass('Enter your password: ')

# %% [markdown]
# https://resources.marine.copernicus.eu/product-detail/GLOBAL_MULTIYEAR_PHY_001_030/INFORMATION

# %%
DATASET_ID = "cmems_mod_glo_phy_my_0.083_P1D-m"

# cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M
ds = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD), drop_variables=["bottomT", "sithick", "siconc", "usi", "vsi"])
ds

# %% [markdown]
# ## A single drifter segment

# %%
gdp = xr.open_dataset("../data/internal/hourly_GPS_1.04_evenly_segmented.nc")
gdp

# %%
seg = gdp.isel(segment=12430)
seg

# %%
fig, axs = plt.subplots(3, 1, figsize=(16, 20))
axs[0].plot(seg.lon, seg.lat, "-o", ms=3)
axs[0].plot(seg.lon[0], seg.lat[0], "go", ms=10)
axs[0].plot(seg.lon[-1], seg.lat[-1], "ro", ms=10)
axs[1].plot(seg.time, seg.u, label="u")
axs[1].plot(seg.time, seg.v, label="v")
axs[1].legend()
axs[2].plot(seg.time, seg.u_err, label="u_err")
axs[2].plot(seg.time, seg.v_err, label="v_err")



# %%
# Convert lon
lon = seg.lon.data
lon[lon > 180] -= 360
seg["lon"] = ("sample", lon)

# %%
dl = 0.25
dt = np.timedelta64(12, "h")
lllon = seg.lon.min().data - dl
lllat = seg.lat.min().data - dl
urlon = seg.lon.max().data + dl
urlat = seg.lat.max().data + dl
tmin = seg.time[0].data - dt
tmax = seg.time[-1].data + dt

# %%
ds_ = ds.sel(depth=15, method="nearest").sel(longitude=slice(lllon, urlon), latitude=slice(lllat, urlat), time=slice(tmin, tmax))
ds_

# %% [markdown]
# Strain?

# %%
import utils

# %%
dudx, dudy = utils.spherical_polar_gradient_ts(ds_.uo.data, ds_.longitude.data, ds_.latitude.data)
dvdx, dvdy = utils.spherical_polar_gradient_ts(ds_.vo.data, ds_.longitude.data, ds_.latitude.data)

ds_["dudx"] = (["time", "latitude", "longitude"], dudx)
ds_["dudy"] = (["time", "latitude", "longitude"], dudy)
ds_["dvdx"] = (["time", "latitude", "longitude"], dvdx)
ds_["dvdx"] = (["time", "latitude", "longitude"], dvdx)
ds_["nstrain"] = (["time", "latitude", "longitude"], dudx - dvdy)
ds_["sstrain"] = (["time", "latitude", "longitude"], dvdx + dudy)


# %%
ds_

# %%
dsi = ds_.interp(dict(longitude=seg.lon, latitude=seg.lat, time=seg.time))

# %%
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(dsi.time, dsi.uo, "C0")
ax.plot(seg.time, seg.u, "C0:")
ax.plot(dsi.time, dsi.vo, "C1")
ax.plot(seg.time, seg.v, "C1:")

# %%
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(dsi.time, dsi.nstrain, "C0", label="normal strain")
ax.plot(dsi.time, dsi.sstrain, "C1", label="shear strain")
ax.legend()

# %%
