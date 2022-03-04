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

# %% [markdown]
# Help site: https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python#h_33df7ebcce
#
# Copernicus login info.

# %%
USERNAME = 'jcusack1'
PASSWORD = getpass.getpass('Enter your password: ')

# %% [markdown]
# ## Multi-obs UV inc. Ekman
#
# https://resources.marine.copernicus.eu/product-detail/MULTIOBS_GLO_PHY_REP_015_004/INFORMATION

# %%
DATASET_ID = "dataset-uv-rep-hourly"
ds1 = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD))
ds1

# %%
lon_slice = slice(-85, -55)
lat_slice = slice(30, 50)
depth_level = 0
lon_step = 1
lat_step = 1
date = "2011-09-05T12:00"

ds_ = ds1.isel(depth=depth_level, longitude=slice(0, None, lon_step), latitude=slice(0, None, lat_step))
ds_ = ds_.sel(time=date, method="nearest").sel(longitude=lon_slice, latitude=lat_slice)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))

# make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_extent((lon_slice.start, lon_slice.stop, lat_slice.start, lat_slice.stop))
ax.coastlines()
ax.gridlines(draw_labels=True)

ax.quiver(ds_.longitude, ds_.latitude, ds_.uo, ds_.vo, transform=ccrs.PlateCarree())

# %% [markdown]
# ## Sea level derived geostrophic UV
#
# https://resources.marine.copernicus.eu/product-detail/SEALEVEL_GLO_PHY_L4_MY_008_047/INFORMATION

# %%
DATASET_ID = "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D" # "cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1D"

# cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.25deg_P1M
ds2 = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD))
ds2

# %%
lon_slice = slice(-85, -55)
lat_slice = slice(30, 50)
lon_step = 1
lat_step = 1
date = "2011-09-05T12:00"

ds_ = ds2.isel(longitude=slice(0, None, lon_step), latitude=slice(0, None, lat_step))
ds_ = ds_.sel(time=date, method="nearest").sel(longitude=lon_slice, latitude=lat_slice)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))

# make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_extent((lon_slice.start, lon_slice.stop, lat_slice.start, lat_slice.stop))
ax.coastlines()
ax.gridlines(draw_labels=True)

ax.quiver(ds_.longitude, ds_.latitude, ds_.ugos, ds_.vgos, transform=ccrs.PlateCarree())


# %% [markdown]
# ## Cartopy example

# %%
def sample_data(shape=(20, 30)):
    """
    Return ``(x, y, u, v, crs)`` of some vector data
    computed mathematically. The returned crs will be a rotated
    pole CRS, meaning that the vectors will be unevenly spaced in
    regular PlateCarree space.

    """
    crs = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)

    x = np.linspace(311.9, 391.1, shape[1])
    y = np.linspace(-23.6, 24.8, shape[0])

    x2d, y2d = np.meshgrid(x, y)
    u = 10 * (2 * np.cos(2 * np.deg2rad(x2d) + 3 * np.deg2rad(y2d + 30)) ** 2)
    v = 20 * np.cos(6 * np.deg2rad(x2d))

    return x, y, u, v, crs



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(-10, 45))

ax.add_feature(cfeature.OCEAN, zorder=0)
ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')

ax.set_global()
ax.gridlines()

x, y, u, v, vector_crs = sample_data()
ax.quiver(x, y, u, v, transform=vector_crs)

# plt.show()
