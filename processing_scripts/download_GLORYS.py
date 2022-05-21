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
# # Download GLORYS data
#
# https://resources.marine.copernicus.eu/product-detail/GLOBAL_MULTIYEAR_PHY_001_030/INFORMATION

# %%
import xarray as xr
import getpass
import copernicusmarine as cm
import numpy as np

# %%
USERNAME = 'jcusack1'
PASSWORD = getpass.getpass('Enter your password: ')

# %%
DATASET_ID = "cmems_mod_glo_phy_my_0.083_P1D-m"

# Drop as many unnecessary variables as possible
ds = xr.open_dataset(cm.copernicusmarine_datastore(DATASET_ID, USERNAME, PASSWORD), drop_variables=["bottomT", "sithick", "siconc", "usi", "vsi", "mlotst", "zos", "thetao", "so"])
ds = ds.sel(time=slice("2007-01-01", None)).sel(depth=15, method="nearest")
ds

# %% [markdown]
# According to an error message when doing `ds.uo.data.size`:
#
# `message = "Request too big=131498.52768 Mbytes, max=1500.0";`
#
# Which is 131 GB for just one variable.

# %%
print(f"{4320*2041*3804*32 / 8 / 2**30:.1f} GiB")

# %%
from tqdm import tqdm
import os
save_dir = "../data/external/cmems_mod_glo_phy_my_0.083_P1D-m/uv15m/"

# %%
for month in tqdm(np.unique(np.datetime_as_string(ds.time, "M"))): 
    dsi = ds.sel(time=month)
    filename = f"uv15m_{month}.nc"
    path = os.path.join(save_dir, filename)
    if not os.path.isfile(path):
        dsi.to_netcdf(path)
