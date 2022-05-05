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
# # First look at drifer data

# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import utils
import gsw

ds = xr.open_dataset("../data/driftertrajGPS_1.03.nc")
ds = ds.assign_coords({"BID": ds.ID})

# %%
ds

# %%
IDs = np.unique(ds.ID.values[np.isfinite(ds.ID.values)])
fig, axs = plt.subplots(1, 2)
axs[0].plot(IDs, '.')
axs[1].hist(IDs[IDs > 5e7]);

print(IDs[IDs > 5e7][::20])

# %%
ds_ = ds.isel(TIME=ds.BID==64734970.)
ds_

# %%
fcor = gsw.f(ds_.LAT.min())  # minimum coriolis
Tcor = 2*np.pi/np.abs(fcor)
Ts = 60.*60.  # Sampling period
Tlpf = 1.1*Tcor  # Low pass period

ul = utils.butter_filter(ds_.U, 1/Tlpf, 1/Ts)
vl = utils.butter_filter(ds_.V, 1/Tlpf, 1/Ts)

fig, ax = plt.subplots(1, 1)
ax.plot(ds_.LON, ds_.LAT)
ax.plot(ds_.LON[0], ds_.LAT[0], "go")
ax.plot(ds_.LON[-1], ds_.LAT[-1], "ro")

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
axs[0].plot(ds_.TIME, ds_.U)
axs[0].plot(ds_.TIME, ul)
axs[1].plot(ds_.TIME, ds_.V)
axs[1].plot(ds_.TIME, vl)
axs[2].plot(ds_.TIME, ds_.U_ERR, ".")

fig.autofmt_xdate()

# %%
# zoom in on a chunk
i0 = 100
i1 = 500

cut = slice(i0, i1)

fig, ax = plt.subplots(1, 1)
ax.plot(ds_.LON[cut], ds_.LAT[cut])
ax.plot(ds_.LON[cut][0], ds_.LAT[cut][0], "go")
ax.plot(ds_.LON[cut][-1], ds_.LAT[cut][-1], "ro")

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
axs[0].plot(ds_.TIME[cut], ds_.U[cut])
axs[0].plot(ds_.TIME[cut], ul[cut])
axs[1].plot(ds_.TIME[cut], ds_.V[cut])
axs[1].plot(ds_.TIME[cut], vl[cut])
axs[2].plot(ds_.TIME[cut], ds_.U_ERR[cut], ".")

fig.autofmt_xdate()

# %%
