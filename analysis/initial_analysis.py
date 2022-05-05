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
import scipy.signal as sig

T_M2 = 24/12.4206012  # M2 tidal period in cycles per day

ds = xr.open_dataset("../data/internal/hourly_GPS_1.04.nc")
ds

# %%
fig, ax = plt.subplots()
ax.plot(ds.lon, ds.lat, '.', ms=0.2)

# %%
ds_ = ds.sel(ID = ds.ID_unique[100].data)

print(ds_.time[0].data)
print(ds_.time[-1].data)

ds_

# %%
fcor = gsw.f(ds_.lat.min())  # minimum coriolis
Tcor = 2*np.pi/np.abs(fcor)
Ts = 60.*60.  # Sampling period
Tlpf = 2*Tcor  # Low pass period

ul = utils.butter_filter(ds_.u, 1/Tlpf, 1/Ts)
vl = utils.butter_filter(ds_.v, 1/Tlpf, 1/Ts)

fig, ax = plt.subplots(1, 1)
ax.plot(ds_.lon, ds_.lat)
ax.plot(ds_.lon[0], ds_.lat[0], "go")
ax.plot(ds_.lon[-1], ds_.lat[-1], "ro")

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
axs[0].plot(ds_.time, ds_.u)
axs[0].plot(ds_.time, ul)
axs[1].plot(ds_.time, ds_.v)
axs[1].plot(ds_.time, vl)
axs[2].plot(ds_.time, ds_.u_err, ".")

fig.autofmt_xdate()

# %%
# zoom in on a chunk
i0 = 3000
i1 = i0 + 256

cut = slice(i0, i1)

fig, ax = plt.subplots(1, 1)
ax.plot(ds_.lon[cut], ds_.lat[cut], ls="-", marker=".")
ax.plot(ds_.lon[cut][0], ds_.lat[cut][0], "go")
ax.plot(ds_.lon[cut][-1], ds_.lat[cut][-1], "ro")

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 9))
axs[0].plot(ds_.time[cut], ds_.u[cut])
axs[0].plot(ds_.time[cut], ul[cut])
axs[1].plot(ds_.time[cut], ds_.v[cut])
axs[1].plot(ds_.time[cut], vl[cut])
axs[2].plot(ds_.time[cut], ds_.u_err[cut], ".")

fig.autofmt_xdate()

print(f"Days in chunk: {(i1 - i0)/24}")

# %%
nperseg = 256

freq, Puu = sig.csd(ds_.u[i0:i1], ds_.u[i0:i1], fs=24, nperseg=nperseg)
_, Pvv = sig.csd(ds_.v[i0:i1], ds_.v[i0:i1], fs=24, nperseg=nperseg)
_, Puv = sig.csd(ds_.u[i0:i1], ds_.v[i0:i1], fs=24, nperseg=nperseg)
freq2, Puu2 = sig.csd(ds_.u[i0:i1], ds_.u[i0:i1], fs=24, nperseg=nperseg/2)
_, Pvv2 = sig.csd(ds_.v[i0:i1], ds_.v[i0:i1], fs=24, nperseg=nperseg/2)
_, Puv2 = sig.csd(ds_.u[i0:i1], ds_.v[i0:i1], fs=24, nperseg=nperseg/2)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(5, 10))
axs[0].loglog(freq, Puu.real, "C0:")
axs[0].loglog(freq, Pvv.real, "C1:")
axs[0].loglog(freq2, Puu2.real, "C0")
axs[0].loglog(freq2, Pvv2.real, "C1")
axs[1].semilogx(freq, Puv.real, "C0:")
axs[1].semilogx(freq2, Puv2.real, "C0")

f_mean = gsw.f(ds_.lat[i0:i1].mean())*86400/(np.pi*2)

for ax in axs:
    ax.axvline(f_mean, color="k")
    ax.axvline(1/T_M2, color="gray")

# %%
