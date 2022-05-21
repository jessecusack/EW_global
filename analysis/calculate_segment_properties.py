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
# # Quantities (velocity, stress, etc) on drifter segments

# %% [markdown]
# Eddy wave decomposition is done spectrally, e.g. 
#
# \begin{equation}
# Q = \int_{\omega_0}^{\omega_1} S(\omega) d\omega
# \end{equation}
#
# Where Q is the integrated quantity (e.g. energy in a given frequency band), computed by integrating the spectrum (e.g energy spectrum) over the frequencies associated with the process of interest. In the case of eddies $\omega_0 = 0$ and $\omega_1 = 0.25$ cpd. For near-inertial waves, we can take $\omega_0 = 0.75f$ and $\omega_1 = 1.5f$, for example.

# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
import cartopy
import cartopy.crs as ccrs
from scipy import stats
import scipy.signal as sig
import gsw
import scipy.integrate as itgr
import utils
from tqdm import tqdm

M2 = 24/12.4206012  # M2 frequency in cycles per day
K1 = 24/23.93447213  # K1 frequency

# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# Global parameters

# %%
fhi = 1.25  # coriolis upper cut off
flo = 1/fhi  # Coriolis lower cut off
fe = 0.25  # Eddy cut off

# %% [markdown]
# ## Even segment analysis

# %%
dss = xr.open_dataset("../data/internal/hourly_GPS_1.04_evenly_segmented.nc")

# %% [markdown]
# Mean total energy...

# %%
us = dss.u.data
vs = dss.v.data

KE_av = np.average(0.5*(us**2 + vs**2), axis=1, weights=sig.hann(dss.sample.size))

dss["KE_mean"] = ("segment", KE_av)

# %%
trans = ccrs.PlateCarree()

fig, ax = plt.subplots(1, 1, figsize=(25, 8), subplot_kw=dict(projection=ccrs.Robinson()))
scm = ax.scatter(dss.lon_mean, dss.lat_mean, s=2, c=dss.KE_mean, transform=trans, vmin=0, vmax=0.25)
fig.colorbar(scm, ax=ax)
ax.coastlines()
ax.add_feature(cartopy.feature.LAND)

# %%
dl = 2
lon_bins = np.arange(0, 360 + dl, dl)
lat_bins = np.arange(-80, 90 + dl, dl)

KE_mean_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.KE_mean, bins=[lon_bins, lat_bins])
# KE_mean_bin[:, np.abs(utils.mid(lat_bins)) < 10] = np.nan

counts, _, _ = np.histogram2d(dss.lon_mean, dss.lat_mean, [lon_bins, lat_bins], density=False)
counts[counts == 0] = np.nan

# %%
trans = ccrs.PlateCarree()

fig, axs = plt.subplots(2, 1, figsize=(25, 16), subplot_kw=dict(projection=ccrs.Robinson()))
pccm = axs[0].pcolormesh(lon_bins, lat_bins, KE_mean_bin.T, transform=trans, vmin=0, vmax=0.25)
fig.colorbar(pccm, ax=axs[0])
axs[0].coastlines()
axs[0].add_feature(cartopy.feature.LAND)


cmap, norm = from_levels_and_colors([1, 2, 5, 10, 20, 50], ['C0', 'C1', 'C2', 'C3', 'C4', 'C5'], extend="max")
pccm = axs[1].pcolormesh(lon_bins, lat_bins, counts.T, transform=trans, norm=norm, cmap=cmap)
fig.colorbar(pccm, ax=axs[1])
axs[1].coastlines()
axs[1].add_feature(cartopy.feature.LAND)

# %% [markdown]
# ### Compute spectra
#
# Energy integrated in frequency...

# %%
welch_kwargs = dict(fs=24, window="hann", nperseg=dss.sample.size//2, detrend=False, axis=-1)

freq, Puu = sig.welch(us, **welch_kwargs)
_, Pvv = sig.welch(vs, **welch_kwargs)
_, Cuv = sig.csd(us, vs, **welch_kwargs)
Cuv = Cuv.real

# %%
dss["freq"] = freq
dss["Suu"] = (["segment", "freq"], Puu)
dss["Svv"] = (["segment", "freq"], Pvv)
dss["Cuv"] = (["segment", "freq"], Cuv)
dss["SKE"] = (["segment", "freq"], 0.5*(Puu + Pvv))

# %%
i = 84939

seg = dss.isel(segment=i)

fcor = np.abs(seg.f_mean)*86400/(np.pi*2)

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlabel("Frequency [cpd]")
ax.set_ylabel("KE spectral density [J kg$^{-1}$ cpd$^{-1}$]")

# ax.loglog(freq, Puu[i, :])
# ax.loglog(freq, Pvv[i, :])
ax.loglog(seg.freq, seg.SKE)
ax.axvline(fcor, color="k", label="$f$")
ax.axvline(M2, color="r", label="$M_2$")
ax.axvline(K1, color="b", label="$K_1$")
ax.set_ylim(1e-6, 1)
ax.set_xlim(1e-2, 10)

ax.fill_betweenx(ax.get_ylim(), 2*[flo*fcor], 2*[fhi*fcor], alpha=0.3, color="gray")
ax.fill_betweenx(ax.get_ylim(), 2*[0], 2*[fe], alpha=0.3, color="gray")

ax.legend()

step = 10
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection=ccrs.Mercator()))
ax.plot(seg.lon, seg.lat, "-o", transform=ccrs.PlateCarree(), ms=1)
ax.quiver(seg.lon[::step], seg.lat[::step], seg.u[::step], seg.v[::step], transform=ccrs.PlateCarree(), color="k")
ax.quiver(seg.lon[::step], seg.lat[::step], seg.ug[::step], seg.vg[::step], transform=ccrs.PlateCarree(), color="r")
ax.gridlines(draw_labels=True, linewidth=0.5, color='k', alpha=0.5, linestyle='--')

fig, ax = plt.subplots()
ax.plot(np.diff(seg.time).astype("timedelta64[h]"), "o", ms=1)
ax.set_ylabel("Diff of time")

# Integrate eddy
freqi = np.linspace(0, fe, 50)
Ei = np.interp(freqi, seg.freq, seg.SKE)
Ee = itgr.simpson(Ei, freqi)

# Integrate NI
freqi = np.linspace(flo*fcor, fhi*fcor, 50)
Ei = np.interp(freqi, seg.freq, seg.SKE)
ENI = itgr.simpson(Ei, freqi)

print(f"ID {seg.ID.data}")
print(f"Time gap detected? {seg.time_flag.data}")
print(f"Mean position = {seg.lat_mean.data:.1f} N, {seg.lon_mean.data:.1f} E")
print(f"Eddy energy = {Ee:.6f} J kg-1")
print(f"NI energy = {ENI:.6f} J kg-1")
print(f"Ratio Ee/EIW = {Ee/ENI:.4f}")

# %% [markdown]
# ### Kinetic energy in frequency bands
#
# Near inertial KE and eddy KE for all spectra...

# %%
SKE = dss.SKE.data
fcor = np.abs(dss.f_mean)*86400/(np.pi*2)
freq = dss.freq.data
freqie = np.linspace(0, fe, 50)

KEe = []
KEIW = []

for i in tqdm(dss.segment.data):
    # Integrate eddy
    Ei = np.interp(freqie, freq, SKE[i])
    KEe.append(itgr.simpson(Ei, freqi))

    # Integrate NI
    freqi = np.linspace(flo*fcor[i], fhi*fcor[i], 50)
    Ei = np.interp(freqi, freq, SKE[i])
    KEIW.append(itgr.simpson(Ei, freqi))
    
dss["KEe"] = ("segment", KEe)
dss["KEIW"] = ("segment", KEIW)

# %%
dl = 2
lon_bins = np.arange(0, 360 + dl, dl)
lat_bins = np.arange(-90, 90 + dl, dl)

mask = np.abs(utils.mid(lat_bins)) < 7.5

KEe_mean_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.KEe, bins=[lon_bins, lat_bins])
KEe_mean_bin[:, mask] = np.nan

KEIW_mean_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.KEIW, bins=[lon_bins, lat_bins])
KEIW_mean_bin[:, mask] = np.nan

G_mean_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.KEe/dss.KEIW, bins=[lon_bins, lat_bins])
G_mean_bin[:, mask] = np.nan

# %%
trans = ccrs.PlateCarree()
fontsize = 14

fig, axs = plt.subplots(3, 1, figsize=(25, 24), subplot_kw=dict(projection=ccrs.Robinson()))

pccm = axs[0].pcolormesh(lon_bins, lat_bins, KEe_mean_bin.T, transform=trans, vmin=0, vmax=3e-1)
cb = fig.colorbar(pccm, ax=axs[0])
cb.set_label("Low frequency KE [J kg$^{-1}$]", fontsize=fontsize)

pccm = axs[1].pcolormesh(lon_bins, lat_bins, KEIW_mean_bin.T, transform=trans, vmin=0, vmax=2e-2)
cb = fig.colorbar(pccm, ax=axs[1])
cb.set_label("Near inertial KE [J kg$^{-1}$]", fontsize=fontsize)

pccm = axs[2].pcolormesh(lon_bins, lat_bins, np.log10(G_mean_bin.T), transform=trans, vmin=-2, vmax=2, cmap="RdBu_r")
cb = fig.colorbar(pccm, ax=axs[2])
cb.set_label("$\log_{10}(E_{eddy}/E_{NI})$", fontsize=fontsize)

for ax in axs:

    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND)
    
    
fig.savefig("../figures/low_freq_to_NI_energy_ratio.pdf", dpi=180, bbox_inches="tight", pad_inches=0.01)

# %% [markdown]
# ### Stress
#
# Integrate stresses in NI band:
#
# Shear stress
#
# \begin{equation}
# \overline{u^\prime v^\prime} = \int_{0.75}^{1.5f} C_{uv} (\omega) d \omega
# \end{equation}
#
# Normal stress
#
# \begin{equation}
# \overline{u^\prime u^\prime} - \overline{v^\prime v^\prime}  = \int_{0.75f}^{1.5f} P_{uu} - P_{vv} (\omega) d \omega
# \end{equation}

# %%
Suu = dss.Suu.data
Svv = dss.Svv.data
Cuv = dss.Cuv.data
freq = dss.freq.data
fcor = np.abs(dss.f_mean)*86400/(np.pi*2)

Suui = []
Svvi = []
Cuvi = []

for i in tqdm(dss.segment.data):
    # Integrate NI
    freqi = np.linspace(flo*fcor[i], fhi*fcor[i], 50)
    Suui.append(itgr.simpson(np.interp(freqi, freq, Suu[i]), freqi))
    Svvi.append(itgr.simpson(np.interp(freqi, freq, Svv[i]), freqi))
    Cuvi.append(itgr.simpson(np.interp(freqi, freq, Cuv[i]), freqi))
    
dss["Suui"] = ("segment", Suui)
dss["Svvi"] = ("segment", Svvi)
dss["Cuvi"] = ("segment", Cuvi)

dss["nstress"] = dss.Suui - dss.Svvi
dss["sstress"] = -2*dss.Cuvi

# %% [markdown]
# What does this stress look like?

# %%
trans = ccrs.PlateCarree()

fig, axs = plt.subplots(2, 1, figsize=(25, 16), subplot_kw=dict(projection=ccrs.Robinson()))
scm = axs[0].scatter(dss.lon_mean, dss.lat_mean, s=2, c=dss.nstress, transform=trans, vmin=-1e-3, vmax=1e-3, cmap="PiYG")
fig.colorbar(scm, ax=axs[0])
scm = axs[1].scatter(dss.lon_mean, dss.lat_mean, s=2, c=dss.sstress, transform=trans, vmin=-1e-3, vmax=1e-3, cmap="PiYG")
fig.colorbar(scm, ax=axs[1])

for ax in axs:
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND)

# %%
lon_bins = np.arange(0, 362, 2)
lat_bins = np.arange(-90, 92, 2)

mask = np.abs(utils.mid(lat_bins)) < 7.5

nstress_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.nstress, bins=[lon_bins, lat_bins])
nstress_bin[:, mask] = np.nan

sstress_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.sstress, bins=[lon_bins, lat_bins])
sstress_bin[:, mask] = np.nan

# %%
trans = ccrs.PlateCarree()

fig, axs = plt.subplots(2, 1, figsize=(25, 16), subplot_kw=dict(projection=ccrs.Robinson()))
pccm = axs[0].pcolormesh(lon_bins, lat_bins, nstress_bin.T, transform=trans, vmin=-1e-3, vmax=1e-3, cmap="PiYG")
cb = fig.colorbar(pccm, ax=axs[0])
cb.set_label("Normal stress [N m$^{-2}$]", fontsize=fontsize)

pccm = axs[1].pcolormesh(lon_bins, lat_bins, sstress_bin.T, transform=trans, vmin=-1e-3, vmax=1e-3, cmap="PiYG")
cb = fig.colorbar(pccm, ax=axs[1])
cb.set_label("Shear stress [N m$^{-2}$]", fontsize=fontsize)

for ax in axs:
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND)

# %% [markdown]
# ## Strain rates from GLORYS
#

# %%
i = 91943
seg = dss.isel(segment=i)

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
ax.plot(seg.time, seg.ug, "C0")
ax.plot(seg.time, seg.u, "C0:")
ax.plot(seg.time, seg.vg, "C1")
ax.plot(seg.time, seg.v, "C1:")

# %%
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(seg.time, seg.nstrain, "C0", label="normal strain")
ax.plot(seg.time, seg.sstrain, "C1", label="shear strain")
ax.legend()

# %% [markdown]
# # Energy Transfer

# %%
nstrain = dss.nstrain.data
sstrain = dss.sstrain.data

nstrain_mean = np.average(nstrain, axis=1, weights=sig.hann(dss.sample.size))
sstrain_mean = np.average(sstrain, axis=1, weights=sig.hann(dss.sample.size))

dss["nstrain_mean"] = ("segment", nstrain_mean)
dss["sstrain_mean"] = ("segment", sstrain_mean)

dss["Fn"] = dss.nstrain_mean*dss.nstress
dss["Fs"] = dss.sstrain_mean*dss.sstress
dss["F"] = dss.Fn + dss.Fs

# %%
lon_bins = np.arange(0, 362, 2)
lat_bins = np.arange(-90, 92, 2)

mask = np.abs(utils.mid(lat_bins)) < 7.5

Fn_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.Fn, bins=[lon_bins, lat_bins])
Fn_bin[:, mask] = np.nan

Fs_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.Fs, bins=[lon_bins, lat_bins])
Fs_bin[:, mask] = np.nan

F_bin, _, _, binnumber = stats.binned_statistic_2d(dss.lon_mean, dss.lat_mean, dss.F, bins=[lon_bins, lat_bins])
F_bin[:, mask] = np.nan

# %%
trans = ccrs.PlateCarree()
fontsize = 14

fig, axs = plt.subplots(3, 1, figsize=(25, 24), subplot_kw=dict(projection=ccrs.Robinson()))

pccm = axs[0].pcolormesh(lon_bins, lat_bins, Fn_bin.T, transform=trans, vmin=-1e-9, vmax=1e-9, cmap="RdBu_r")
cb = fig.colorbar(pccm, ax=axs[0])
# cb.set_label("Low frequency KE [J kg$^{-1}$]", fontsize=fontsize)

pccm = axs[1].pcolormesh(lon_bins, lat_bins, Fs_bin.T, transform=trans, vmin=-1e-9, vmax=1e-9, cmap="RdBu_r")
cb = fig.colorbar(pccm, ax=axs[1])
# cb.set_label("Near inertial KE [J kg$^{-1}$]", fontsize=fontsize)

pccm = axs[2].pcolormesh(lon_bins, lat_bins, F_bin.T, transform=trans, vmin=-1e-9, vmax=1e-9, cmap="RdBu_r")
cb = fig.colorbar(pccm, ax=axs[2])
# cb.set_label("$\log_{10}(E_{eddy}/E_{NI})$", fontsize=fontsize)

for ax in axs:

    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND)
    
    
# fig.savefig("../figures/low_freq_to_NI_energy_ratio.pdf", dpi=180, bbox_inches="tight", pad_inches=0.01)

fig, ax = plt.subplots()
ax.plot(np.nanmean(F_bin, axis=0), utils.mid(lat_bins))

# %% [markdown]
# ## Uneven segment analysis (slow)

# %%
# ds = xr.open_dataset("../data/internal/hourly_GPS_1.04.nc")
# segs = xr.open_dataset("../data/internal/segments_uneven_hourly_GPS_1.04.nc")

# %% [markdown]
# Loop method for estimating segment data

# %%
# Estimate segment data via loop

# from tqdm import tqdm

# segs["lon_mean"] = ("segment", np.full(segs.segment.shape, np.NaN, "float32"), dict(long_name="mean longitude"))
# segs["lat_mean"] = ("segment", np.full(segs.segment.shape, np.NaN, "float32"), dict(long_name="mean latitude"))
# # segs["time"] = ("segment", np.full(segs.segment.shape, np.NaN, "float32"), dict(long_name="mean time"))
# # segs["KE_mean"] = ("segment", np.full(segs.segment.shape, np.NaN, "float32"), dict(long_name="mean total KE"))


# for i in tqdm(segs.segment.data[:10]):
#     dsi = ds.isel(ID=slice(segs.start_idx[i].data, segs.start_idx[i].data + segs.length[i].data))
#     segs.lon_mean[i] = dsi.lon.mean().data
#     segs.lat_mean[i] = dsi.lat.mean().data


# %% [markdown]
# Reduce the data for tests...

# %%
# nsegmax = 1000

# ds = ds.set_coords(["segment"])
# insegment = (ds.segment.data > -1) #& (ds.segment.data < nsegmax)
# # drogued = ds.DROGUE.data
# dsr = ds.isel(ID=insegment)

# %%
# KE = 0.5*(dsr.u**2 + dsr.v**2)

# %%
# dsm = dsr.groupby("segment").mean("ID")

# %%
# KE_tot_mean = KE.groupby("segment").mean("ID")

# %%
# KE_seg_mean = 0.5*(dsm.u**2 + dsm.v**2)

# %%
# trans = ccrs.PlateCarree()

# fig, axs = plt.subplots(3, 1, figsize=(25, 20), subplot_kw=dict(projection=ccrs.Robinson()))

# scm = axs[0].scatter(dsm.lon, dsm.lat, s=2, c=KE_tot_mean, transform=trans)
# fig.colorbar(scm, ax=axs[0])
# scm = axs[1].scatter(dsm.lon, dsm.lat, s=2, c=KE_seg_mean, transform=trans)
# fig.colorbar(scm, ax=axs[1])
# scm = axs[2].scatter(dsm.lon, dsm.lat, s=2, c=KE_tot_mean - KE_seg_mean, transform=trans, cmap="RdBu_r", vmin=-5e-2, vmax=5e-2)
# fig.colorbar(scm, ax=axs[2])

# for ax in axs:
#     # ax.set_global()
#     ax.coastlines()
    
#     # ax.stock_img()

