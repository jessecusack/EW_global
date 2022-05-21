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
# # Interpolate GLORYS properties to drifter tracks (the fast way)

# %%
import xarray as xr
import numpy as np
from tqdm import tqdm
import scipy.interpolate as itpl
import utils
import glob

# %% [markdown]
# Load drifter and reanalysis data.

# %%
gdp = xr.open_dataset("../data/internal/hourly_GPS_1.04.nc") # Drifters
files = sorted(glob.glob("../data/external/cmems_mod_glo_phy_my_0.083_P1D-m/uv15m/*.nc"))
glr = xr.open_mfdataset(files)  # GLORYS

# %% [markdown]
# The code below takes about 2 hrs to run on my 2021 M1 chip macbook pro.

# %%
ndat = 240  # Number of hours of data for each interpolation chunk

re = 6.3781e6  # Earth radius [m]
dt = np.max(np.diff(glr.time))  # Additional time bounds
dl = 2*np.maximum(np.max(np.diff(glr.latitude)), np.max(np.diff(glr.longitude)))  # Additional lon/lat bounds

# Float start indexes
idxs = np.hstack((0, np.nonzero(np.diff(gdp.ID.data))[0] + 1, gdp.ID.size))

# Store GLORYS quantites
u = np.zeros_like(gdp.u)
v = np.zeros_like(gdp.u)
dudx = np.zeros_like(gdp.u)
dudy = np.zeros_like(gdp.u)
dvdx = np.zeros_like(gdp.u)
dvdy = np.zeros_like(gdp.u)
successful_interp = np.full_like(gdp.u, False, dtype=bool)

for i in tqdm(range(gdp.ID_unique.size)):
    ID = gdp.ID_unique.data[i]
    # selecting with a slice into isel is by far the fasted way of getting data
    d = gdp.isel(ID=slice(idxs[i], idxs[i+1]))
    
    if d.time[0] < glr.time[0]:
        # print(f"Drifter {ID} existed before our reanalysis data. Skipping.")
        continue
        
    # Convert to same coordinate system as GLORYS
    lon = d.lon.data
    lon[lon > 180] -= 360
    d["lon"] = ("ID", lon)
    
    # Even smaller chunks of data are needed for accurate interpolation
    # Find indices for small chunks
    jdxs = np.arange(0, d.time.size + ndat, ndat)
    
    for j in range(jdxs.size - 1):
        lonc = d.lon.data[jdxs[j]:jdxs[j+1]].copy()
        latc = d.lat.data[jdxs[j]:jdxs[j+1]].copy()
        timec = d.time.data[jdxs[j]:jdxs[j+1]].copy()
    
        # Do we cross the 180 deg line?
        cross180 = np.any(np.abs(np.diff(lonc)) > 359.)
        # Does the GLORYS chunk need to cross the dateline?
        if not cross180:
            if lonc.min() < glr.longitude.data[0]:
                cross180 = True
            elif lonc.max() > glr.longitude.data[-1]:
                cross180 = True
        
        if cross180:
            # print(f"Drifter {ID}, index {i} at chunk {j} crosses (or gets close to) the 180 deg line")
            # Change our coordinate system so that the dateline is 0 deg
            west_side = lonc < 0
            lonc[west_side] += 180.
            lonc[~west_side] -= 180.
            
            lllon = lonc.min() - dl
            urlon = lonc.max() + dl
            
            lllat = latc.min() - dl
            urlat = latc.max() + dl
            tmin = timec[0] - dt
            tmax = timec[-1] + dt
            
            # gl in normal coordinates, slice in latitude and time first
            gln = glr.sel(latitude=slice(lllat, urlat), time=slice(tmin, tmax))
            
            # Now split and recombine with new coordinate system
            # The left side is the eastern most data
            dle = gln.sel(longitude=slice(180 + lllon, 181))
            dlw = gln.sel(longitude=slice(-181, urlon - 180))
            
            # Change coordinates
            dle = dle.assign_coords(longitude=dle.longitude - 180.)
            dlw = dlw.assign_coords(longitude=dlw.longitude + 180.)

            # Combine
            gl = xr.concat([dle, dlw], "longitude")

        elif not cross180:
            # Select from GLORYS around the drifter chunk
            lllon = lonc.min() - dl  # Lower left lon
            lllat = latc.min() - dl
            urlon = lonc.max() + dl  # Upper right lon
            urlat = latc.max() + dl
            tmin = timec[0] - dt
            tmax = timec[-1] + dt
            gl = glr.sel(longitude=slice(lllon, urlon), latitude=slice(lllat, urlat), time=slice(tmin, tmax))

        # Mid lon and lat for coordinate transformation
        mlon = 0.5*(urlon + lllon)
        mlat = 0.5*(urlat + lllat)

        # Convert to distance in m, approximately and time in seconds
        x = re * np.deg2rad(gl.longitude.data - mlon) * np.cos(np.deg2rad(mlat))
        y = re * np.deg2rad(gl.latitude.data - mlat)
        t = (gl.time.data - tmin).astype("timedelta64[s]").astype(float)

        xd = re * np.deg2rad(lonc - mlon) * np.cos(np.deg2rad(mlat))
        yd = re * np.deg2rad(latc - mlat)
        td = (timec - tmin).astype("timedelta64[s]").astype(float)

        # Calculate spherical gradients
        dudx_, dudy_ = utils.spherical_polar_gradient_ts(gl.uo.data, gl.longitude.data, gl.latitude.data, r=re)
        dvdx_, dvdy_ = utils.spherical_polar_gradient_ts(gl.vo.data, gl.longitude.data, gl.latitude.data, r=re)

        # Interpolate velocity and its gradients
        sl = slice(idxs[i] + jdxs[j], np.minimum(idxs[i] + jdxs[j+1], idxs[i+1]))

        # For some reason, this is faster than interpn
        u[sl] = itpl.RegularGridInterpolator((t, y, x), gl.uo.data.compute())((td, yd, xd))
        v[sl] = itpl.RegularGridInterpolator((t, y, x), gl.vo.data.compute())((td, yd, xd))
        dudx[sl] = itpl.RegularGridInterpolator((t, y, x), dudx_.compute())((td, yd, xd))
        dudy[sl] = itpl.RegularGridInterpolator((t, y, x), dudy_.compute())((td, yd, xd))
        dvdx[sl] = itpl.RegularGridInterpolator((t, y, x), dvdx_.compute())((td, yd, xd))
        dvdy[sl] = itpl.RegularGridInterpolator((t, y, x), dvdy_.compute())((td, yd, xd))

        # If we made it here, great!
        successful_interp[sl] = True


# %%
data_vars = {
    "u": ("ID", u),
    "v": ("ID", v),
    "dudx": ("ID", dudx),
    "dudy": ("ID", dudy),
    "dvdx": ("ID", dvdx),
    "dvdy": ("ID", dvdy),
    "time": ("ID", gdp.time.data),
    "interp_flag": ("ID", successful_interp),
}

coords = dict(ID = ("ID", gdp.ID.data))

ds = xr.Dataset(data_vars, coords)
ds.to_netcdf("../data/internal/hourly_GPS_1.04_GLORYS_variables.nc")

# %%
ds

# %%
print(f"{successful_interp.sum()/successful_interp.size:.3%} of data sucessfully interpolated")
