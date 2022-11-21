#!/usr/bin/env python 
# Convert external/hourly_GPS_1.04.mat to a netcdf file

import numpy as np
from scipy import io
import xarray as xr


dat = io.loadmat("external/hourly_GPS_1.04.mat", simplify_cells=True)

use = np.isfinite(dat["ID"])

ID = dat["ID"][use].astype("int32")
TIME = np.datetime64("1979-01-01") + dat["TIME"][use].astype(int).astype("timedelta64[h]")

data_vars = dict(
    DROGUE=("ID", dat["DROGUE"][use].astype(bool)),
    time=("ID", TIME),
    lon=("ID", dat["LON"][use].astype(np.float32)),
    lat=("ID", dat["LAT"][use].astype(np.float32)),
    lon_err=("ID", 1e-5*dat["LON_ERR"][use].astype(np.float32)),
    lat_err=("ID", 1e-5*dat["LAT_ERR"][use].astype(np.float32)),
    u=("ID", dat["U"][use].astype(np.float32)),
    v=("ID", dat["V"][use].astype(np.float32)),
    u_err=("ID", dat["U_ERR"][use].astype(np.float32)),
    v_err=("ID", dat["V_ERR"][use].astype(np.float32)),
    gap=("ID", dat["GAP"][use].astype(np.float32)),
    ID_unique=np.unique(ID)
)

coords = dict(
    ID=ID,
)

ds = xr.Dataset(data_vars, coords).to_netcdf("internal/hourly_GPS_1.04.nc")
