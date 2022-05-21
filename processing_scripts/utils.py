# -*- coding: utf-8 -*-


import numpy as np
import gsw


def ceil_to_2_power_x(y):
    return 2**int(np.ceil(np.log2(y)))


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
        n = ceil_to_2_power_x(T/dt)
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


def spherical_polar_gradient(f, lon, lat, r=6371000.0):
    """Extension of the np.gradient function to spherical polar coordinates.
    Only gradients on a surface of constant radius (i.e. 2 dimensional) are
    currently supported. The grid must be evenly spaced in latitude and
    longitude.

    Important
    ---------
    For f(i, j), this function assumes i denotes position in latitude and j the
    position in longitude.

    Parameters
    ----------
    f : 2d array
        Scalar to calculate gradient.
    lon : 1d array
        Longitude points. [Degrees]
    lat : 1d array
        Latitude points. [Degrees]
    r : float
        Radius of sphere, defaults to Earth radius, 6371000 m.


    Returns
    -------
    dflon: 2d array
        Derivative in longitudinal direction.
    dflat: 2d array
        Derivative in the latitudinal direction.

    """
    nr, nc = f.shape
    if (nr != len(lat)) or (nc != len(lon)):
        raise ValueError(
            "Latitude and longitude are expected to be rows and" "columns respectively"
        )

    lon, lat = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    dfi = np.gradient(f, axis=0)
    dfj = np.gradient(f, axis=1)
    dlon = np.gradient(lon, axis=1)
    dlat = np.gradient(lat, axis=0)

    # Cosine because latitude from -90 to 90. Not 0 to pi.
    dfdlon = dfj / (r * dlon * np.cos(lat))
    dfdlat = dfi / (r * dlat)

    return dfdlon, dfdlat


def spherical_polar_gradient_ts(f, lon, lat, r=6371000.0):
    """Gradient of a two dimensional time series.

    Important
    ---------
    For f(i, j, k), this function assumes i denotes time, j latitude and k
    longitude.

    See spherical_polar_gradient for details.

    """
    nt, nr, nc = f.shape
    if (nr != len(lat)) or (nc != len(lon)):
        raise ValueError

    lon, lat = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    dfi = np.gradient(f, axis=1)
    dfj = np.gradient(f, axis=2)
    dlon = np.gradient(lon, axis=1)
    dlat = np.gradient(lat, axis=0)

    dlon = np.tile((r * dlon * np.cos(lat))[np.newaxis, ...], (nt, 1, 1))
    dlat = np.tile((r * dlat)[np.newaxis, ...], (nt, 1, 1))

    dfdlon = dfj / dlon
    dfdlat = dfi / dlat

    return dfdlon, dfdlat


def spherical_polar_area(r, lon, lat):
    """Calculates the area bounding an array of latitude and longitude points.

    Parameters
    ----------
    r : float
        Radius of sphere.
    lon : 1d array
        Longitude points. [Degrees]
    lat : 1d array
        Longitude points. [Degrees]

    Returns
    -------
    areas: 2d array

    """

    mid_dlon = (lon[2:] - lon[:-2]) / 2.0
    s_dlon = lon[1] - lon[0]
    e_dlon = lon[-1] - lon[-2]
    dlon = np.hstack((s_dlon, mid_dlon, e_dlon))

    mid_dlat = (lat[2:] - lat[:-2]) / 2.0
    s_dlat = lat[1] - lat[0]
    e_dlat = lat[-1] - lat[-2]
    dlat = np.hstack((s_dlat, mid_dlat, e_dlat))

    dlon, dlat = np.deg2rad(dlon), np.deg2rad(dlat)

    gdlon, gdlat = np.meshgrid(dlon, dlat)

    solid_angle = gdlon.T * gdlat.T * np.cos(np.deg2rad(lat))

    return solid_angle.T * r ** 2