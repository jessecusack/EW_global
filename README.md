# Eddy-Internal Wave Interactions From Historical Moorings and Drifters

Internal wave stress estimated from moorings and drifters are coupled with strain from reanalysis products. 

Most of the python scripts contained here can be converted to notebooks with [jupytext](https://github.com/mwouts/jupytext). 

## Requirements

* Linux/macOS
* [conda](https://docs.conda.io/en/latest/miniconda.html) package manager
* MATLAB
* jupytext (optional)

Clone or download this repository.

Install the environment with the script.

```bash
./install_environment.sh
```

## Data access

### Mooring data (GMACMD)

http://stockage.univ-brest.fr/~scott/GMACMD/gmacmd.html

Download the data from the Google Drive. Follow the setup instructions that come with the data and create a symbolic link here, e.g.

```bash
ln -s path/to/GMACMD GMACMD
```

### Drifter analysis

Download the hourly interpolated drifter data from [NOAA](https://www.aoml.noaa.gov/phod/gdp/hourly_data.php) and convert to netcdf with the scripts.

```bash
conda activate ewglobal
cd data
./get_drifter_data.sh
python convert_drifter_mat_to_nc.py
```

### Reanalysis data

I download the GLORYS reanalysis data at 15 m depth to speed up the interpolation. (It is accessible via opendap, but this is really slow)

```bash
conda activate ewglobal
cd processing_scripts
python download_GLORYS.py
```

## Data preparation

The downloaded data is not just ready-to-use unfortunately. A series of processing steps are required.

### Drifter preparation

Interpolate the reanalysis output to the drifter location. This takes about 2 hrs and is done by

```bash
conda activate ewglobal
cd processing_scripts
python interpolate_GLORYS_to_tracks.py
```

Then the drifter data needs to be split into segments for Fourier analysis.

```bash
python create_drifter_segments.py
```

## Scientific analysis

After the data prep comes the actual science! This all takes place in notebooks in the `analysis/` directory.

## Other things

I test opendap access to various reanalysis/altimetry products in `processing_scripts/opendap_check.py`