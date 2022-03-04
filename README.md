# Eddy-Internal Wave Interactions From Historical Moorings

Internal wave stress estimated from moorings are coupled with mesoscale eddy strain from altimetry. 

## Requirements

* Linux/macOS
* conda package manager
* MATLAB

Clone or download this repository. 

### Moorings - GMACMD data access

Download the data from the Google Drive. Follow the setup instructions and link here, e.g.

```bash
ln -s path/to/GMACMD GMACMD
```

### Altimetry -  Copernicus Marine data access

Create an account: https://marine.copernicus.eu/

## How to run

1) Install the conda environment

```bash
./install_environment.sh
```

2) Run some preliminary mooring processing scripts. 