#!/usr/bin/env bash
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda env create -f environment.yml
conda activate ewglobal && python -m ipykernel install --user --name ewglobal