#!/usr/bin/env bash
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate ewglobal && jupyter kernelspec uninstall ewglobal && conda deactivate
conda remove --name ewglobal --all