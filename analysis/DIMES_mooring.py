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
#     display_name: dimes-eiw
#     language: python
#     name: dimes-eiw
# ---

# %%
import numpy as np
from scipy.io import loadmat

# %%
dat = loadmat("../data/external/DIMES_mooring.mat", squeeze_me=True)

# %%
dat.keys()

# %%
c = dat["c"]

# %%
c["Temp"].item()

# %%
