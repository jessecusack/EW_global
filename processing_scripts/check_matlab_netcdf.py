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

# %%
import xarray as xr

# %%
ds = xr.open_dataset("../analysis/test.nc")

# %%
ds

# %%
ds.u.plot(yincrease=False)

# %%
