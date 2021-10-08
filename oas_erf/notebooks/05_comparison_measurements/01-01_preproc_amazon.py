# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
from useful_scit.imps import (pd,np, xr,)


# %% [markdown]
# ## Data from ATTO tower level 3 (should be 60 m)
#
#

# %% [markdown]
# Downloaded from http://ftp.lfa.if.usp.br/ftp/public/LFA_Processed_Data/T0a_ATTO/Level3/SMPS_2014toNov2020_ATTO_60m_InstTower/

# %%
# !wget -O tmp_data/ATTO-SMPS-clean_stp-transm-2014-202009_60m_correction_2021.dat http://ftp.lfa.if.usp.br/ftp/public/LFA_Processed_Data/T0a_ATTO/Level3/SMPS_2014toNov2020_ATTO_60m_InstTower/ATTO-SMPS-clean_stp-transm-2014-202009_60m_correction_2021.dat

# %% [markdown]
# Need to 

# %%
from pathlib import Path

import pandas as pd

# %%
fn = Path('tmp_data/ATTO-SMPS-clean_stp-transm-2014-202009_60m_correction_2021.dat')

# %%
df = pd.read_table(fn).set_index('Date').to_csv(fn.with_suffix('.csv'))

# %%
fn = Path('tmp_data/ATTO-SMPS-clean_stp-transm-2014-202009_60m_correction_2021.csv')

# %%
df = pd.read_csv(fn, index_col=0)#, delimiter=';', index_col=0)
df

# %% [markdown]
# Make index datetime

# %%
df.index = pd.to_datetime(df.index)


# %%
df_h = df.resample('1h').mean()
df_h

# %%
df_h = df_h.drop(['doy','dec.t'], axis=1)

# %%
df_h.resample('1D').mean()['diam385.4'].plot()

# %%
df_h.resample('1M').mean()['diam385.4'].plot()

# %%
_diam = df_h.columns

# %%
diameters = [float(d[4:]) for d in _diam]
diameters
diameters_dic = {d:float(d[4:]) for d in _diam}

# %%
df_hc = df_h.rename(diameters_dic, axis=1)

# %%
df_hc.index

# %%
df_hc.to_xarray()

# %%
da = df_hc.to_xarray().to_array(dim='diameter')


# %%
ds = da.to_dataset(name='sizedist')

# %%
ds['Date']=df_hc.index.values#.#.plot()

# %%
ds['Date']

# %%
ds['sizedist'].plot(robust = True, yscale='log')

# %%
ds

# %%
ds.to_netcdf('tmp_data/amazon_ds.nc')

# %%
ds['sizedist'].median('Date').plot(xscale='log')#x='diameter',robust = True)

# %%
ds['sizedist'].mean('Date').plot(xscale='log')#x='diameter',robust = True)
