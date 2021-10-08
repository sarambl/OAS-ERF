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
import pandas as pd

# %%
lat_CHC = -16.35
lon_CHC = -68.13
lat_ATTO = -2.1457
lon_ATTO = - 59.0048
_dic = dict(Amazon={'lat': lat_ATTO, 'lon': lon_ATTO},
            Chacataya={'lat':lat_CHC, 'lon':lon_CHC}
           )
collocate_locations = pd.DataFrame.from_dict(_dic)
collocate_locations

# %%
from oas_erf.util.Nd.sizedist_class_v2 import SizedistributionSurface, SizedistributionStation

# %%
from oas_erf.util.plot import plot_profiles

from oas_erf.util.naming_conventions import var_info

from oas_erf.util.imports import get_averaged_fields

from IPython import get_ipython
from useful_scit.imps import (plt)
from matplotlib.lines import Line2D
import seaborn as sns

from oas_erf.data_info.simulation_types import get_diff_by_type, get_casen_by_type_mod
from oas_erf.util.imports import get_averaged_fields
from oas_erf.util.naming_conventions.var_info import get_fancy_var_name, get_fancy_unit_xr
from oas_erf.util.plot.colors import get_case_col
from oas_erf.util.plot.plot_maps import plot_map_diff, plot_map

from oas_erf.constants import get_plotpath, path_data_info
from oas_erf.util.practical_functions import make_folders
import cartopy.crs as ccrs
from matplotlib import gridspec
from matplotlib import colors

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass
# %%

# %%
p_level = 1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013., 900., 800., 700., 600.]  # used if not avg

# %%

# %%
model = 'NorESM'

startyear = '0004-01'
endyear = '0008-12'

# %%
model = 'NorESM'

startyear = '2008-01'
endyear = '2014-12'

# %% [markdown]
# ## Cases

# %%
cases_sec = [
    'SECTv21_ctrl_koagD',
    #'NF1850_SECT_ctrl',
    #'NF1850_aeroxid2014_SECT_ctrl',
]
cases_orig = [
    'noSECTv21_default_dd',
    'noSECTv21_ox_ricc_dd',

]

cases_pd = cases_orig + cases_sec

# %% [markdown]
# ## Cases

# %%
cases = cases_pd #+ cases_pi

# %%
from pathlib import Path

# %%
version = 'pd_amazon_chc'
plot_path = get_plotpath('measurement_comp')
filen_base = Path(plot_path + '/%s' % version)
filen_base.mkdir(exist_ok=True)
# print(plot_path)
make_folders(plot_path)

# %%
from oas_erf.util.slice_average.avg_pkg import yearly_mean_dic

# %%
varl = ['NCONC01','N50','N60','Z3']#,'N60']

# %%
cases

# %%
from oas_erf.util.collocate.collocate import CollocateModel

# %%
from oas_erf.data_info.variable_info import sized_varListNorESM,sized_varlist_SOA_SEC, sized_varlist_SO4_SEC

# %%
vl = sized_varListNorESM['NCONC']+sized_varListNorESM['NMR']+ sized_varListNorESM['SIGMA']

# %%
vl_sec = sized_varlist_SO4_SEC + sized_varlist_SOA_SEC

# %% [markdown]
# ## Create collocated datasets:

# %%
dic_collds = {}

# %%
for case in cases_pd:
    print(case)
    _vl = vl
    isSec = (case in cases_sec)
    if isSec:
        _vl = vl + vl_sec
    collmod = CollocateModel(case,
                             startyear, 
                             endyear,
                             isSec,
                             'month',
                             space_res='locations',
                             locations = collocate_locations
               #[5,39.6],
               #False,
              )
    collmod.locations = collocate_locations
    try:
        _ds = collmod.collocate_dataset_vars(_vl)
    except:
        collmod.load_raw_ds(_vl)
        _ds = collmod.collocate_dataset_vars(_vl)
    dic_collds[case] = _ds#.copy()
#collmod.get

# %%
import xarray as xr

# %%
dic_sdist = {}
dic_ds={}

for case in cases_pd:
    print(case)
    _vl = vl
    isSec = (case in cases_sec)
    if isSec:
        _vl = vl + vl_sec
    sdist = SizedistributionStation.SizedistributionStation(
        case,
        startyear, 
        endyear,
        [5,39.6],
        isSec,
        'month',
        locations=collocate_locations,
        
    )
    sdist.get_collocated_dataset(vl)
    
    ds = sdist.compute_sizedist_tot()
    for var in ['dNdlogD_mode01','dNdlogD_mode04','dNdlogD_mode02',
            'dNdlogD_mode05','dNdlogD_mode06','dNdlogD_mode07',
            'dNdlogD_mode08','dNdlogD_mode09','dNdlogD_mode10'
           ,'dNdlogD_mode12','dNdlogD_mode14']:
        ds[var] = sdist.compute_sizedist_mod_var(var)[var]
    dic_ds[case] = ds.copy()
    dic_sdist[case] = sdist

# %%
import xarray as xr

# %%
yscale = 'linear'

for case in cases_orig:
    ds= dic_ds[case]
    ds['dNdlogD_mod'].isel(location=1,lev=-5).mean('time').plot(xscale='log', yscale=yscale, label=case)
for case in cases_sec:
    ds= dic_ds[case]
    (ds['dNdlogD_sec']+ds['dNdlogD_mod']).isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale=yscale, label=case)
plt.legend()
plt.ylim([1e-1,5e3])

# %%
import numpy as np

# %%
yscale = 'linear'

for case in cases_orig:
    ds= dic_ds[case]
    (ds['dNdlogD_mod']*np.log(10)).isel(location=1,lev=-1).mean('time').plot(xscale='log', yscale=yscale, label=case)
for case in cases_sec:
    ds= dic_ds[case]
    ((ds['dNdlogD_sec']+ds['dNdlogD_mod'])*np.log(10)).isel(location=1,lev=-1).mean('time').plot(xscale='log', yscale=yscale, label=case)
plt.legend()
plt.ylim([1e1,8e3])
plt.xlim([1e1,5e2])

# %%

# %%
from useful_scit.imps import *

# %%
from oas_erf.constants import project_base_path, project_name
from pathlib import Path

# %%
base_dir = Path(project_base_path+project_name)

# %%
path_data = base_dir /'data/EBAS'
path_data.exists()

# %%
fl = list(path_data.glob('*.1y.*.nc'))

# %%
fl

# %% [markdown]
# fl  = path_data / 'BO0001R.*.1y.1h.BO01L_SMPS_UMSA-LAP_CHACALTAYA_BOLIVIA.DE08L_TROPOS_CL_SMPS.lev2.nc'
# #'BO0001R.20120101000000.20190917174536.smps.particle_number_size_distribution.aerosol.1y.1h.BO01L_SMPS_UMSA-LAP_CHACALTAYA_BOLIVIA.DE08L_TROPOS_CL_SMPS.lev2.nc'

# %%
ds_obs = xr.open_mfdataset(fl, combine='by_coords')

# %%
ds_obs['particle_number_size_distribution_amean']

# %%
ds_obs#['particle_number_size_distribution_amean']

# %%
logD = np.log10(ds_obs['D'])

# %%
dlogD = logD.isel(D=slice(1,None)).values -logD.isel(D=slice(0,-1))

# %%
dlogD

# %%
N = (ds_obs['particle_number_size_distribution_amean']*dlogD).sum('D')#.plot()

# %%
N.plot()

# %% jupyter={"outputs_hidden": true}
ds_obs

# %%
ds_obs['particle_number_size_distribution_amean'].isel(D=10,metadata_time=2 ).plot()

# %%
ds_obs['particle_number_size_distribution_amean'].isel(D=10).plot()

# %%
N.isel(metadata_time=0)#.plot()

# %%
ds_obs['particle_number_size_distribution_amean'].mean('metadata_time').mean('time').plot(xscale='log')

ds_obs['particle_number_size_distribution_perc8413'].mean('metadata_time').mean('time').plot(xscale='log')

ds_obs['particle_number_size_distribution_prec1587'].mean('metadata_time').mean('time').plot(xscale='log')

# %%
from oas_erf.data_info import get_nice_name_case

# %%
import numpy as np

# %%
np.log(10)

# %%
yscale = 'linear'

cndic = dict(noSECTv21_default_dd='OsloAero$_{def}$',
            noSECTv21_ox_ricc_dd = 'OsloAero$_{imp}$')

for case in cases_orig:
    ds= dic_ds[case]
    (np.log(10)*ds['dNdlogD_mod']).isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale=yscale, label=cndic[case], c = get_case_col(cndic[case]))
for case in cases_sec:
    ds= dic_ds[case]
    (np.log(10)*(ds['dNdlogD_sec']+ds['dNdlogD_mod'])).isel(location=0,lev=-1).mean('time').plot(xscale='log', c = get_case_col('OsloAeroSec'),
                                                                                    yscale=yscale, label='OsloAeroSec')
plt.ylim([5,9e3])
plt.xlim([10,500])
#plt.yscale('log')
#ds_obs['sizedist'].mean('Date').plot(xscale='log', c='k', label='OBS: ATTO tower')#x='diameter',robust = True)
plt.legend()


# %%
yscale = 'linear'

cndic = dict(noSECTv21_default_dd='OsloAero$_{def}$',
            noSECTv21_ox_ricc_dd = 'OsloAero$_{imp}$')

for case in cases_orig:
    ds= dic_ds[case]
    (np.log(10)*ds['dNdlogD_mod']).isel(location=0,lev=-1).sel(time = slice('2014-01-01','2015-01-01')).mean('time').plot(xscale='log', yscale=yscale, label=cndic[case], c = get_case_col(cndic[case]))
for case in cases_sec:
    ds= dic_ds[case]
    (np.log(10)*(ds['dNdlogD_sec']+ds['dNdlogD_mod'])).isel(location=0,lev=-1).sel(time = slice('2014-01-01','2015-01-01')).mean('time').plot(xscale='log', c = get_case_col('OsloAeroSec'),
                                                                                    yscale=yscale, label='OsloAeroSec')
plt.ylim([5,9e3])
plt.xlim([10,500])
#plt.yscale('log')
#ds_obs['sizedist'].sel(Date = slice('2014-01-01','2015-01-01')).mean('Date').plot(xscale='log', c='k', label='OBS: ATTO tower')#x='diameter',robust = True)
plt.legend()


# %%
ds['dNdlogD_mod']


# %%
ds['dNdlogD_mod'].sel(location='Chacataya',lev=slice(970,1000)).mean('lev').mean('time')

# %%
for c in cases:
    ds= dic_ds[c]

    ds['dNdlogD_mod'].sel(location='Chacataya',lev=slice(970,1000)).mean('lev').mean('time').plot()

    ds['dNdlogD_mod'].sel(location='Chacataya',).isel(lev=-1).mean('time').plot()
    ds['dNdlogD_mod'].sel(location='Chacataya',).isel(lev=-2).mean('time').plot()
    ds['dNdlogD_mod'].sel(location='Chacataya',).isel(lev=-3).mean('time').plot()
    plt.xscale('log')
    plt.xlim([10,700])
    plt.ylim([10,8000])
    plt.ylim([10,8000])
    plt.yscale('log')
    plt.show()


# %%
def mean_for_comp(ds, isSec=False, vmon = 'dNdlogD_mod', vsec='dNdlogD_sec', p0=970,p1=1000):
    log10 = np.log(10)
    _da = ds[vmon]
    if isSec:
        _da = _da + ds[vsec]
    da_dlog10 = log10*_da
    da_sel = da_dlog10.sel(location='Chacataya',lev = slice(p0,p1)).mean('lev')
    da_mean = da_sel.mean('time')
    return da_mean


# %%
yscale = 'linear'
p0=650
p1=900
cndic = dict(noSECTv21_default_dd='OsloAero$_{def}$',
            noSECTv21_ox_ricc_dd = 'OsloAero$_{imp}$')
f, ax = plt.subplots(figsize=[4.5,3], dpi=100)
for case in cases_orig:
    ds= dic_ds[case]
    da = mean_for_comp(ds, p0=p0,p1=p1)
    da.plot(xscale='log', yscale=yscale, label=cndic[case], c = get_case_col(cndic[case]))
for case in cases_sec:
    ds= dic_ds[case]
    da = mean_for_comp(ds, isSec=True, p0=p0,p1=p1)
    
    da.plot(xscale='log', c = get_case_col('OsloAeroSec'),yscale=yscale, label='OsloAeroSec')
plt.ylim([20,16e3])
plt.xlim([8,600])
plt.yscale('log')
da_obs = ds_obs['particle_number_size_distribution_amean'].mean('metadata_time')
da_obs.mean('time').plot(xscale='log', c='k', label='Chacataya')#x='diameter',robust = True)
#da_obs.groupby(da_obs.time.dt.month).mean().mean('month').plot(xscale='log', c='k', label='ATTO Tower')#x='diameter',robust = True)

plt.legend(frameon=False)
sns.despine(f)
ax.set_ylabel('dN/dlog$_{10}$D$_p$ [cm$^{-3}$]')
ax.set_xlabel('diameter [nm]')
f.tight_layout()

fn = filen_base /'CHC_sizedist_mean2'
f.savefig(fn.with_suffix('.png'))
f.savefig(fn.with_suffix('.pdf'))



# %%
ds_obs


# %%
def mean_for_comp_month(ds, isSec=False, vmon = 'dNdlogD_mod', vsec='dNdlogD_sec', p0=650, p1=900):
    log10 = np.log(10)
    _da = ds[vmon]
    if isSec:
        _da = _da + ds[vsec]
    da_dlog10 = log10*_da
    da_sel = da_dlog10.isel(location=1).sel(lev = slice(p0,p1)).mean('lev')
    da_mean = da_sel.groupby(da_sel.time.dt.month).mean()
    return da_mean


# %%
from matplotlib.colors import LogNorm

# %%

p0=650
p1=900
yscale = 'log'
xscale='linear'
f,axs = plt.subplots(1,4, figsize=[9 ,4], sharex=True, sharey=True, dpi=100)
cndic = dict(noSECTv21_default_dd='OsloAero$_{def}$',
            noSECTv21_ox_ricc_dd = 'OsloAero$_{imp}$',
            SECTv21_ctrl_koagD = 'OsloAeroSec')
norm = LogNorm(vmin=50, vmax=12000)

plt_sett = dict(
    norm=norm,
    cmap='cividis',
    yscale=yscale,
    xscale=xscale,
    ylim=[10,500],
    add_colorbar=False,
    rasterized=True
)


for case,ax in zip(cases_orig,axs):
    ds= dic_ds[case]
    da = mean_for_comp_month(ds, p0=p0, p1=p1)
    da.plot(x='month',**plt_sett,ax=ax)#, label=cndic[case], c = get_case_col(cndic[case]))
    ax.set_title(cndic[case])
for case, ax in zip(cases_sec, axs[len(cases_orig):]):
    ds= dic_ds[case]
    da = mean_for_comp_month(ds, isSec=True, p0=p0, p1=p1)
    da.plot(x='month',**plt_sett,ax=ax)#, label=cndic[case], c = get_case_col(cndic[case]))
    ax.set_title(cndic[case])

ax = axs[-1]
da_obs = ds_obs['particle_number_size_distribution_amean'].mean('metadata_time')
#da_obs.mean('time').plot(xscale='log', c='k', label='ATTO Tower')#x='diameter',robust = True)
im = da_obs.groupby(da_obs.time.dt.month).mean().plot(x='month',**plt_sett, ax = ax)#x='diameter',robust = True)


for ax in axs[1:].flatten():
    ax.set_ylabel('')
#for ax in axs_ba[-1,:].flatten():
#    ax.set_xlabel('Latitude [$^\circ$ N]')
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.83, 0.25, 0.015, 0.5])

cb = f.colorbar(im, cax=cbar_ax, extend = 'both', label= 'dN/dlog$_{10}$D$_p$ [cm$^{-3}$]' )
ax.set_title('Chacataya')

fn = filen_base /'CHC_sizedist_month'
f.savefig(fn.with_suffix('.png'), bbox_extra_artists=(cb,))
f.savefig(fn.with_suffix('.pdf'), bbox_extra_artists=(cb,))
plt.show()


# %% [markdown]
# ## Extra: 

# %%
for case in cases_pd:
    ds = dic_ds[case]
    ds['dNdlogD_mod'].isel(location=1,lev=-1).mean('time').plot(xscale='log', yscale='log')

    for var in [
        'dNdlogD_mode01',
        'dNdlogD_mode02',
        'dNdlogD_mode04',
        'dNdlogD_mode05',
        'dNdlogD_mode06',
        'dNdlogD_mode07',
        'dNdlogD_mode08',
        'dNdlogD_mode09',
        'dNdlogD_mode10',
        'dNdlogD_mode12',
        'dNdlogD_mode14'
    ]:
        ds[var].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log', label=var)
#ds['dNdlogD_mode02'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
#ds['dNdlogD_mode04'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
    plt.legend()
    #plt.yscale('linear')
    plt.ylim([10e-1,1e4])
    plt.title(case)
    plt.show()

# %%
for case in cases_pd:
    ds = dic_ds[case]
    #ds['dNdlogD_mod'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')

    for var in [#'dNdlogD_mode01',
                'dNdlogD_mode02','dNdlogD_mode04',
            #'dNdlogD_mode05','dNdlogD_mode06','dNdlogD_mode07',
            #'dNdlogD_mode08',
            #    'dNdlogD_mode09','dNdlogD_mode10'
           #'dNdlogD_mode12','dNdlogD_mode14'
    ]:
        ds[var].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log', label=var)
#ds['dNdlogD_mode02'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
#ds['dNdlogD_mode04'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
    plt.legend()
    #plt.yscale('linear')
    plt.ylim([10e-1,1e4])
    plt.title(case)
    #plt.show()

# %%
ds['dNdlogD_mod'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')

for var in ['dNdlogD_mode01','dNdlogD_mode04','dNdlogD_mode02',
            'dNdlogD_mode05','dNdlogD_mode06','dNdlogD_mode07',
            'dNdlogD_mode08','dNdlogD_mode09','dNdlogD_mode10'
           ,'dNdlogD_mode12','dNdlogD_mode14']:
    ds[var].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log', label=var)
#ds['dNdlogD_mode02'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
#ds['dNdlogD_mode04'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
plt.legend()
plt.yscale('linear')
plt.ylim([10e-1,7e2])

# %%
ds['dNdlogD_mod'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
ds['dNdlogD_mode01'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
ds['dNdlogD_mode02'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
ds['dNdlogD_mode04'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
plt.ylim([10e-1,10e4])

# %%
ds['dNdlogD_mod'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')
ds['dNdlogD_mode01'].isel(location=0,lev=-1).mean('time').plot(xscale='log', yscale='log')

# %%
f, ax = plt.subplots(dpi=150)
(ds['dNdlogD_mod']+ds['dNdlogD_sec']).sel(location='Chacataya').isel(lev=-1).plot(yscale='log', x='time', robust=True)

# %%
