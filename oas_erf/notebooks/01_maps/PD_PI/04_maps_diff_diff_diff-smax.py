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

# %% [markdown]
# ## Differencene plots for output with Smax. 

# %% [markdown]
# ### Div imports:

# %%
# load and autoreload
from IPython import get_ipython

from oas_erf.data_info.simulation_types import get_diff_by_type
from oas_erf.util.imports import get_averaged_fields

from oas_erf.constants import get_plotpath
from oas_erf.util.practical_functions import make_folders
from matplotlib import colors
from oas_erf.data_info.simulation_types import get_abs_by_type

from oas_erf.data_info.simulation_types import get_casen_by_type_mod

from oas_erf.util.slice_average.significance import load_and_plot_sign
import cartopy.crs as ccrs

from oas_erf.util.plot.maps_PIPD import abs_diffs_PI_PD_sep, diffs_PI_PD_sep

from scipy.stats import lognorm
import matplotlib.pyplot as plt
import xarray as xr 
import numpy as np

from IPython.display import clear_output

# noinspection PyBroadException,DuplicatedCode
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass

# %% [markdown]
# ### Div settings: 

# %%
p_level = 1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013., 900., 800., 700., 600.]  # used if not avg

# %%
model = 'NorESM'

startyear = '0004-01'
endyear = '0005-12'

# %%
cases_sec = [
    'NF1850_SECT_ctrl_smax',
    'NF1850_aeroxid2014_SECT_ctrl_smax'
]
cases_orig = [
    'NF1850_noSECT_def_smax',
    'NF1850_aeroxid2014_noSECT_def_smax',
]

cases = cases_orig + cases_sec

# %% [markdown]
# ### For output names:

# %%
version = 'pi_pd_diff_smax'
plot_path = get_plotpath('maps')
filen_base = plot_path + '/%s' % version
# print(plot_path)
make_folders(plot_path)

# %% [markdown]
# ### Variables to load:

# %%
varl = ['NCONC01', 'NMR01', 'N_AER', 'ACTNL_incld', 'ACTREL_incld', 'CDNUMC', 'cb_NA',
        'cb_SOA_NA', 'cb_SO4_NA', 'AWNC_incld', 'AREL_incld', 'NACT01', 'NACT02', 'NACT04',
        'NACT06', 'NACT08', 'NACT14', 'Smax_w', 'Smax', 'SIGMA01', 'NMR01', 'NCONC04',
        'NCONC06', 'NCONC08', 'NCONC14', 'NACT_FRAC01', 'NACT_FRAC04', 'NACT_FRAC06',
        'NACT_FRAC08', 'NACT_FRAC14', 'SWCF_Ghan', 'LWCF_Ghan', 'NCFT_Ghan', 'NCONC01']
subfig_size = 2.9
asp_ratio = .9
print(varl)
# %% [markdown]
# ## Load data from file: 
# The following algortithm computes and loades map averages. 

# %%
case_dic = get_averaged_fields.get_maps_cases(cases,
                                              varl,
                                              startyear,
                                              endyear,
                                              avg_over_lev=avg_over_lev,
                                              pmin=pmin,
                                              pressure_adjust=pressure_adjust,
                                              p_level=p_level,
                                              )

# %% [markdown]
# ## Calculate various variables:

# %%
for case in cases:
    _ds = case_dic[case]
    # Fraction of particles from NPF
    _ds['NPF_frac'] = _ds['NCONC01'] / _ds['N_AER'] * 100
    _ds['NPF_frac'].attrs['units'] = '%'
    _ds['NPF_frac'].attrs['long_name'] = 'frac of particles from NPF'
    
    # noinspection DuplicatedCode
    _ds['NACT1_4'] = _ds['NACT01'] + _ds['NACT04']
    _ds['NACT1,4,6,8,14'] = _ds['NACT01'] + _ds['NACT04'] + _ds['NACT06'] + _ds['NACT08'] + _ds['NACT14']
    if _ds['Smax_w'].attrs['units'] != '%':
        _ds['Smax_w'] = _ds['Smax_w'] * 100.
        _ds['Smax_w'].attrs['units'] = '%'
        _ds['Smax_w'].attrs['long_name'] = 'S$_{max}$'
if 'NPF_frac' not in varl:
    varl.append('NPF_frac')

# %% [markdown]
# ### Get difference from PI to PD

# %%
relative = False
dic_diff = get_diff_by_type(case_dic, varl, ctrl='PI_smx', case_types=['PI_smx', 'PD_smx'],
                            mod_types=['OsloAeroSec', 'OsloAero$_{def}$'],
                            relative=relative)

dic_diff.keys()
di_dic = dic_diff['PD_smx-PI_smx']
di_dic.keys()

# %% [markdown]
# ### Organize data in easy to use format:

# %%

dic_abs = get_abs_by_type(case_dic,
                          case_types=['PI_smx', 'PD_smx'],
                          mod_types=['OsloAeroSec', 'OsloAero$_{def}$'],
                          )

# %%
dic_abs.keys()


# %%
dic_abs['PI_smx'].keys()

# %% [markdown]
# ## Set which case should act as ctr: 

# %%
ctrl = 'OsloAeroSec'

cases_oth = ['OsloAero$_{def}$']

# %% [markdown]
# ## Plots:

# %%
var = 'Smax_w'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   norm_abs=colors.LogNorm(vmin=.1, vmax=2),
                                   case_types=['PI_smx', 'PD_smx'],
                                   cases_oth=cases_oth,
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True,
                                   )

for ct in ['PI_smx', 'PD_smx']:
    ax_di = axs_dic[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,
                           ci=.95,
                           groupby=None,
                           dims=('lev',),
                           area='Global',
                           avg_dim='time',

                           hatches=['...', ''],
                           hatch_lw=.6,

                           transform=ccrs.PlateCarree(),
                           reverse=False)
clear_output()

fn = filen_base + f'{var}_PIPD_sep_rel{relative}_pmin{str(int(pmin))}.pdf'
fig.savefig(fn, dpi=300)

plt.show()

# %%
print(fn)

# %%
var = 'NCONC04'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True
                                   )
for ct in ['PI_smx', 'PD_smx']:
    ax_di = axs_dic[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,
                           ci=.95,
                           groupby=None,
                           dims=('lev',),
                           area='Global',
                           avg_dim='time',

                           hatches=['...', ''],
                           hatch_lw=.6,

                           transform=ccrs.PlateCarree(),
                           reverse=False)
clear_output()

# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NCONC01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   ctrl='OsloAero$_{def}$',
                                   cases_oth=['OsloAeroSec'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'ACTNL_incld'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'AWNC_incld'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT1,4,6,8,14'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT1,4,6,8,14'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT_FRAC04'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
for ct in ['PI_smx', 'PD_smx']:
    ax_di = axs_dic[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,
                           ci=.95,
                           groupby=None,
                           dims=('lev',),
                           area='Global',
                           avg_dim='time',

                           hatches=['...', ''],
                           hatch_lw=.6,

                           transform=ccrs.PlateCarree(),
                           reverse=False)
clear_output()

fig.show()

# %%
var = 'NACT_FRAC04'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT_FRAC01'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%

# %%

daf1 = dic_abs['PD_smx']['OsloAero$_{def}$']
a = lognorm.ppf(1 - daf1['NACT_FRAC01'], s=np.log(daf1['SIGMA01']), scale=daf1['NMR01'])
daf1['test'] = xr.DataArray(a, dims=daf1['NACT_FRAC01'].dims, coords=daf1['NACT_FRAC01'].coords)

# %%

daf = dic_abs['PD_smx']['OsloAeroSec']
a = lognorm.ppf(1 - daf['NACT_FRAC01'], s=np.log(daf['SIGMA01']), scale=daf['NMR01'])
daf['test'] = xr.DataArray(a, dims=daf['NACT_FRAC01'].dims, coords=daf['NACT_FRAC01'].coords)

# %%
(daf['test'] - daf1['test']).plot(robust=True)

# %%
((daf['test']) * 2).plot(robust=True)

# %%
((daf['test'] - daf1['test']) / daf1['test'] * 100).plot(robust=True)

# %%
(2 * daf['test']).plot(robust=True)

# %%
(2 * daf1['test']).plot(robust=True)

# %%
(2 * daf['NMR01']).plot(robust=True)

# %%
var = 'NACT_FRAC04'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT_FRAC01'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT_FRAC01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%

# %%
var = 'NACT1_4'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT01'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT04'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT04'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT06'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
var = 'NACT08'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   case_types=['PI_smx', 'PD_smx'],
                                   # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                   cases_oth=['OsloAero$_{def}$'],
                                   type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
fig.show()

# %%
for var in ['Smax_w', 'Smax', 'NACT_FRAC01', 'NACT_FRAC04', 'NACT_FRAC06', 'NACT_FRAC08', 'NACT_FRAC14']:
    relative = False
    fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                       var,
                                       relative=relative,
                                       case_types=['PI_smx', 'PD_smx'],
                                       # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                                       cases_oth=['OsloAero$_{def}$'],
                                       type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'}
                                       )
    fig.show()

# %%
# noinspection DuplicatedCode
norm_dic = dict(
    SOA_LV=colors.SymLogNorm(vmin=-5e-1, vmax=5e-1, linthresh=.01, base=10, linscale=.4),
    H2SO4=colors.SymLogNorm(vmin=-5e-1, vmax=5e-1, linthresh=.01, base=10, linscale=.4),
    NCONC01=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10, base=10, linscale=.4),
    NMR01=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=1, base=10),  # linscale=.5),
    AWNC_incld=colors.SymLogNorm(vmin=-50, vmax=50, linthresh=1, base=10),
    ACTNL_incld=colors.SymLogNorm(vmin=-40, vmax=40, linthresh=1, linscale=0.4, base=10),
    AREL_incld=colors.SymLogNorm(vmin=-5, vmax=5, linthresh=.1, base=10),
    ACTREL_incld=colors.SymLogNorm(vmin=-7, vmax=7, linthresh=.1, base=10, linscale=0.5),
    CDNUMC=None,
    SWCF_Ghan=colors.Normalize(vmin=-2, vmax=2),
    LWCF_Ghan=colors.Normalize(vmin=-2, vmax=2),
    NCFT_Ghan=colors.Normalize(vmin=-2, vmax=2),

)
norm_dic_rel = dict(
    SOA_LV=colors.Normalize(vmin=-50, vmax=50),
    H2SO4=colors.Normalize(vmin=-50, vmax=50),
    NCONC01=colors.Normalize(vmin=-250, vmax=250),
    NMR01=colors.Normalize(vmin=-50, vmax=50),
    AWNC_incld=colors.Normalize(vmin=-50, vmax=50),
    ACTNL_incld=colors.Normalize(vmin=-13, vmax=13),
    AREL_incld=colors.Normalize(vmin=-10, vmax=10),
    ACTREL_incld=colors.Normalize(vmin=-7, vmax=7),
    CDNUMC=colors.Normalize(vmin=-12, vmax=12),

    SWCF_Ghan=colors.Normalize(vmin=-2, vmax=2),
    LWCF_Ghan=colors.Normalize(vmin=-2, vmax=2),
    NCFT_Ghan=colors.Normalize(vmin=-2, vmax=2),

)

norm_abs = norm_dic.copy()
norm_abs['SWCF_Ghan'] = colors.Normalize(vmin=-5, vmax=5)
norm_abs['LWCF_Ghan'] = colors.Normalize(vmin=-3, vmax=3)
norm_abs['NCFT_Ghan'] = colors.Normalize(vmin=-5, vmax=5)

# %%
norm_diff_dic = dict(
    ACTNL_incld=colors.SymLogNorm(vmin=-50, vmax=50, linthresh=1, base=10),
    N50=colors.Normalize(vmin=-45, vmax=45),
    N100=colors.Normalize(vmin=-20, vmax=20),
    N150=colors.Normalize(vmin=-10, vmax=10),  # colors.Normalize(vmin=-5, vmax=5),
    N200=colors.Normalize(vmin=-15, vmax=15),  # colors.Normalize(vmin=-5, vmax=5),
    NACT01=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=.1, base=10),
    NACT04=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=.1, base=10),
    NACT1_4=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=.1, base=10),
    NCONC01=colors.SymLogNorm(vmin=-100, vmax=100, linthresh=10, base=10),
    NCONC04=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=1, base=10)
)

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NACT1_4', 'NACT01', 'NACT04']  # ,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = False

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

cc1 = ctrl.replace('$', '').replace('{', '').replace('}', '')
cc2 = case_oth.replace('$', '').replace('{', '').replace('}', '')
fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
fg.savefig(fn, dpi=300)
fg.show()

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NACT1_4', 'NACT01', 'NACT04']  # ,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = True

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          # norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

cc1 = ctrl.replace('$', '').replace('{', '').replace('}', '')
cc2 = case_oth.replace('$', '').replace('{', '').replace('}', '')
fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
fg.savefig(fn, dpi=300)
fig.show()

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NACT_FRAC01', 'NACT_FRAC04']  # ,'NACT04']#,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = False

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          # norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

cc1 = ctrl.replace('$', '').replace('{', '').replace('}', '')
cc2 = case_oth.replace('$', '').replace('{', '').replace('}', '')
fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
fg.savefig(fn, dpi=300)
fg.show()

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NACT_FRAC01', 'NACT_FRAC04']  # ,'NACT04']#,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = True

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          # norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

cc1 = ctrl.replace('$', '').replace('{', '').replace('}', '')
cc2 = case_oth.replace('$', '').replace('{', '').replace('}', '')
fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
fg.savefig(fn, dpi=300)
fig.show()

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NCONC01', 'NCONC04']  # ,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = True

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          # norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

# cc1 = ctrl.replace('$','').replace('{','').replace('}','')
# cc2 = case_oth.replace('$','').replace('{','').replace('}','')
# fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
# fg.savefig(fn, dpi=300)
fg.show()

# %%
varl = ['ACTNL_incld', 'N100', 'N150', 'N200']  # ,'N250'] #'N50',
varl = ['ACTNL_incld', 'NACT01', 'NACT04']  # ,'N200']#,'N250'] #'N50',

ctrl = 'OsloAero$_{def}$'
case_oth = 'OsloAeroSec'
relative = True

fg, axs = diffs_PI_PD_sep(dic_abs,
                          varl,
                          case_types=['PI_smx', 'PD_smx'],
                          ctrl=ctrl,
                          case_oth=case_oth,
                          sfg_size=3.5,
                          asp_rat=.5,
                          relative=relative,
                          norm_diff_dic=norm_diff_dic,
                          type_nndic={'PI_smx': 'Pre-industrial', 'PD_smx': 'Present day'},
                          height_ratios=None
                          )
vl = '_'.join(varl)

# cc1 = ctrl.replace('$','').replace('{','').replace('}','')
# cc2 = case_oth.replace('$','').replace('{','').replace('}','')
# fn = filen_base + f'{vl}_{cc2}_{cc1}_sep_rel{relative}.pdf'
# fg.savefig(fn, dpi=300)
fg.show()
