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
# %% [markdown]
#

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
endyear = '0005-12'

# %% [markdown]
# ## Cases

# %%
cases_sec = [
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl_smax',
    'NF1850_SECT_paas',
    'NF1850_SECT_depT',
    'NF1850_aeroxid2014_SECT_paas',
    'NF1850_aeroxid2014_SECT_depT',
]
cases_orig = [
    'NF1850_noSECT_def_smax',
    'NF1850_aeroxid2014_noSECT_def_smax',
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_noSECT_ox_ricc',
    'NF1850_aeroxid2014_noSECT_ox_ricc',
    'NF1850_noSECT_ox_ricc_depT',
    'NF1850_aeroxid2014_noSECT_ox_ricc_depT',
]

cases = cases_orig + cases_sec

# %%
version = 'pi_pd_prof_review'
plot_path = get_plotpath('prof')
filen_base = plot_path + '/%s' % version
# print(plot_path)
make_folders(plot_path)

# %%
# %%
varl = ['NCONC01', 'N_AER', 'AWNC_incld', 'AREL_incld','SOA_NA','SO4_NA']#,'N50','N100','N150','N200','N250']#, 'CDNUMC',
        #'SWCF_Ghan', 'LWCF_Ghan', 'NCFT_Ghan']
subfig_size = 2.9
asp_ratio = .9
print(varl)
# %%
from useful_scit.imps import *

area='Global'
log.ger.setLevel(log.log.DEBUG)

prof_dic = get_averaged_fields.get_profiles(cases, varl, startyear, endyear, area=area,
                                            pressure_adjust=pressure_adjust)

# %%
import pandas as pd
from pathlib import Path
case_types = ['PI', 'PD']
mod_types = ['OsloAeroSec$_{paas}$','OsloAeroSec$_{depT}$', 'OsloAero$_{depT}$','OsloAeroSec',   'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]
cdic = {key: get_case_col(key) for key in mod_types}  # , ['r','g','b'])}


# %%
def plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title=''):
    #var='NCONC01'
    #xscale='log'
    #yscale='log'
    #ylim=[1e3, 100]
    #pressure_coords=True
    #title=''
    if ylim is None:
        ylim = [1e3, 100]
    if ax is None:
        fig, ax = plt.subplots(1, figsize=[3.5,3.5], dpi=200)
    hndls=[]
    for mty in mod_types:
        case_PI= get_casen_by_type_mod('PI',mty)
        case_PIPDaer= get_casen_by_type_mod('PD',mty)
        daPI = prof_dic[case_PI][var]
        daPIPD = prof_dic[case_PIPDaer][var]
        h = plt.fill_betweenx(daPI.lev,daPI, daPIPD, color = cdic[mty],
                        label=mty,
                          alpha=.4
                              )
        hndls.append(h)
    for mty in mod_types:
        
        case_PI= get_casen_by_type_mod('PI',mty)
        case_PIPDaer= get_casen_by_type_mod('PD',mty)
        daPI = prof_dic[case_PI][var]
        daPIPD = prof_dic[case_PIPDaer][var]

        daPI.plot(  y='lev',ax=ax, color=cdic[mty])
        daPIPD.plot(y='lev',ax=ax, linestyle='dashed',color=cdic[mty],)

        #plot_profile(da, ax, xscale='log', yscale='log', label='', ylim=[1e3, 100], pressure_coords=True, kwargs={},
        #             title='')
    ax.set_ylim(ylim)
    if len(title) > 0:
        ax.set_title(title)
    if pressure_coords:
        ax.set_ylabel('Pressure [hPa]')
    xlabel = get_fancy_var_name(var) + ' [%s]' % var_info.get_fancy_unit_xr(daPI,
                                                                                                              var)
    ax.set_xlabel(xlabel)
    #ax.grid(True, which='both')
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    first_leg= plt.legend(handles=hndls, loc=1, frameon=False)
    ax_l = plt.gca().add_artist(first_leg)
    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle='dashed'),]
    plt.legend(custom_lines, ['PI', 'PD'], loc='center right', frameon=False)
    sns.despine(ax.get_figure(), ax)
    plot_profiles.set_scalar_formatter(ax)
    
    
    
    return ax

# %%
var = 'NCONC01'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'NCONC01'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'NCONC01'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N50'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N100'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N150'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N200'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N250'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'NCONC01'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
var = 'N_AER'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
                        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
# %%
var = 'AWNC_incld'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
                        title=var)
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
# %%
var = 'AREL_incld'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
                        title=var)
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
# %%


for case in cases:
    if 'aeroxid2014' in case:
        prof_dic[case]['NCONC01'].plot(y='lev', yscale='log',ylim=[1e3,100],
                                     linestyle='dashed',
                                       label=case)
    else:
        prof_dic[case]['NCONC01'].plot(y='lev', yscale='log',ylim=[1e3,100],
                                       label=case
                                       )#,linestyle='dashed')

plt.legend()
plt.show()

# %%
# %%
var = 'SOA_NA'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
                        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%

# %%
# %%
var = 'SO4_NA'
ax = plt_pipd_diff_prof(var, prof_dic, ax=None, xscale='linear', yscale='log', ylim=None, pressure_coords=True,
                        title='')
fin = f'{filen_base}_{var}'
plt.tight_layout()
plt.savefig(fin+'.pdf', dpi=300)
plt.show()

# %%
