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
import matplotlib.colors as colors
# load and autoreload
from IPython import get_ipython
from oas_erf.util.plot.plot_levlat import plot_levlat_abs, get_cbar_label
from useful_scit.imps import (np, plt, sns)
import numpy as np
from oas_erf.data_info import get_nice_name_case
from oas_erf.util.imports import get_averaged_fields
from oas_erf.util.imports.get_fld_fixed import get_field_fixed
from oas_erf.util.plot.plot_levlat import plot_levlat_diff, get_cbar_eq_kwargs, make_cbar_kwargs
from oas_erf.constants import get_plotpath
from oas_erf.util.practical_functions import make_folders
from oas_erf.util.naming_conventions.var_info import get_fancy_var_name, get_fancy_unit_xr

from oas_erf.data_info.simulation_types import get_diff_by_type, get_abs_by_type

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass

# %%
from oas_erf.util.slice_average.significance import load_and_plot_sign

from useful_scit.plot.fig_manip import subp_insert_abc

# %%
model = 'NorESM'

startyear = '0004-01'
endyear = '0005-12'
p_level = 1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013., 900., 800., 700., 600.]  # used if not avg

# %%
cases_sec = [
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_svoc_smax',
    'NF1850_aeroxid2014_SECT_svoc_smax',
    'NF1850_SECT_elvoc_smax',
    'NF1850_aeroxid2014_SECT_elvoc_smax'

    ]
cases_orig = [
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_aeroxid2014_noSECT_ox_ricc',
    'NF1850_noSECT_ox_ricc'
    ]

cases = cases_orig + cases_sec

# %%
norm_dic = dict(
    SOA_LV = colors.SymLogNorm(vmin=-1e-1, vmax=1e-1, linthresh=.01),
    H2SO4 = colors.SymLogNorm(vmin=-1e-1, vmax=1e-1, linthresh=.01),
    NCONC01=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10),
    N_AER=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10),
    N=colors.SymLogNorm(vmin=-2e2, vmax=2e2, linthresh=1),
    NMR01=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=1),# linscale=.5),
    AWNC_incld=colors.SymLogNorm(vmin=-50, vmax=50, linthresh=1),
    AREL_incld=colors.SymLogNorm(vmin=-5, vmax=5, linthresh=.1)
)


# %%
def abs_diffs(di_dic, ctrl, cases_oth, varl):
    fig, axs = plt.subplots(4,len(varl),
                        gridspec_kw={'height_ratios': [4,3,3, .3]},
                        figsize=[subfig_size*len(varl),subfig_size*3*asp_ratio])
    axs_diff = axs[1:-1,:]
    axs_diff_cb = axs[-1,:]
    #ctrl = 'OsloAeroSec'
    #cases_oth = ['OsloAero$_{imp}$','OsloAero$_{def}$']
    for i, var in enumerate(varl):
        print(i,var)
        saxs = axs_diff[:,i]
        ax = saxs[0]
        for case_oth, ax in zip(cases_oth, saxs.flatten()):
            _, im = plot_levlat_diff(var, ctrl,case_oth,
                                     di_dic,
                            cbar_orientation='horizontal',
                            #title=None,
                            ax=ax,
                            #ylim=None,
                            #figsize=None,
                            cmap='RdBu_r',
                            #use_ds_units=True,
                            #add_colorbar=True,
                            norm = norm_dic[var],
                            add_colorbar=False
                            )




        #ax.set_title(f'{key}: PIaerPD-PI')
        lab = f'$\Delta${get_cbar_label(di_dic[ctrl][var], var, diff=True)}'
        plt.colorbar(im, cax = axs_diff_cb[i],label=lab,  orientation='horizontal')

    for i, var in enumerate(varl):
        print(i,var)
        ax = axs[0,i]
        _, im = plot_levlat_abs(var, ctrl,
                                     di_dic,
                                     cbar_orientation='horizontal',
                                     #title=None,
                                     ax=ax,
                                     #ylim=None,
                                     #figsize=None,
                                     cmap='PuOr_r',
                                     #use_ds_units=True,
                                     #add_colorbar=True,
                                     norm = norm_dic[var],
                                     add_colorbar=False
                                     )

        #ax.set_title(f'{key}: PIaerPD-PI')
        lab = get_cbar_label(di_dic[ctrl][var], var, diff=True)
        plt.colorbar(im, ax = ax,label=lab,  orientation='horizontal')

    for i in range(len(axs_diff[:,0])):
        for j in range(len(axs_diff[0,:])):
            ax = axs_diff[i,j]
            if i<(len(axs_diff[:,0])-1):
                ax.set_xlabel('')
                plt.setp(ax.get_xticklabels(), visible=False)
            if i ==(len(axs_diff[:,0])-1):
                ax.set_xlabel('Latitude [$^\circ$N]')
            if j>0:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_ylabel('')
    for ax in axs[0,:]:
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axs[0,1:]:
        ax.set_ylabel('')
        plt.setp(ax.get_yticklabels(), visible=False)

    fig.tight_layout()
    return axs, fig

# %%
# %%
version = 'diffs_svoc'
plot_path = get_plotpath('levlat')
filen_base = plot_path + '/_%s' % version
# print(plot_path)
make_folders(plot_path)
pressure_adjust=True

varl = ['NCONC01', 'AWNC_incld', 'AREL_incld','HYGRO01','CLDLIQ', 'FREQL','N_AER']# 'N100','N250',
subfig_size = 2.9
asp_ratio = .9
print(varl)
# %%
case_dic = get_averaged_fields.get_levlat_cases(cases, varl, startyear, endyear,
                                                pressure_adjust=pressure_adjust)

# %%
relative=False
dic_diff = get_diff_by_type(case_dic, varl, ctrl ='PI', case_types=['PI', 'PIaerPD'],  #mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                            relative=relative)

# %%
from oas_erf.util.plot.levlat_PIPD import abs_diffs_PI_PD_sep
from oas_erf.data_info.simulation_types import get_abs_by_type
dic_abs =get_abs_by_type(case_dic,
                    case_types=['PI', 'PD'],
                    mod_types=['OsloAeroSec','OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$','OsloAero$_{imp}$', 'OsloAero$_{def}$'])

# %%
var='NCONC01'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                        var,
                        case_types=None,
                        ctrl='OsloAeroSec',
                        cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$'],#,'OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                        sfg_size=2.9,
                        asp_rat=.9,
                        relative=False,
                    norm_diff=colors.SymLogNorm(vmin=-5e2, vmax=5e2, linthresh=10),#norm_dic[var], 
                    switch_diff=True
                   )
#plt.tight_layout()
fn = filen_base +f'{var}_PI_diff_PD_diff'
plt.savefig(fn+'.pdf', dpi=300)
print(fn)
plt.show()

# %%
var='N_AER'
abs_diffs_PI_PD_sep(dic_abs,
                    var,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$'],#'OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                    sfg_size=2.9,
                    asp_rat=.9,
                    relative=False,
                    norm_diff=norm_dic[var],
                    switch_diff=True
                   )
#plt.tight_layout()
fn = filen_base +f'{var}_PI_diff_PD_diff'
plt.savefig(fn+'.pdf', dpi=300)
print(fn)
plt.show()

# %%

var='AWNC_incld'
    

abs_diffs_PI_PD_sep(dic_abs,
                    var,
                    case_types=None,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$'],#'OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                    sfg_size=2.9,
                    asp_rat=.9,
                    relative=False,
                    norm_diff=colors.SymLogNorm(vmin=-20, vmax=20, linthresh=1),#norm_dic[var],
                    switch_diff=True
                   )
#plt.tight_layout()
fn = filen_base +f'{var}_PI_diff_PD_diff'
plt.savefig(fn+'.pdf', dpi=300)
print(fn)
plt.show()

# %%

var='FREQL'
abs_diffs_PI_PD_sep(dic_abs,
                        var,
                        case_types=None,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$','OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                        sfg_size=2.9,
                        asp_rat=.9,
                        relative=False,
                    norm_diff=colors.Normalize(vmin=-.005,vmax=.005),
                    switch_diff=True
                   )
#plt.tight_layout()
fn = filen_base +f'{var}_PI_diff_PD_diff'
plt.savefig(fn+'.pdf', dpi=300)
print(fn)
plt.show()

# %%
var='CLDLIQ'
abs_diffs_PI_PD_sep(dic_abs,
                        var,
                        case_types=None,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$','OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                        sfg_size=2.9,
                        asp_rat=.9,
                        relative=True,
                    norm_diff=colors.Normalize(vmin=-50,vmax=50),
                    switch_diff=True
                   )
#plt.tight_layout()
fn = filen_base +f'{var}_PI_diff_PD_diff'
plt.savefig(fn+'.pdf', dpi=300)
print(fn)
plt.show()

# %%
var='HYGRO01'
abs_diffs_PI_PD_sep(dic_abs,
                        var,
                        case_types=None,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$'],#'OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                        sfg_size=2.9,
                        asp_rat=.9,
                        relative=True,
                    
                    switch_diff=True
                    #norm_diff=norm_dic[var]
                   )
plt.show()
# %%

var='NCONC01'
abs_diffs_PI_PD_sep(dic_abs,
                    var,
                    case_types=None,
                    ctrl='OsloAeroSec',
                    cases_oth=['OsloAeroSec$_{svoc}$','OsloAeroSec$_{elvoc}$','OsloAero$_{imp}$', 'OsloAero$_{def}$'],
                    sfg_size=2.9,
                    asp_rat=.9,
                    relative=True,
                    switch_diff=True,
                    norm_diff=colors.Normalize(vmin=-60,vmax=60)
                   )
plt.show()
# %%

var='CLDLIQ'
abs_diffs_PI_PD_sep(dic_abs,
                        var,
                        case_types=None,
                        ctrl=None,
                        cases_oth=None,
                        sfg_size=2.9,
                        asp_rat=.9,
                        relative=True,
                    switch_diff=True                    
                    #norm_diff=norm_dic[var]
                   )
plt.show()
