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
# ## Levlat plots with abs, diff for PI and PD separately

# %% [markdown]
# ### div imports

# %%
import matplotlib.colors as colors
# load and autoreload
from IPython import get_ipython
import matplotlib.pyplot as plt

from oas_erf.constants import get_plotpath
from oas_erf.data_info.simulation_types import get_diff_by_type
from oas_erf.util.imports import get_averaged_fields
from oas_erf.util.plot.plot_levlat import plot_levlat_abs, get_cbar_label
from oas_erf.util.plot.plot_levlat import plot_levlat_diff
from oas_erf.util.practical_functions import make_folders
from oas_erf.util.plot.levlat_PIPD import abs_diffs_PI_PD_sep
from oas_erf.data_info.simulation_types import get_abs_by_type
from IPython.display import clear_output

from oas_erf.data_info.simulation_types import get_casen_by_type_mod

from oas_erf.util.slice_average.significance import load_and_plot_sign

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass

# %%

# %% [markdown]
# ## Div settings:

# %%
model = 'NorESM'

startyear = '0004-01'
endyear = '0008-12'
p_level = 1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013., 900., 800., 700., 600.]  # used if not avg

# %%
cases_sec = [
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl'
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
    SOA_LV=colors.SymLogNorm(vmin=-1e-1, vmax=1e-1, linthresh=.01, base=10),
    H2SO4=colors.SymLogNorm(vmin=-1e-1, vmax=1e-1, linthresh=.01, base=10),
    NCONC01=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10, base=10),
    N_AER=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10, base=10),
    N=colors.SymLogNorm(vmin=-2e2, vmax=2e2, linthresh=1, base=10),
    NMR01=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=1, base=10),  # linscale=.5),
    AWNC_incld=colors.SymLogNorm(vmin=-50, vmax=50, linthresh=1, base=10),
    AREL_incld=colors.SymLogNorm(vmin=-5, vmax=5, linthresh=.1, base=10)
)


# %% [markdown]
# ## Plot function:

# %%
def abs_diffs(di_dic, ctrl, cases_oth, varl):
    fig, axs = plt.subplots(4, len(varl),
                            gridspec_kw={'height_ratios': [4, 3, 3, .3]},
                            figsize=[subfig_size * len(varl), subfig_size * 3 * asp_ratio])
    axs_diff = axs[1:-1, :]
    axs_diff_cb = axs[-1, :]
    # ctrl = 'OsloAeroSec'
    # cases_oth = ['OsloAero$_{imp}$','OsloAero$_{def}$']
    for i, var in enumerate(varl):
        print(i, var)
        saxs = axs_diff[:, i]
        ax = saxs[0]
        for case_oth, ax in zip(cases_oth, saxs.flatten()):
            _, im = plot_levlat_diff(var, ctrl, case_oth,
                                     di_dic,
                                     cbar_orientation='horizontal',
                                     # title=None,
                                     ax=ax,
                                     # ylim=None,
                                     # figsize=None,
                                     cmap='RdBu_r',
                                     # use_ds_units=True,
                                     # add_colorbar=True,
                                     norm=norm_dic[var],
                                     add_colorbar=False
                                     )

        # ax.set_title(f'{key}: PIaerPD-PI')
        lab = f'$\Delta${get_cbar_label(di_dic[ctrl][var], var, diff=True)}'
        plt.colorbar(im, cax=axs_diff_cb[i], label=lab, orientation='horizontal')

    for i, var in enumerate(varl):
        print(i, var)
        ax = axs[0, i]
        _, im = plot_levlat_abs(var, ctrl,
                                di_dic,
                                cbar_orientation='horizontal',
                                # title=None,
                                ax=ax,
                                # ylim=None,
                                # figsize=None,
                                cmap='PuOr_r',
                                # use_ds_units=True,
                                # add_colorbar=True,
                                norm=norm_dic[var],
                                add_colorbar=False
                                )

        # ax.set_title(f'{key}: PIaerPD-PI')
        lab = get_cbar_label(di_dic[ctrl][var], var, diff=True)
        plt.colorbar(im, ax=ax, label=lab, orientation='horizontal')

    for i in range(len(axs_diff[:, 0])):
        for j in range(len(axs_diff[0, :])):
            ax = axs_diff[i, j]
            if i < (len(axs_diff[:, 0]) - 1):
                ax.set_xlabel('')
                plt.setp(ax.get_xticklabels(), visible=False)
            if i == (len(axs_diff[:, 0]) - 1):
                ax.set_xlabel('Latitude [$^\circ$N]')
            if j > 0:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_ylabel('')
    for ax in axs[0, :]:
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axs[0, 1:]:
        ax.set_ylabel('')
        plt.setp(ax.get_yticklabels(), visible=False)

    fig.tight_layout()
    return axs, fig


# %% [markdown]
# ## Savepaths:

# %%
version = 'diffs'
plot_path = get_plotpath('levlat')
filen_base = plot_path + '/_%s' % version
# print(plot_path)
make_folders(plot_path)
# %% [markdown]
# ## Variables to load: 

# %%
varl = ['NCONC01', 'AWNC_incld', 'AREL_incld', 'HYGRO01', 'CLDLIQ', 'N100', 'N250', 'FREQL', 'N50', 'N150', 'N200',
        'N_AER']
subfig_size = 2.9
asp_ratio = .9
print(varl)
# %% [markdown]
# ## Load data:

# %%
case_dic = get_averaged_fields.get_levlat_cases(cases, varl, startyear, endyear,
                                                pressure_adjust=pressure_adjust)

# %% [markdown]
# ## Get PI to PD difference

# %%
relative = False
dic_diff = get_diff_by_type(case_dic, varl, ctrl='PI', case_types=['PI', 'PIaerPD'],
                            # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                            relative=relative)

# %% [markdown]
# ### Organize data according to PI/PD and model version:

# %%

dic_abs = get_abs_by_type(case_dic,
                          case_types=['PI', 'PD'],
                          mod_types=None)

# %%

# %%
cases_oth = ['OsloAero$_{imp}$', 'OsloAero$_{def}$']
ctrl = 'OsloAeroSec'

# %% [markdown]
# ## Plots:

# %%

var = 'NCONC01'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=norm_dic[var],
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        t, T = load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                                  avg_over_lev=avg_over_lev,

                                  ci=.95,
                                  groupby=None,
                                  dims=('lon',),
                                  area='Global',
                                  avg_dim='time',
                                  hatches=['...', ''], hatch_lw=.3,
                                  transform=None,
                                  reverse=False)
clear_output()
# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'N_AER'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=norm_dic[var],
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%

var = 'AWNC_incld'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=norm_dic[var],
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'N50'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=norm_dic['NCONC01'],
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%

var = 'N100'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=norm_dic['N'],
                                   switch_diff=True,
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'N150'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=colors.SymLogNorm(vmin=-7e1, vmax=7e1, linthresh=1),
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%

var = 'N200'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=colors.SymLogNorm(vmin=-1e1, vmax=1e1, linthresh=1),
                                   switch_diff=True
                                   # norm_diff=colors.Normalize(vmin=-5, vmax=5)#, linthresh=1)

                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'N250'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=colors.SymLogNorm(vmin=-7, vmax=7, linthresh=1),
                                   switch_diff=True
                                   # norm_dic['N']
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%

var = 'FREQL'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=False,
                                   norm_diff=colors.Normalize(vmin=-.005, vmax=.005),
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'CLDLIQ'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=True,
                                   norm_diff=colors.Normalize(vmin=-50, vmax=50),
                                   switch_diff=True
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

# plt.tight_layout()
fn = filen_base + f'{var}_PI_diff_PD_diff'
plt.savefig(fn + '.pdf', dpi=300)
print(fn)
plt.show()

# %%
var = 'HYGRO01'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=True,

                                   switch_diff=True
                                   # norm_diff=norm_dic[var]
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

plt.show()
# %%

var = 'NCONC01'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=True,
                                   switch_diff=True
                                   # norm_diff=norm_dic[var]
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

plt.show()
# %%

var = 'CLDLIQ'
fg, axs_dict = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   case_types=None,
                                   ctrl=None,
                                   cases_oth=None,
                                   sfg_size=2.9,
                                   asp_rat=.9,
                                   relative=True,
                                   switch_diff=True
                                   # norm_diff=norm_dic[var]
                                   )
for ct in ['PI', 'PD']:
    ax_di = axs_dict[ct]
    for case_oth in cases_oth:
        ax = ax_di[case_oth]
        cs_to = get_casen_by_type_mod(ct, ctrl)
        cs_from = get_casen_by_type_mod(ct, case_oth)

        load_and_plot_sign(cs_to, [cs_from], [ax], var, startyear, endyear, pressure_adjust=pressure_adjust,
                           avg_over_lev=avg_over_lev,

                           ci=.95,
                           groupby=None,
                           dims=('lon',),
                           area='Global',
                           avg_dim='time',
                           hatches=['...', ''], hatch_lw=.3,
                           transform=None,
                           reverse=False)
clear_output()

plt.show()
