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
# ## PI to PD difference plots

# %% [markdown]
# ### Div imports:

# %%
# load and autoreload
from IPython import get_ipython
from matplotlib import colors

import matplotlib.pyplot as plt

from oas_erf.constants import get_plotpath
from oas_erf.data_info.simulation_types import get_diff_by_type
from oas_erf.util.imports import get_averaged_fields
from oas_erf.util.practical_functions import make_folders
from oas_erf.data_info.simulation_types import get_abs_by_type
from oas_erf.util.slice_average.avg_pkg import yearly_mean_dic
import cartopy.crs as ccrs

from IPython.display import clear_output

from oas_erf.data_info.simulation_types import get_casen_by_type_mod

from oas_erf.util.slice_average.significance import calc_significance_field, hatch_area_sign

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except ImportError:
    pass

from oas_erf.util.plot.maps_PIPD import abs_diffs

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
endyear = '0008-12'

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

# %% [markdown]
# ### For output names:

# %%
version = 'pi_pd_diff'
plot_path = get_plotpath('maps')
filen_base = plot_path + '/%s' % version
# print(plot_path)
make_folders(plot_path)

# %%
print(filen_base)

# %% [markdown]
# ### Variables to load

# %%
varl = ['NCONC01', 'NMR01', 'N_AER', 'NCONC08', 'NCONC09', 'NMR08', 'NMR09', 'ACTNL_incld', 'ACTREL_incld', 'CDNUMC',
        'cb_NA',
        'cb_SOA_NA', 'cb_SO4_NA', 'AWNC_incld', 'AREL_incld',
        'SWCF_Ghan', 'LWCF_Ghan', 'NCFT_Ghan', 'N50', 'N100', 'N250', 'N150', 'N200', 'SIGMA01', 'NMR01', 'NCONC01']
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
    _ds['NPF_frac'] = _ds['NCONC01'] / _ds['N_AER'] * 100
    _ds['NPF_frac'].attrs['units'] = '%'
if 'NPF_frac' not in varl:
    varl.append('NPF_frac')

# %% [markdown]
# ### Get difference from PI to PD

# %%
relative = False
dic_diff = get_diff_by_type(case_dic, varl, ctrl='PI', case_types=['PI', 'PIaerPD'],
                            # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                            relative=relative)

dic_diff.keys()
di_dic = dic_diff['PIaerPD-PI']
di_dic.keys()

# %% [markdown]
# ### Organize data in easy to use format:

# %%
relative = False

dic_abs = get_abs_by_type(case_dic,
                          case_types=['PI', 'PD'],
                          mod_types=None)

# %%
dic_abs.keys()

# %%
dic_abs['PI'].keys()

# %% [markdown]
# ## Plots:

# %%
norm_dic = dict(
    SOA_LV=colors.SymLogNorm(vmin=-5e-1, vmax=5e-1, linthresh=.01, base=10, linscale=.4),
    H2SO4=colors.SymLogNorm(vmin=-5e-1, vmax=5e-1, linthresh=.01, base=10, linscale=.4),
    NCONC01=colors.SymLogNorm(vmin=-1e3, vmax=1e3, linthresh=10, base=10, linscale=.6),
    NMR01=colors.SymLogNorm(vmin=-10, vmax=10, linthresh=5, base=10),  # linscale=.5),
    AWNC_incld=colors.SymLogNorm(vmin=-50, vmax=50, linthresh=1, base=10),
    ACTNL_incld=colors.SymLogNorm(vmin=-40, vmax=40, linthresh=1, linscale=0.5, base=10),
    AREL_incld=colors.SymLogNorm(vmin=-5, vmax=5, linthresh=.1, base=10),
    ACTREL_incld=colors.SymLogNorm(vmin=-7, vmax=7, linthresh=.1, base=10, linscale=0.6),
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
    ACTNL_incld=colors.Normalize(vmin=-20, vmax=20),
    N50=colors.Normalize(vmin=-45, vmax=45),
    N100=colors.Normalize(vmin=-20, vmax=20),
    N150=colors.Normalize(vmin=-10, vmax=10),  # colors.Normalize(vmin=-5, vmax=5),
    N200=colors.Normalize(vmin=-15, vmax=15),  # colors.Normalize(vmin=-5, vmax=5),
)

# %%
norm_diff_dic = dict(
    ACTNL_incld=colors.Normalize(vmin=-12, vmax=12),
    N50=colors.Normalize(vmin=-45, vmax=45),
    N100=colors.Normalize(vmin=-10, vmax=10),
    N150=colors.Normalize(vmin=-4, vmax=4),  # colors.Normalize(vmin=-5, vmax=5),
    N200=colors.Normalize(vmin=-5, vmax=5),  # colors.Normalize(vmin=-5, vmax=5),
)

# %%
norm_dic['NCFT_Ghan'] = colors.Normalize(vmin=-1.8, vmax=1.8)
norm_dic['ACNTL_incld'] = colors.SymLogNorm(vmin=-40, vmax=40, linthresh=1, linscale=0.4, base=10)

# %% [markdown]
# ## Set case for ctrl and other

# %%
cases_oth = ['OsloAero$_{imp}$', 'OsloAero$_{def}$']
ctrl = 'OsloAeroSec'

# %%
var = 'NCONC01'


# dic_means_yr = yearly_mean_dic(['NCONC01'], cases, startyear, endyear, pmin, pressure_adjust,
#                                   avg_over_lev=avg_over_lev, groupby=None, dims=('lev',), area='Global')
def add_significance_diff_diff(ax, dic_means_yr, var, case_oth, ctrl='OsloAeroSec', case_types=None):
    if case_types is None:
        case_types = ['PI', 'PD']
    PI = case_types[0]
    PD = case_types[1]
    c1 = get_casen_by_type_mod(PI, ctrl)
    c2 = get_casen_by_type_mod(PD, ctrl)
    da1 = dic_means_yr[c2][var] - dic_means_yr[c1][var]
    co1 = get_casen_by_type_mod(PI, case_oth)
    co2 = get_casen_by_type_mod(PD, case_oth)

    da2 = dic_means_yr[co2][var] - dic_means_yr[co1][var]

    case1 = 'Delta_PD-PI_' + ctrl
    case2 = 'Delta_PD-PI_' + case_oth
    T, data4comp, sig_map, t = calc_significance_field(case1, case2, da1, da2, 'time', 0.95)
    hatch_area_sign(t, T, ax, transform=ccrs.PlateCarree(),
                    hatches=['...', ''],
                    hatch_lw=.7

                    )
    return T, data4comp, sig_map, t


def add_significance_diff(ax, dic_means_yr, var, ctrl='OsloAeroSec', case_types=None):
    if case_types is None:
        case_types = ['PI', 'PD']
    PI = case_types[0]
    PD = case_types[1]
    c1 = get_casen_by_type_mod(PI, ctrl)
    c2 = get_casen_by_type_mod(PD, ctrl)
    da1 = dic_means_yr[c2][var]
    da2 = dic_means_yr[c1][var]  # dic_means_yr[co2][var]-dic_means_yr[co1][var]

    case1 = 'Delta_PD-PI_' + ctrl
    case2 = 'Delta_PD-PI_' + case_oth
    T, data4comp, sig_map, t = calc_significance_field(case1, case2, da1, da2, 'time', .95)
    hatch_area_sign(t, T, ax, transform=ccrs.PlateCarree(),
                    # hatches=hatches,
                    hatches=['...', ''],
                    hatch_lw=.7
                    # hatch_lw = hatch_lw,
                    # reverse=reverse
                    )
    return T, data4comp, sig_map, t


# %%
vl = ['NCONC01', 'ACTNL_incld', 'NCFT_Ghan']
fig, axs_dic = abs_diffs(di_dic,
                         vl,
                         norm_abs=norm_abs,
                         norm_dic=norm_dic,
                         invert_diff=True,

                         )

dic_means_yr = yearly_mean_dic(vl + ['ACTNL'], cases, startyear, endyear, pmin, pressure_adjust,
                               avg_over_lev=avg_over_lev, groupby=None, dims=('lev',), area='Global')

for var in vl:
    for case_oth in cases_oth:
        ax = axs_dic[var][case_oth]

        add_significance_diff_diff(ax, dic_means_yr, var, case_oth, ctrl='OsloAeroSec', case_types=None)
    ax = axs_dic[var][ctrl]
    add_significance_diff(ax, dic_means_yr, var, ctrl='OsloAeroSec', case_types=None)
clear_output()
fn = filen_base + 'aer_cld_frc.pdf'
fig.savefig(fn, dpi=300)
plt.show()

# %%
print(fn)

# %%

# %%
vl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan']  # ['NCONC01', 'ACTNL_incld', 'NCFT_Ghan']
fig, axs_dic = abs_diffs(di_dic,
                         vl,
                         norm_abs=norm_abs,
                         norm_dic=norm_dic,
                         invert_diff=True,

                         )

dic_means_yr = yearly_mean_dic(vl, cases, startyear, endyear, pmin, pressure_adjust,
                               avg_over_lev=avg_over_lev, groupby=None, dims=('lev',), area='Global')

for var in vl:
    for case_oth in cases_oth:
        ax = axs_dic[var][case_oth]

        add_significance_diff_diff(ax, dic_means_yr, var, case_oth, ctrl='OsloAeroSec', case_types=None)
    ax = axs_dic[var][ctrl]
    add_significance_diff(ax, dic_means_yr, var, ctrl='OsloAeroSec', case_types=None)
clear_output()
fp = filen_base + 'ERF_aci_SWLW.pdf'
fig.savefig(fp, dpi=300)

plt.show()

# %%
print(fp)
