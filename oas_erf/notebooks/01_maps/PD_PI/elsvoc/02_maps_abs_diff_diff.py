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
# # Plots abs, diff for PI and PD separate

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

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass

from oas_erf.util.plot.maps_PIPD import abs_diffs_PI_PD_sep

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
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_svoc_smax',
    'NF1850_aeroxid2014_SECT_svoc_smax',
    'NF1850_SECT_elvoc_smax',
    'NF1850_aeroxid2014_SECT_elvoc_smax',

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
version = 'pi_pd_diff-elsvoc'
plot_path = get_plotpath('maps')
filen_base = plot_path + '/%s' % version
# print(plot_path)
make_folders(plot_path)

# %%
print(filen_base)

# %% [markdown]
# ### Variables to load:

# %%
varl = ['NCONC01', 'NMR01', 'N_AER', 'ACTNL_incld', 'ACTREL_incld', 'CDNUMC', 'cb_NA',
        'cb_SOA_NA', 'cb_SO4_NA', 'AWNC_incld', 'AREL_incld',
        'SWCF_Ghan', 'LWCF_Ghan', 'NCFT_Ghan']
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

# %%
relative = False
dic_diff = get_diff_by_type(case_dic, varl, ctrl='PI', case_types=['PI', 'PIaerPD'],
                            # mod_types=['OsloAeroSec','OsloAero$_{def}$'],
                            relative=relative)

dic_diff.keys()
di_dic = dic_diff['PIaerPD-PI']
di_dic.keys()

# %%

# %% [markdown]
# ### Organize data in easy to use format:

# %%
relative = False

dic_abs = get_abs_by_type(case_dic,
                          case_types=['PI', 'PIaerPD'],
                          mod_types=['OsloAeroSec', 'OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$', 'OsloAero$_{imp}$',
                                     'OsloAero$_{def}$'], )

# %% [markdown]
# ## Plots:

# %%
var = 'NPF_frac'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$', 'OsloAero$_{imp}$',
                                              'OsloAero$_{def}$'],
                                   # , 'ACTNL_incld', 'ACTREL_incld'],
                                   # norm_abs=norm_abs,
                                   # norm_dic=norm_dic
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'}
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

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
norm_dic = dict()
norm_dic['NCFT_Ghan'] = colors.Normalize(vmin=-1.8, vmax=1.8)
norm_dic['ACNTL_incld'] = colors.SymLogNorm(vmin=-40, vmax=40, linthresh=1, linscale=0.4, base=10)

# %%
var = 'NCFT_Ghan'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=False,
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$'],
                                   # ,'OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                                   norm_diff=norm_dic['NCFT_Ghan'],
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
# axs =[axs_dic['PI'][c] for c in ['OsloAeroSec','OsloAero$_{imp}$', 'OsloAero$_{def}$']]#[cases_oth]
# subp_insert_abc(np.array(axs), pos_x=0.01,pos_y=1.0)
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'SWCF_Ghan'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=False,
                                   norm_diff=norm_dic['NCFT_Ghan'],
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$'],
                                   # ,'OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'LWCF_Ghan'
relative = False
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   norm_diff=norm_dic['NCFT_Ghan'],
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$'],
                                   # ,'OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                                   )
fn = filen_base + f'{var}_PIPD_sep_rel{relative}.pdf'
fig.savefig(fn, dpi=300)
plt.show()

# %%
print(fn)
# %%
var = 'NCONC01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=False,
                                   norm_diff=colors.Normalize(vmin=-100, vmax=100),
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$'],
                                   # 'OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'ACTNL_incld'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=False,
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$'],
                                   # 'OsloAero$_{imp}$', 'OsloAero$_{def}$'],

                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'CDNUMC'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=False,
                                   cases_oth=['OsloAeroSec$_{svoc}$', 'OsloAeroSec$_{elvoc}$', 'OsloAero$_{imp}$',
                                              'OsloAero$_{def}$'],

                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
# fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N50'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   # norm_diff=colors.Normalize(vmin=-30, vmax=30)
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N100'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N150'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N200'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   norm_diff=colors.Normalize(vmin=-13, vmax=+13)
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N250'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   norm_diff=colors.Normalize(vmin=-13, vmax=+13)
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'AWNC_incld'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True,
                                   # norm_diff=colors.Normalize(vmin=-100, vmax=100)
                                   )
fn = filen_base + f'{var}_PIPD_sep_rel{relative}.pdf'
fig.savefig(fn, dpi=300)
plt.show()

# %%

print(fn)
# %%
var = 'SIGMA01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,

                                   # , 'ACTNL_incld', 'ACTREL_incld'],
                                   # norm_abs=norm_abs,
                                   # norm_dic=norm_dic
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'}
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NMR01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,

                                   # , 'ACTNL_incld', 'ACTREL_incld'],
                                   # norm_abs=norm_abs,
                                   # norm_dic=norm_dic
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'}
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NMR01'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NCONC08'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True

                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NMR08'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NMR09'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'NCONC09'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%
var = 'N_AER'
relative = True
fig, axs_dic = abs_diffs_PI_PD_sep(dic_abs,
                                   var,
                                   relative=relative,
                                   type_nndic={'PI': 'Pre-industrial', 'PIaerPD': 'Present day'},
                                   switch_diff=True
                                   )
fig.savefig(filen_base + f'{var}_PIPD_sep_rel{relative}.pdf', dpi=300)
plt.show()

# %%

# %%
