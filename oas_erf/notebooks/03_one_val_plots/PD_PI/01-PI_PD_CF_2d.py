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
# # Bar plots global means

# %% [markdown]
# ## Div imports
#
# load and autoreload

# %%
from IPython import get_ipython

from oas_erf.constants import get_plotpath
from oas_erf.util.naming_conventions.var_info import get_fancy_var_name
from oas_erf.util.practical_functions import make_folders
from oas_erf.util.slice_average import one_val_tab
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# %%
from oas_erf.data_info.simulation_types import get_abs_by_type

# noinspection PyBroadException

try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass
# %% [markdown]
# ## Filenames

# %%
from oas_erf.util.plot.colors import get_case_col

version = ''
plt_path = get_plotpath('one_value')


def create_filename(name):
    fn = f'{plt_path}_2d_{version}{name}.'
    make_folders(fn)
    return fn


# %% [markdown]
# ## Div settings 

# %%
startyear = '0004-01'
endyear = '0008-12'
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev

# %%
varl = ['N_AER', 'NCONC01', 'TGCLDCWP', 'CDNUMC', 'NCFT_Ghan', 'DIR_Ghan', 'LWDIR_Ghan', 'SWDIR_Ghan', 'SWCF_Ghan',
        'LWCF_Ghan', 'cb_SO4_NA', 'cb_SOA_NA', 'cb_NA', 'SOA_NA', 'SO4_NA',
        'ACTNL_incld', 'ACTREL_incld', 'SFisoprene', 'SFmonoterp']
varl_ex = ['FSNT', 'FSNT_DRF', 'FLNT', 'FLNT_DRF', 'FSNTCDRF']
varl = varl + varl_ex
cases = ['SECTv21_ctrl_koagD', 'SECTv21_incY', 'SECTv21_decY', 'noSECTv21_ox_ricc_dd', 'noSECTv21_ox_ricc_decY',
         'noSECTv21_ox_ricc_incY']
cases_sec = [
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl',
]
cases_nsec = [
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_aeroxid2014_noSECT_ox_ricc',
    'NF1850_noSECT_ox_ricc',

]
# 'noSECTv21_ox_ricc_decY','noSECTv21_ox_ricc_dd','noSECTv21_ox_ricc_incY',
#               'noSECTv21_def_decY','noSECTv21_default_dd','noSECTv21_def_incY']
# %% [markdown]
# ## Import yearly means

# %%
# varl = ['N_AER']
df2, dic_vals = one_val_tab.get_tab_yearly_mean(varl,
                                                cases_sec + cases_nsec,
                                                startyear,
                                                endyear,
                                                pmin=pmin,
                                                pressure_adjust=pressure_adjust,
                                                average_over_lev=avg_over_lev,
                                                groupby='time.year',  # 'time',
                                                dims=None,
                                                area='Global',
                                                invert_dic=True
                                                )

# %%

# %%
di = get_abs_by_type(dic_vals, case_types=['PI', 'PD'])
di.keys()

# %% [markdown]
# ## Concatinate data in one dataframe

# %%
ls = []
for ct in di.keys():
    _di = di[ct]
    for cn in _di.keys():
        print(cn)
        _df = _di[cn]
        _df['case'] = cn
        _df['case_type'] = ct

        ls.append(_df.reset_index())

df_tot = pd.concat(ls)
df_tot.head()

# %% [markdown]
# ## Reorganize data

# %%
di = get_abs_by_type(dic_vals, case_types=['PI', 'PD'])
_df = pd.DataFrame(columns=['val', 'type', 'var', 'model'])
di_var = {}
for t in di.keys():
    for m in di[t].keys():
        di[t][m] = di[t][m].mean()

for v in varl:
    _df_v = pd.DataFrame(columns=['val', 'type', 'model'])
    for t in di.keys():
        for m in di[t].keys():
            _df = _df.append(pd.DataFrame([di[t][m][v], t, v, m], index=['val', 'type', 'var', 'model']).transpose(),
                             ignore_index=True)
            _df_v = _df_v.append(pd.DataFrame([di[t][m][v], t, m], index=['val', 'type', 'model']).transpose(),
                                 ignore_index=True)
            _df_v['val'] = pd.to_numeric(_df_v['val'])
    di_var[v] = _df_v.copy()

di_2 = di_var[v].set_index(['type', 'model']).rename({'val': v}, axis=1)  # .plot.bar()
for v in varl:
    di_2[v] = di_var[v].set_index(['type', 'model']).rename({'val': v}, axis=1)

# %%
di_2.reset_index()

# %%


case_types = ['PI', 'PIaerPD']
model_types = ['OsloAeroSec', 'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]

cdic = {key: get_case_col(key) for key in model_types}  # , ['r','g','b'])}
# %% [markdown]
# ## Plots

# %%
y = 'NCFT_Ghan'
x = 'cb_NA'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[4, 3], dpi=150)
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0] + 0.05, sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1] + 0.05, sdf['type'].iloc[1], color=c, size=12)

for mod in model_types:
    c = get_case_col(mod)

    _df = df_tot[df_tot['case'] == mod][[x, y]]
    plt.scatter(_df[x], _df[y], alpha=.4, edgecolor=c, facecolors='none', )


def label(v):
    _n = get_fancy_var_name(v)
    _u = df2[v]['unit']
    return f'{_n} [{_u}]'


plt.xlabel(label(x))  # 'Col. burden N$_{NPF}$ [kg m$^{-2}$]')
plt.ylabel(label(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %% [markdown]
# ### EXTRA:

# %%
y = 'NCFT_Ghan'
x = 'cb_NA'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[4, 3])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0] + 0.05, sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1] + 0.05, sdf['type'].iloc[1], color=c, size=12)


def label(v):
    _n = get_fancy_var_name(v)
    _u = df2[v]['unit']
    return f'{_n} [{_u}]'


plt.xlabel(label(x))
plt.ylabel(label(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
# plt.savefig(fn+'pdf', dpi=300)


# %%
y = 'NCFT_Ghan'
x = 'cb_NA'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%
y = 'NCFT_Ghan'
x = 'NCONC01'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%
y = 'NCFT_Ghan'
x = 'cb_NA'
x = 'cb_NA'
y = 'CDNUMC'

_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%

y = 'NCFT_Ghan'
x = 'N_AER'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%

y = 'ACTNL_incld'
x = 'NCONC01'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%

x = 'ACTNL_incld'
y = 'N_AER'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)

# %%

x = 'ACTNL_incld'
y = 'NCFT_Ghan'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[5, 4])
for mod in model_types:
    sdf = _df[_df['model'] == mod]
    c = get_case_col(mod)
    plt.plot(sdf[x], sdf[y], label=mod, marker='o', c=c)

    plt.text(sdf[x].iloc[0], sdf[y].iloc[0], sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1], sdf[y].iloc[1], sdf['type'].iloc[1], color=c, size=12)
plt.xlabel(get_fancy_var_name(x))
plt.ylabel(get_fancy_var_name(y))

plt.legend(frameon=False)
sns.despine(fig)
fn = create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn + 'pdf', dpi=300)
