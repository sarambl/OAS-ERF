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
# from IPython import get_ipython

# %%
from oas_erf.constants import get_plotpath
from oas_erf.util.naming_conventions.var_info import get_fancy_var_name
from oas_erf.util.practical_functions import make_folders
from oas_erf.util.slice_average import one_val_tab
# noinspection PyBroadException
from oas_erf.util.slice_average.one_val_tab import get_mean_std_by_type
from IPython import get_ipython

try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import numpy as np

from oas_erf.data_info import simulation_types

# %% [markdown]
# ## Filenames

# %%
from oas_erf.util.plot.colors import get_case_col

version = 'clean'
plt_path = get_plotpath('one_value')


def create_filename(name):
    fn = f'{plt_path}_bar_{version}_{name}.'
    make_folders(fn)
    return fn


# %%
startyear = '0004-01'
endyear = '0008-12'
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev

# %% [markdown]
# ## Models and variables to include:

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
# %%
model_types = ['OsloAeroSec', 'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]

cdic = {key: get_case_col(key) for key in model_types}

# %%
trans_dic = {v: get_fancy_var_name(v) for v in varl}

# %%
rn_dic = {
    'NCRE$_{Ghan}$': 'ERF$_{aci}$',
    'SWCRE$_{Ghan}$': 'ERF$_{aci,SW}$',
    'LWCRE$_{Ghan}$': 'ERF$_{aci,LW}$',
    'DRE$_{Ghan}$': 'ERF$_{ari}$'

}

# %% [markdown]
# ## Import data:

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

# %% [markdown]
# ## Get PI to PD difference

# %%
relative = False
dic_diff = simulation_types.get_diff_by_type(dic_vals, varl, case_types=['PI', 'PD'],
                                             relative=relative,
                                             mod_types=model_types,
                                             ctrl='PI'
                                             )['PD-PI']

ls = []
for key in dic_diff.keys():
    print(key)
    _df = dic_diff[key]

    _df['case'] = key
    print(_df.keys())
    ls.append(_df.reset_index())

df_tot = pd.concat(ls)
df_tot.head()

# %%

svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan', 'DIR_Ghan']

df1 = df_tot[[*svarl, 'case']]
df1 = df1.rename(trans_dic, axis=1)
df1 = df1.rename(rn_dic, axis=1)
df2 = pd.melt(df1, id_vars='case')

df2.groupby(['variable', 'case']).mean()

# %%
cols = [cdic[c] for c in df_tot['case'].unique()]
print(cols)

# %%

figsize = [4., 4.7]
f, ax = plt.subplots(figsize=figsize, dpi=150)

pts = np.linspace(0, np.pi * 2, 24)
circ = np.c_[np.sin(pts) / 2, -np.cos(pts) / 2]
vert = np.r_[circ, circ[::-1] * .7]
open_circle = mpl.path.Path(vert)
sns.barplot(y='value', x='variable', hue='case', data=df2, errcolor='.55', ci=90,
            palette=cols, alpha=.6)
g = sns.stripplot(y='value', x='variable', hue='case', data=df2,
                  dodge=True, marker=open_circle, palette=['none'] * 3,
                  jitter=.2, label='_nolegend_')  # , add_=False)#label='_nolabel_')

handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
l = plt.legend(handles[3:], labels[3:], frameon=False, loc='center right')

ax.axhline(0, linewidth=0.4, c='k')
sns.despine(f, bottom=True, trim=True, offset=10)
# ax.legend(frameon=False)
plt.tick_params(
    axis='x',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    bottom=False,  # ticks along the bottom edge are off
    top=False,  # ticks along the top edge are off
    labeltop=True,
    labelbottom=False
)  # labels along the bottom edge are off
ax.set_ylabel('[Wm$^{-2}$]')
ax.set_xlabel('')

fn = create_filename('forcing')
plt.tight_layout()
f.savefig(fn + 'pdf', dpi=300)
print(fn)
plt.show()

print(fn)

# %% [markdown]
# ## Compute relative difference

# %%

relative = True
dic_diff = simulation_types.get_diff_by_type(dic_vals, varl, case_types=['PI', 'PD'],
                                             relative=relative,
                                             mod_types=model_types,
                                             ctrl='PI'
                                             )['(PD-PI)/PI']

ls = []
for key in dic_diff.keys():
    print(key)
    _df = dic_diff[key]

    _df['case'] = key
    print(_df.keys())
    ls.append(_df.reset_index())

df_tot = pd.concat(ls)
df_tot.head()

# %% [markdown]
# ### Rename and extract

# %%

# svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']
svarl = ['N_AER', 'NCONC01', 'cb_SOA_NA', 'cb_SO4_NA', 'cb_NA']
rn_dic = {
    'c.b. SOA$_{NPF}$+SO4$_{NPF}$': 'c.b. SOA$_{NPF}$ \n+ SO4$_{NPF}$',
    'c.b. SOA$_{NPF}$': 'c.b. \n SOA$_{NPF}$',
    'c.b. SO4$_{NPF}$': 'c.b. \n SO4$_{NPF}$',
}
df1 = df_tot[[*svarl, 'case']]
df1 = df1.rename(trans_dic, axis=1)
df1 = df1.rename(rn_dic, axis=1)
df2 = pd.melt(df1, id_vars='case')

# %%
figsize = [5, 3]

figsize = [4.2, 4.95]
f, ax = plt.subplots(figsize=figsize, dpi=150)

pts = np.linspace(0, np.pi * 2, 24)
circ = np.c_[np.sin(pts) / 2, -np.cos(pts) / 2]
vert = np.r_[circ, circ[::-1] * .7]
open_circle = mpl.path.Path(vert)
sns.barplot(y='value', x='variable', hue='case', data=df2, errcolor='.55',
            ci=95, palette=cols, alpha=.6)
g = sns.stripplot(y='value', x='variable', hue='case', data=df2, dodge=True, marker=open_circle, palette=['none'] * 3,
                  jitter=.2, label='_nolegend_')  # , add_=False)#label='_nolabel_')

handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
l = plt.legend(handles[3:], labels[3:], frameon=False, loc='upper left')

ax.axhline(0, linewidth=0.4, c='k')
sns.despine(f, bottom=True, trim=True, offset=10)
plt.tick_params(
    axis='x',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    bottom=False,  # ticks along the bottom edge are off
    top=False,  # ticks along the top edge are off
    labeltop=True,
    labelbottom=False
)  # labels along the bottom edge are off

ax.set_ylabel('(PDaer-PI)/PI [%]')
fn = create_filename('aerosol')
plt.tight_layout()
ax.set_xlabel('')
f.savefig(fn + 'pdf', dpi=300)
plt.show()

# %%


# svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']
svarl = ['N_AER', 'NCONC01', 'cb_SOA_NA', 'cb_SO4_NA', 'cb_NA']
svarl = ['CDNUMC', 'ACTREL_incld', 'TGCLDCWP']  # 'cb_SOA_NA', 'cb_SO4_NA']

rn_dic = {
    'c.b. SOA$_{NPF}$+SO4$_{NPF}$': 'c.b. SOA$_{NPF}$ \n+ SO4$_{NPF}$',
    'c.b. SOA$_{NPF}$': 'c.b. \n SOA$_{NPF}$',
    'c.b. SO4$_{NPF}$': 'c.b. \n SO4$_{NPF}$',
}
df1 = df_tot[[*svarl, 'case']]
df1 = df1.rename(trans_dic, axis=1)
df1 = df1.rename(rn_dic, axis=1)
df2 = pd.melt(df1, id_vars='case')

# %%

svarl = ['CDNUMC', 'ACTREL_incld', 'TGCLDCWP']  # 'cb_SOA_NA', 'cb_SO4_NA']
relative = True
figsize = [4., 4.8]
f, ax = plt.subplots(figsize=figsize, dpi=150)

pts = np.linspace(0, np.pi * 2, 24)
circ = np.c_[np.sin(pts) / 2, -np.cos(pts) / 2]
vert = np.r_[circ, circ[::-1] * .7]
open_circle = mpl.path.Path(vert)
sns.barplot(y='value', x='variable', hue='case', data=df2, errcolor='.55', ci=95,
            palette=cols, alpha=.6)
g = sns.stripplot(y='value', x='variable', hue='case', data=df2, dodge=True, marker=open_circle, palette=['none'] * 3,
                  jitter=.2, label='_nolegend_')  # , add_=False)#label='_nolabel_')

handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
l = plt.legend(handles[3:], labels[3:], frameon=False, loc='upper right')

ax.axhline(0, linewidth=0.4, c='k')
sns.despine(f, bottom=True, trim=True, offset=10)
plt.tick_params(
    axis='x',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    bottom=False,  # ticks along the bottom edge are off
    top=False,  # ticks along the top edge are off
    labeltop=True,
    labelbottom=False
)  # labels along the bottom edge are off

ax.set_ylabel('(PDaer-PI)/PI [%]')
fn = create_filename('cld_props')
plt.tight_layout()
f.savefig(fn + 'pdf', dpi=300)
plt.show()

# %%


case_types = ['PI', 'PIaerPD']
model_types = ['OsloAeroSec', 'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]

cdic = {key: get_case_col(key) for key in model_types}  # , ['r','g','b'])}
# %%
figsize = [5, 3]
svarl = ['FSNT', 'FSNT_DRF', 'FLNT', 'FLNT_DRF', 'FSNTCDRF']  # '', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']

mean, std, mean_nn, std_nn = get_mean_std_by_type(dic_vals, svarl, case_types=['PI', 'PIaerPD'],
                                                  model_types=model_types, ctrl='PI')

fig, ax = plt.subplots(1, figsize=figsize)
mean_nn.plot.bar(alpha=0.5, ax=ax, color=[cdic[c] for c in pd.DataFrame(mean).columns],
                 yerr=std_nn)  # , colors={'OsloAeroSec':'b'})
ax.axhline(0, linewidth=0.4, c='k')
sns.despine(fig, bottom=True, trim=True, offset=10)
ax.legend(frameon=False)
plt.tick_params(
    axis='x',  # changes apply to the x-axis
    which='both',  # both major and minor ticks are affected
    bottom=False,  # ticks along the bottom edge are off
    top=False,  # ticks along the top edge are off
    labeltop=True,
    labelbottom=False
)  # labels along the bottom edge are off
ax.set_ylabel('PDaer-PI [Wm$^{-2}$]')
fn = create_filename('forcing')
plt.tight_layout()
# fig.savefig(fn+'pdf', dpi=300)
plt.show()
# %%
