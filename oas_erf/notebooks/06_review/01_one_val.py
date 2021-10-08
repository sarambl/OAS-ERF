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

# %%
from oas_erf.constants import get_plotpath
from oas_erf.util.practical_functions import make_folders
from oas_erf.util.slice_average import one_val_tab
# load and autoreload
from IPython import get_ipython
from oas_erf.util.naming_conventions.var_info import get_fancy_var_name, get_fancy_unit_xr

# noinspection PyBroadException
from oas_erf.util.slice_average.one_val_tab import plt_var, get_diff_by_type, get_mean_std_by_type

try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass
# %%
from oas_erf.util.plot.colors import get_case_col
version='review'
plt_path = get_plotpath('one_value')

def create_filename(name):
    fn= f'{plt_path}_2d_{version}{name}.'
    make_folders(fn)
    return fn

# %%
startyear = '0004-01'
endyear = '0008-12'
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True#True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev


# %%
varl = ['N_AER', 'NCONC01', 'TGCLDCWP', 'CDNUMC', 'NCFT_Ghan', 'DIR_Ghan', 'LWDIR_Ghan', 'SWDIR_Ghan', 'SWCF_Ghan',
        'LWCF_Ghan' ,'cb_SO4_NA','cb_SOA_NA','cb_NA', 'SOA_NA','SO4_NA',
        'ACTNL_incld','ACTREL_incld', 'SFisoprene','SFmonoterp']
cases = ['SECTv21_ctrl_koagD', 'SECTv21_incY', 'SECTv21_decY', 'noSECTv21_ox_ricc_dd', 'noSECTv21_ox_ricc_decY',
         'noSECTv21_ox_ricc_incY']
cases_sec = [
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_ctrl',
    #'NF1850_aeroxid2014_SECT_ctrl_smax',
    'NF1850_SECT_paas',
    'NF1850_SECT_depT',
    'NF1850_aeroxid2014_SECT_paas',
    'NF1850_aeroxid2014_SECT_depT',
]
cases_nsec = [
    #'NF1850_noSECT_def_smax',
    #'NF1850_aeroxid2014_noSECT_def_smax',
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_noSECT_ox_ricc',
    'NF1850_aeroxid2014_noSECT_ox_ricc',
    'NF1850_noSECT_ox_ricc_depT',
    'NF1850_aeroxid2014_noSECT_ox_ricc_depT',
]

# %%
varl_ex = ['FSNT','FSNT_DRF','FLNT','FLNT_DRF','FSNTCDRF']
varl =varl+ varl_ex

# %%
case_types = ['PI', 'PD']
model_types = ['OsloAeroSec','OsloAeroSec$_{paas}$','OsloAeroSec$_{depT}$', 'OsloAero$_{depT}$',   'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]
mod_types = ['OsloAeroSec','OsloAeroSec$_{paas}$','OsloAeroSec$_{depT}$', 'OsloAero$_{depT}$',   'OsloAero$_{imp}$', 'OsloAero$_{def}$'][::-1]
cdic = {key: get_case_col(key) for key in mod_types}  # , ['r','g','b'])}

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
from useful_scit.imps import ( sns,pd, plt)

# %%
from oas_erf.data_info.simulation_types import get_abs_by_type

# %%
di = get_abs_by_type(dic_vals, case_types=['PI','PD'], mod_types=model_types)
di.keys()

# %%
ls =[]
for ct in  di.keys():
    _di = di[ct]
    for cn in _di.keys():
        print(cn)
        _df = _di[cn]
        _df['case']=cn
        _df['case_type']=ct
    
        ls.append(_df.reset_index())
    

df_tot = pd.concat(ls)
df_tot

# %%
svarl = ['ACTNL_incld','cb_NA']#,'FSNT_DRF','FLNT','FLNT_DRF','FSNTCDRF']#'', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']
v1 = 'ACTNL_incld'
v2 = 'ACTREL_incld'
di = get_abs_by_type(dic_vals, case_types=['PI','PD'], mod_types=model_types)
_df = pd.DataFrame(columns=['val', 'type','var','model'])
di_var = {}
for t in di.keys():
    for m in di[t].keys():
        di[t][m]=di[t][m].mean()

for v in varl:
    _df_v = pd.DataFrame(columns=['val', 'type','model'])
    for t in di.keys():
        for m in di[t].keys():

            _df = _df.append(pd.DataFrame([di[t][m][v], t, v, m],index=['val', 'type','var','model']).transpose(), ignore_index=True)
            _df_v=_df_v.append(pd.DataFrame([di[t][m][v], t, m],index=['val', 'type','model']).transpose(), ignore_index=True)
            _df_v['val'] = pd.to_numeric(_df_v['val'])
    di_var[v]=_df_v.copy()


# %%
di_2 = di_var[v].set_index(['type','model']).rename({'val':v}, axis=1)#.plot.bar()
for v in varl:
    di_2[v] = di_var[v].set_index(['type','model']).rename({'val':v}, axis=1)

# %%
di_2.reset_index()

# %%

# %%
df_tot.head()

# %%
mod_types

# %%
y='NCFT_Ghan'
x='cb_NA'
_df = di_2.reset_index()
fig, ax = plt.subplots(1, figsize=[6,3], dpi=150)
for mod in model_types:
    print(mod)
    sdf = _df[_df['model']==mod]
    c=get_case_col(mod)
    plt.plot(sdf[x],sdf[y], label=mod, marker='o', c=c)
    
    plt.text(sdf[x].iloc[0],sdf[y].iloc[0]+0.05, sdf['type'].iloc[0], color=c, size=12)
    plt.text(sdf[x].iloc[1],sdf[y].iloc[1]+0.05, sdf['type'].iloc[1], color=c, size=12)

for mod in model_types:
    c=get_case_col(mod)
    
    _df = df_tot[df_tot['case']==mod][[x,y]]
    plt.scatter(_df[x],_df[y], alpha=.4, edgecolor=c,facecolors='none',)

def label(v):
    _n = get_fancy_var_name(v)
    _u = df2[v]['unit']
    return f'{_n} [{_u}]'
plt.xlabel(label(x))#'Col. burden N$_{NPF}$ [kg m$^{-2}$]')
plt.ylabel(label(y))

plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left',)
sns.despine(fig)
fn= create_filename(f'{x}_{y}')
plt.tight_layout()
plt.savefig(fn+'pdf', dpi=300)


# %%
_df = di_2.reset_index()
_df[_df['model']=='OsloAero$_{depT}$']#"=='OsloAeroSec']

# %%
df_tot[df_tot['case']=='OsloAeroSec$_{depT}$'][['SWCF_Ghan','LWCF_Ghan','NCFT_Ghan']]#[df_tot['case']=='OsloAero$_{depT}$']

# %%
df_tot[df_tot['case']=='OsloAero$_{depT}$'][['SWCF_Ghan','LWCF_Ghan','NCFT_Ghan']]#[df_tot['case']=='OsloAero$_{depT}$']

# %%
_df['model']

# %%
trans_dic = {v:get_fancy_var_name(v) for v in varl}

# %%
rn_dic = {
    'NCRE$_{Ghan}$':'ERF$_{aci}$',
    'SWCRE$_{Ghan}$':'ERF$_{aci,SW}$',
    'LWCRE$_{Ghan}$':'ERF$_{aci,LW}$',
    'DRE$_{Ghan}$': 'ERF$_{ari}$'
    
    
}

# %%

svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']

df1 = df_tot[[*svarl,'case']]
df1 = df1.rename(trans_dic, axis=1)
df1 = df1.rename(rn_dic, axis=1)
df2 = pd.melt(df1,id_vars='case')


# %%
cols=[cdic[c] for c in df_tot['case'].unique()]
cols


# %%

# %%
figsize=[5,3]
svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']

mean, std, mean_nn, std_nn = get_mean_std_by_type(dic_vals, svarl, case_types=['PI', 'PIaerPD'], ctrl='PI',model_types=model_types)
mean_nn = mean_nn.rename(rn_dic, axis=0)
std_nn = std_nn.rename(rn_dic, axis=0)
fig, ax = plt.subplots(1, figsize=figsize, dpi=200)
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
ax.set_ylabel('PD-PI [Wm$^{-2}$]')
fn = create_filename('forcing')
plt.tight_layout()
fig.savefig(fn+'pdf', dpi=300)
plt.show()
# %%
fn


# %%

# %%

# %%

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
from useful_scit.imps import (plt, pd, sns)

# %%
from oas_erf.data_info.simulation_types import get_abs_by_type

# %%
from oas_erf.data_info import get_nice_name_case, simulation_types


# %%
relative=False
dic_diff = simulation_types.get_diff_by_type(dic_vals, varl, case_types=['PI','PD'],
                 relative=relative, 
                 mod_types=model_types,
                 ctrl='PI'
                )['PD-PI']

ls =[]
for key in dic_diff.keys():
    print(key)
    _df = dic_diff[key]
    
    _df['case'] = key
    print(_df.keys())
    ls.append(_df.reset_index())
    

df_tot = pd.concat(ls)
df_tot

# %%
trans_dic = {v:get_fancy_var_name(v) for v in varl}

# %%
rn_dic = {
    'NCRE$_{Ghan}$':'ERF$_{aci}$',
    'SWCRE$_{Ghan}$':'ERF$_{aci,SW}$',
    'LWCRE$_{Ghan}$':'ERF$_{aci,LW}$',
    'DRE$_{Ghan}$': 'ERF$_{ari}$'
    
    
}

# %%

svarl = ['NCFT_Ghan', 'SWCF_Ghan', 'LWCF_Ghan','DIR_Ghan']

df1 = df_tot[[*svarl,'case']]
df1 = df1.rename(trans_dic, axis=1)
df1 = df1.rename(rn_dic, axis=1)
df2 = pd.melt(df1,id_vars='case')


# %%
cols=[cdic[c] for c in df_tot['case'].unique()]
cols


# %%
df2.groupby(['variable','case']).mean()

# %%
import matplotlib as mpl
import numpy as np
figsize=[5,3]
figsize=[4.,4.7]
f,ax = plt.subplots(figsize=figsize,dpi=150)

pts = np.linspace(0, np.pi * 2, 24)
circ = np.c_[np.sin(pts) / 2, -np.cos(pts) / 2]
vert = np.r_[circ, circ[::-1] * .7]
open_circle = mpl.path.Path(vert)
sns.barplot(y='value',x='variable', hue='case', data=df2,errcolor='.55',ci=90
           ,palette=cols, alpha=.6)
g = sns.stripplot(y='value',x='variable', hue='case', data=df2,dodge=True,marker=open_circle,palette=['none']*3,jitter=.2,label='_nolegend_')#, add_=False)#label='_nolabel_')


handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
l = plt.legend(handles[3:], labels[3:],  frameon=False, loc='center right')


ax.axhline(0, linewidth=0.4, c='k')
sns.despine(f, bottom=True, trim=True, offset=10)
#ax.legend(frameon=False)
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
f.savefig(fn+'pdf', dpi=300)
print(fn)
plt.show()

print(fn)
