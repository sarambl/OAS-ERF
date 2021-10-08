from oas_erf.util.Nd.sizedist_class_v2.SizedistributionBins import SizedistributionBins


# load and autoreload
from IPython import get_ipython

# noinspection PyBroadException
try:
    _ipython = get_ipython()
    _magic = _ipython.magic
    _magic('load_ext autoreload')
    _magic('autoreload 2')
except:
    pass
# %%
## Edit and specify period (0004-01--0008-12 or 0004-01-0005-12)
startyear = '0004-01'
endyear = '0005-12'


# %%
## Edit cases to compute
cases_sec = [
    'NF1850_SECT_ctrl',
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_ctrl_smax',
    'NF1850_aeroxid2014_SECT_ctrl_smax',

    'NF1850_SECT_elvoc_smax',
    'NF1850_aeroxid2014_SECT_elvoc_smax',

    'NF1850_SECT_paas',
    'NF1850_aeroxid2014_SECT_paas',
    'NF1850_SECT_svoc_smax',
    'NF1850_aeroxid2014_SECT_svoc_smax',
    'NF1850_SECT_depT',
    'NF1850_aeroxid2014_SECT_depT',
]
cases_orig = [
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_noSECT_ox_ricc',
    'NF1850_aeroxid2014_noSECT_ox_ricc',

    'NF1850_noSECT_def_smax',
    'NF1850_aeroxid2014_noSECT_def_smax',

    'NF1850_noSECT_ox_ricc_smax',
    'NF1850_aeroxid2014_noSECT_ox_ricc_smax',


    'NF1850_noSECT_ox_ricc_depT',
    'NF1850_aeroxid2014_noSECT_ox_ricc_depT'

]
cases = cases_sec + cases_orig
# %%

model = 'NorESM'
p_level = 1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True # True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013., 900., 800., 700., 600.]  # used if not avg

# %%

SizedistributionBins
for case_name in cases_sec:
    isSec = True
    sdb = SizedistributionBins(case_name, startyear, endyear, [5., 39.6], isSec, 'month',
                               # raw_data_path=constants.get_input_datapath(),
                               space_res='full',
                               nr_bins=5,
                               print_stat=True,
                               # model_name='NorESM',
                               history_field='.h0.',
                               # locations=constants.locations,
                               # chunks={'diameter': 20},
                               # use_pressure_coords=True,
                               use_eusaar_diam=True
                               )

    sdb.compute_Nd_vars()
for case_name in cases_orig:
    isSec = False
    sdb = SizedistributionBins(case_name, startyear, endyear, [5., 39.6], isSec, 'month',
                               # raw_data_path=constants.get_input_datapath(),
                               space_res='full',
                               nr_bins=5,
                               print_stat=True,
                               # model_name='NorESM',
                               history_field='.h0.',
                               # locations=constants.locations,
                               # chunks={'diameter': 20},
                               # use_pressure_coords=True,
                               use_eusaar_diam=True
                               )

    sdb.compute_Nd_vars()

    # %%
