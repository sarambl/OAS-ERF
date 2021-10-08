# %%

from useful_scit.imps import *

from oas_erf.util.imports import get_averaged_fields
from oas_erf.util.imports.get_fld_fixed import get_field_fixed

log.ger.setLevel(log.log.DEBUG)

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

# Edit and specify period (0004-01--0008-12 or 0004-01-0005-12)

startyear = '0004-01'
endyear = '0008-12'


# %%
# Edit cases to compute
cases_sec = [
    #'NF1850_aeroxid2014_SECT_ctrl_smax',
    #'NF1850_SECT_ctrl_smax',
    #'NF1850_aeroxid2014_SECT_elvoc_smax',
    #'NF1850_SECT_elvoc_smax',
    #'NF1850_aeroxid2014_SECT_svoc_smax',
    #'NF1850_SECT_svoc_smax',
    'NF1850_aeroxid2014_SECT_ctrl',
    'NF1850_SECT_ctrl',
#    'NF1850_aeroxid2014_SECT_ctrl_smax',
#    'NF1850_SECT_paas',
#    'NF1850_SECT_depT',
#    'NF1850_aeroxid2014_SECT_paas',
#    'NF1850_aeroxid2014_SECT_depT',
]
cases_orig = [
    'NF1850_noSECT_def',
    'NF1850_aeroxid2014_noSECT_def',
    'NF1850_noSECT_ox_ricc',
    'NF1850_aeroxid2014_noSECT_ox_ricc',
#    'NF1850_noSECT_def_smax',
#    'NF1850_aeroxid2014_noSECT_def_smax',
#    'NF1850_noSECT_ox_ricc_depT',
#    'NF1850_aeroxid2014_noSECT_ox_ricc_depT',
#    'NF1850_noSECT_def_smax_svoc',
#    'NF1850_noSECT_ox_ricc_smax_svoc',
#    'NF1850_aeroxid2014_noSECT_def_smax_svoc',
#    'NF1850_aeroxid2014_noSECT_ox_ricc_smax_svoc',
    ]

# %%
cases = cases_sec + cases_orig

# %%
p_level=1013.
pmin = 850.  # minimum pressure level
avg_over_lev = True  # True # True#False#True
pressure_adjust = True  # Can only be false if avg_over_lev false. Plots particular hybrid sigma lev
p_levels = [1013.,900., 800., 700., 600.]  # used if not avg
model = 'NorESM'

# %%
varl = [
    'T','Z3',
    'CCN1',
    'CCN2',
    'CCN3',
    'CCN4',
    'CCN5',
    'CCN6',
    'CCN7',
    'ACTNL_incld',
    'ACTREL_incld',
    'TGCLDCWP',
    'TGCLDIWP',
    'TGCLDLWP',
    'NCFT_Ghan',
    'HYGRO01',
    'SOA_NAcondTend',
    'SO4_NAcondTend',
    'cb_SOA_NA',
    'cb_SO4_NA',
    'HYGRO01',
    'cb_SOA_LV',
    'cb_H2SO4',
    'SO2',
    'DMS',
    'T',
    'isoprene',
    'monoterp',
    'N_AER',
    'NCONC01',
    'NMR01',
    'GR',
    'COAGNUCL',
    'NUCLRATE',
    'FORMRATE',
    'H2SO4',
    'SOA_LV',
    'SOA_SV',
    'SOA_NA',
    'SO4_NA',
    'SOA_A1',
    'NCFT_Ghan',
    'SFisoprene',
    'SFmonoterp',
    'SOA_NA_totLossR',
    'SOA_NA_lifetime',
    'SO4_NA_totLossR',
    'SO4_NA_lifetime',
    'cb_SOA_NA_OCW',
    'cb_SO4_NA_OCW',
    'SO4_NA_OCWDDF',
    'SO4_NA_OCWSFWET',
    'SOA_NA_OCWDDF',
    'SOA_NA_OCWSFWET',
    'cb_SOA_A1',
    'cb_SO4_A1',
    'cb_SOA_NA',
    'cb_SO4_NA',
    'cb_NA',
    'SWCF_Ghan',
    'LWCF_Ghan',
    'AWNC_incld',
    'AREL_incld',
    'CLDHGH',
    'CLDLOW',
    'CLDMED',
    'CLDTOT',
    'CDNUMC',
    'CLDLIQ',
    'CLOUD',
    'CLDICE',
    'DIR_Ghan',
    'CDOD550',
    'SWDIR_Ghan',
    'SO4_NA_mixnuc1',
    'SOA_NA_mixnuc1',
    'SO4_NA_OCW_mixnuc1',
    'SOA_NA_OCW_mixnuc1',
    'SO4_NASFWET',
    'SOA_NASFWET',
    'SO4_NA_OCWSFWET',
    'SOA_NA_OCWSFWET',
    'SOA_NADDF',
    'SO4_NADDF',
    'SO4_NA_OCWDDF',
    'SOA_NA_OCWDDF',
    'SOA_NAcoagTend',
    'SO4_NAcoagTend',
    'NCONC01',
    'NCONC02',
    'NCONC03',
    'NCONC04',
    'NCONC05',
    'NCONC06',
    'NCONC07',
    'NCONC08',
    'NCONC09',
    'NCONC10',
    'NCONC11',
    'NCONC12',
    'NCONC13',
    'NCONC14',
    'N50',
    'N100',
    'N150',
    'N200',
    'N250'
    ]
smax_vars = [
    'Smax',
    'NDROPSRC',
    'NDROPMIX',
    'NDROPCOL',
    'WTKE',
    'Smax_supZero',
    'Smax_w',
    'NACT_FRAC01',
    'NACT_FRAC02',
    'NACT_FRAC03',
    'NACT_FRAC04',
    'NACT_FRAC05',
    'NACT_FRAC06',
    'NACT_FRAC07',
    'NACT_FRAC08',
    'NACT_FRAC09',
    'NACT_FRAC10',
    'NACT_FRAC11',
    'NACT_FRAC12',
    'NACT_FRAC13',
    'NACT_FRAC14',
    'NACT01',
    'NACT02',
    'NACT03',
    'NACT04',
    'NACT05',
    'NACT06',
    'NACT07',       
    'NACT08',
    'NACT09',
    'NACT10',                                
    'NACT11',
    'NACT12',
    'NACT13',
    'NACT14',
]
#varl = varl + smax_vars
varl_sec = [
    'nrSOA_SEC_tot',
    'nrSO4_SEC_tot',
    'nrSEC_tot',
    'cb_SOA_SEC01',
    'cb_SOA_SEC02',
    'cb_SOA_SEC03',
    'leaveSecSOA',
    'leaveSecH2SO4',
  ]
# %%
for case in cases:
    get_field_fixed(case,varl, startyear, endyear, #raw_data_path=constants.get_input_datapath(),
                pressure_adjust=True, model = 'NorESM', history_fld='.h0.', comp='atm', chunks=None)
for case in cases_sec:
    get_field_fixed(case, varl_sec, startyear, endyear, #raw_data_path=constants.get_input_datapath(),
                        pressure_adjust=True, model = 'NorESM', history_fld='.h0.', comp='atm', chunks=None)
    # %%
maps_dic = get_averaged_fields.get_maps_cases(cases,varl,startyear, endyear,
                                              avg_over_lev=avg_over_lev,
                                              pmin=pmin,
                                              pressure_adjust=pressure_adjust, p_level=p_level)
maps_dic = get_averaged_fields.get_maps_cases(cases_sec,varl_sec,startyear, endyear,
                                              avg_over_lev=avg_over_lev,
                                              pmin=pmin,
                                              pressure_adjust=pressure_adjust, p_level=p_level)


for period in ['JJA','DJF']:
    maps_dic = get_averaged_fields.get_maps_cases(cases,varl,startyear, endyear,
                                              avg_over_lev=avg_over_lev,
                                              pmin=pmin,
                                              pressure_adjust=pressure_adjust,
                                              p_level=p_level,
                                              time_mask=period)
    maps_dic = get_averaged_fields.get_maps_cases(cases_sec,varl_sec,startyear, endyear,
                                              avg_over_lev=avg_over_lev,
                                              pmin=pmin,
                                              pressure_adjust=pressure_adjust, p_level=p_level,
                                              time_mask=period)
