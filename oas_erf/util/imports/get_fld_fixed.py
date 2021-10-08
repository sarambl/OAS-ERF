import os

import xarray as xr

from oas_erf import constants
from oas_erf.util.filenames import get_filename_ng_field
from oas_erf.util.imports.fix_xa_dataset_v2 import xr_fix
from oas_erf.util.imports.get_pressure_coord_fields import get_pressure_coord_fields
from oas_erf.util.imports.import_fields_xr_v2 import xr_import_NorESM, import_constants


# %%

def get_field_fixed(case, varlist, from_time, to_time, raw_data_path=constants.get_input_datapath(),
                    pressure_adjust=True, model = 'NorESM', history_fld='.h0.', comp='atm', chunks=None,
                    get_constants=True):
    """
    Imports and fixes:
    :param model:
    :param case:
    :param varlist:
    :param from_time:
    :param to_time:
    :param raw_data_path:
    :param pressure_adjust:
    :return:
    """

    # If pressure coordinate --> check if in pressure coordinates, else get_pressure_coordinate etc
    # IF LEV NOT DIM! JUST LOAD AND AVERAGE
    # If not pressure coordinates --> check outpaths['original_coords']=path_outdata + '/computed_fields_ng'
    #raw_data_path=constants.get_input_datapath()
    #pressure_adjust=True; model = 'NorESM'; history_fld='.h0.'; comp='atm'; chunks=None
    if type(from_time) is int: from_time='%s'%from_time
    if type(to_time) is int: to_time='%s'%to_time
    if len(to_time)==4:
        to_time = to_time+'-12'
    if len(from_time)==4: # if only year, add month
        from_time=from_time+'-01'


    if pressure_adjust:
        #log.ger.debug('Calling get_pressure_coord_fields...')
        ds = get_pressure_coord_fields(case,
                                       varlist,
                                       from_time,
                                       to_time,
                                       history_fld,
                                       comp=comp,
                                       model=model)
        if get_constants:
            _ds = import_constants(case)
            ds = xr.merge([ds, _ds])
        return ds
    else:
        if varlist is not None:
            fl = []
            vl_lacking = []
            for var in varlist:
                fn = get_filename_ng_field(var, model, case, from_time, to_time)
                if os.path.isfile(fn):
                    fl.append(fn)
                else:
                    vl_lacking.append(var)
        else:
            vl_lacking=varlist
        # print(vl_lacking)
        ds = xr_import_NorESM(case, vl_lacking, from_time, to_time, path=raw_data_path,
                              model=model,
                              history_fld=history_fld,
                              comp=comp, chunks=chunks)
        ds = xr_fix(ds, model_name=model)
        if len(fl)>0:
            ds_f_file = xr.open_mfdataset(fl, combine='by_coords')
            ds = xr.merge([ds, ds_f_file])
        return ds