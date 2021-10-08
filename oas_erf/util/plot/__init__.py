import xarray as xr


def calc_vmin_vmax(list_xr, quant=0.01):
    """
    Calculates vmin/vmax for list of xr DataArrays. Convenient for plotting.
    :param list_xr: list of xarray DataArrays.
    :param quant: quantile to calculate vmin/vmax from
    :return: vmin, vmax
    """
    max_v = -9e99
    min_v = 9e99
    for xa in list_xr:
        xa:xr.DataArray
        xa.load()
        max_v = max(max_v, xa.compute().quantile((1 - quant), keep_attrs=False))
        min_v = min(min_v, xa.compute().quantile(quant, keep_attrs=False))
    return min_v, max_v
