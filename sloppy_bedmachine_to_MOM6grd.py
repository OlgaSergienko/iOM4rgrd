import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
from cartopy import crs as ccrs

from sloppy.distrib import compute_block
from sloppy.distrib import compute_block_brute


def proj_xy(lon, lat, PROJSTRING):
    """ """
    from pyproj import CRS, Transformer

    # create the coordinate reference system
    crs = CRS.from_proj4(PROJSTRING)
    # create the projection from lon/lat to x/y
    proj = Transformer.from_crs(crs.geodetic_crs, crs)
    # compute the lon/lat
    xx, yy = proj.transform(lon, lat, direction="FORWARD")
    return xx, yy


PROJSTRING = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#---------------- bedmachine + reduction
bedmachine = xr.open_dataset(
    "/net2/rnd/BedMachineAntarctica_2020-07-15_v02.nc"
)


xx_bm_full, yy_bm_full = np.meshgrid(bedmachine["x"].values, bedmachine["y"].values)

#----------------- read lon/lat MOM6 grid corners until 60S

#gridname = "73E82S_025deg"
gridname = "SP_025deg"
hgrid = xr.open_dataset(f"/work/ovs/iOM4/{gridname}.nc")
iOM4dir = "/work/ovs/iOM4/"
j60s=602

lon_model = hgrid["x"][0:j60s:2, 0::2]
lat_model = hgrid["y"][0:j60s:2, 0::2]

xx_model, yy_model = proj_xy(lon_model, lat_model, PROJSTRING)



### remapping

out = xr.Dataset()

mask = compute_block(
    xx_model,
    yy_model,
    bedmachine["mask"].values,
    xx_bm_full,
    yy_bm_full,
    is_stereo=False,
    is_carth=True,
    PROJSTRING=PROJSTRING,
    residual=True,
)

out = xr.Dataset()
out["mask"] = xr.DataArray(data=firn[0, :, :], dims=("y", "x"))
out.to_netcdf(f"{iOM4dir}/mask_bedmachine_remapped_iOM4_{gridname}.nc")


thk = compute_block(
    xx_model,
    yy_model,
    bedmachine["thickness"].values,
    xx_bm_full,
    yy_bm_full,
    is_stereo=False,
    is_carth=True,
    PROJSTRING=PROJSTRING,
    residual=True,
)

out["thickness"] = xr.DataArray(data=thk[0, :, :], dims=("y", "x"))
out.to_netcdf(f"thickness_bedmachine_remapped_iOM4_{gridname}.nc")

firn = compute_block(
    xx_model,
    yy_model,
    bedmachine["firn"].values,
    xx_bm_full,
    yy_bm_full,
    is_stereo=False,
    is_carth=True,
    PROJSTRING=PROJSTRING,
    residual=True,
)

out = xr.Dataset()
out["firn"] = xr.DataArray(data=firn[0, :, :], dims=("y", "x"))
out.to_netcdf(f"firn_bedmachine_remapped_iOM4_{gridname}.nc")

surf = compute_block(
    xx_model,
    yy_model,
    bedmachine["surface"].values,
    xx_bm_full,
    yy_bm_full,
    is_stereo=False,
    is_carth=True,
    PROJSTRING=PROJSTRING,
    residual=True,
)

out = xr.Dataset()
out["surface"] = xr.DataArray(data=surf[0, :, :], dims=("y", "x"))
out.to_netcdf(f"surface_bedmachine_remapped_iOM4_{gridname}.nc")


bed = compute_block(
    xx_model,
    yy_model,
    bedmachine["bed"].values,
    xx_bm_full,
    yy_bm_full,
    is_stereo=False,
    is_carth=True,
    PROJSTRING=PROJSTRING,
    residual=True,
)

out = xr.Dataset()
out["bed"] = xr.DataArray(data=bed[0, :, :], dims=("y", "x"))
out["h2"] = xr.DataArray(data=bed[3, :, :], dims=("y", "x"))

out.to_netcdf(f"bed_bedmachine_remapped_iOM4_{gridname}.nc")

#
#
#
#
#
##plt.figure()
##plt.pcolormesh(xx_model, yy_model, thk10x[0, :, :], vmax=5000)
##plt.colorbar()
##plt.title("10x downsampled - thickness")
 #
 #plt.show()
#
##subplot_kws = dict(
##    projection=ccrs.SouthPolarStereo(central_longitude=0.0), facecolor="grey"
##)
##plt.figure(figsize=[10, 8])
##ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=0.0))
### ax.stock_img()
##plt.pcolormesh(
##    lon_model,
##    lat_model,
##    thk5x[0, :, :],
##    shading="auto",
##    cmap="jet",
##    transform=ccrs.PlateCarree(),
##), plt.clim(-0, 4000), plt.colorbar()
##ax.set_extent([-180, 180, -55, -90], ccrs.PlateCarree())
##ax.gridlines(color="black", alpha=0.5, linestyle="--")
##
##plt.show()
#
#
