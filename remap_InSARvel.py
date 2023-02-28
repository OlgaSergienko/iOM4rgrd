import xarray as xr
from mom_in_stereo.MAR import remap_MAR_velocities_to_MOM6
from mom_in_stereo.coords import add_lon_lat
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

MOMgrids = "/work1/ovs/iOM4"
gridname = "SP_025deg"

# === Load Data


#velA=xr.open_dataset('/work1/ovs/invAntelmer/BM2_grnd5km_MuBeta.nc')
velA=xr.open_dataset('/work1/ovs/Measures_vel/vel_Ant_fld_cntr_Peninsfull.nc')


# load MOM6 grid
#hgrid = xr.open_dataset(f"{MOMgrids}/sosc91W80S_05deg.nc")
#hgrid = xr.open_dataset(f"{MOMgrids}/{gridname}.nc")
hgridSC = xr.open_dataset(f"{MOMgrids}/{gridname}.ncSC.nc")
hgridSO = xr.open_dataset(f"{MOMgrids}/{gridname}.ncSO.nc")
hgridMerc = xr.open_dataset(f"{MOMgrids}/{gridname}.ncMerc.nc")
hgrid =  xr.Dataset()
hgrid["y"]=xr.concat([hgridSC.y,hgridSO.y,hgridMerc.y[:98,:]],dim="nyp")
hgrid["x"]=xr.concat([hgridSC.x,hgridSO.x,hgridMerc.x[:98,:]],dim="nyp")
hgrid["angle_dx"]=xr.concat([hgridSC.angle_dx,hgridSO.angle_dx,hgridMerc.angle_dx[:98,:]],dim="nyp")
hgrid = hgrid.set_coords(["x", "y"])
# === Remap velocities

uvmar_remapped = remap_MAR_velocities_to_MOM6(velA.isel(x=slice(0,-1,10),y=slice(0,-1,10)), hgrid,u="VX",v="VY",units_grid="m")
uvmar_remapped.to_netcdf(f"{MOMgrids}/INPUT/Mesures_vel_mom6_{gridname}.nc4",format='NETCDF3_CLASSIC')

velA=xr.open_dataset('/work1/ovs/invAntelmer/BM2_grnd5km_MuBeta.nc')
uvmar_remapped = remap_MAR_velocities_to_MOM6(velA, hgrid,u="uobs 1",v="uobs 2",units_grid="m")
uvmar_remapped.to_netcdf(f"{MOMgrids}/INPUT/obs_vel_mom6_{gridname}.nc4",format='NETCDF3_CLASSIC')


# === Plot results

## create a coarser dataset
#uvmar_coarse = (
#    velA.isel(x=slice(0, -1), y=slice(0, -1))
#    .coarsen(x=100, y=100,boundary="pad")
#    .mean()
#)
#uvmar_remapped_coarse = (
#    uvmar_remapped.isel(xh=slice(0, -1), yh=slice(0, -1)).coarsen(xh=10, yh=2,boundary="pad").mean()
#)
#
#
## quiver plot original data
#plt.figure(figsize=[12, 8])
##plt.contour(uvmar["x"], uvmar["y"], uvmar["GROUND"], [99.0], colors="k")
#plt.quiver(uvmar_coarse["x"], uvmar_coarse["y"], uvmar_coarse["VX"], uvmar_coarse["VY"])
#
#
## quiver plot in MOM6 coords
#plt.figure(figsize=[12, 8])
#spproj = ccrs.SouthPolarStereo(central_longitude=0.0)
#ax = plt.axes(projection=spproj)
#uvmar_remapped_coarse.plot.quiver(
#    "lon", "lat", "VX", "VY", ax=ax, transform=ccrs.PlateCarree()
#)
#plt.contour(hgrid["x"],hgrid["y"],hgrid["angle_dx"],levels=[-20,0,20],colors="r")
#ax.coastlines()
#ax.set_extent([-180, 180, -50, -90], ccrs.PlateCarree())
#
#
#plt.show()
