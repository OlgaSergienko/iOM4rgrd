import xarray as xr
from mom_in_stereo.generic_scalar import remap_scalar_to_MOM6


PROJSTRING = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#infile = "/work1/ovs/invAntelmer/BM2_grnd5km_MuBeta_nu.nc"
infile = "/archive/ovs/elmerAnt_inv/inv_alpha_crs.nc"
MOMgrids = "/work1/ovs/iOM4"
# === Load Data

ds = xr.open_dataset(infile)

# load MOM6 grid
hgrid = xr.open_dataset(f"{MOMgrids}/SP_025deg.nc")
hgrid = hgrid.set_coords(["x", "y"])

j60s=602

#hgrid_SO=hgrid[0:j60s]

#lon_model = hgrid["x"][0:j60s:2, 0::2]
#lat_model = hgrid["y"][0:j60s:2, 0::2]


# === Remap

remapped = remap_scalar_to_MOM6(ds, hgrid, projection=PROJSTRING)

# === Write to file

#remapped.to_netcdf(f"{MOMgrids}/inv_MuBetaNu_mom6_91W80S_5deg.nc")
remapped.to_netcdf(f"{MOMgrids}/inv_Mucrs_SP_025deg.nc")

# === Load Data

infile = "/archive/ovs/elmerAnt_inv/inv_beta_crs.nc"

# === Remap
ds = xr.open_dataset(infile)
remapped = remap_scalar_to_MOM6(ds, hgrid, projection=PROJSTRING)

# === Write to file
remapped.to_netcdf(f"{MOMgrids}/inv_Betacrs_SP_025deg.nc")
