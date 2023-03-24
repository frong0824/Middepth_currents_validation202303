import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from matplotlib.font_manager import FontProperties
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.util import add_cyclic_point
import matplotlib.ticker as mticker
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import xarray as xr
from time import strftime, localtime
import math
from pyproj import Geod


Times = FontProperties(fname="Times.ttf", size=25) 
config = {
    "mathtext.fontset":'stix',
}
mpl.rcParams.update(config)
mpl.rcParams['axes.unicode_minus']=False

fig=plt.figure(figsize=(60,30))
axe=plt.subplot(1,1,1,projection=ccrs.PlateCarree(-180))
axe.set_extent([-179.999, 179.999, -80, 80], crs=ccrs.PlateCarree())
axe.add_feature(cfeat.LAND, linewidth=0.1,color='grey')
axe.set_xticks(np.arange(-180.0, 180.0, 60), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(-60, 80, 30), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k',length=8)
labels = axe.get_xticklabels() + axe.get_yticklabels()
[label.set_fontproperties(Times) for label in labels]

current_ID = 0
scale_str = '1dg' 
info_map = [
        {
            'name': 'diff_azim', # difference azimuth between real and simu.
            'range': [0, 150, 10],  # min, max, step
            'map_title': 'Azimuth differences between Argo and simulation', 
            'unit': ' ° ',
            'isPercentage': '0',
            'fig_name': 'diff_azimuth'
        },
        #azim_real, azim_simu
        {
            'name': 'azim_real', # difference azimuth between real and simu.
            'range': [-180-22.5, 180+22.5+45, 45],  # min, max, step
            'map_title': 'Argo Azimuths', 
            'unit': ' ° ',
            'isPercentage': '0',
            'fig_name': 'azim_real'
        },
        {
            'name': 'azim_simu', # difference azimuth between real and simu.
            'range': [-180-22.5, 180+22.5+45, 45],  # min, max, step
            'map_title': 'Simulation Azimuths', 
            'unit': ' ° ',
            'isPercentage': '0',
            'fig_name': 'azim_simu'
        }
    ]


ds = xr.open_dataset('file')
lons = ds.lon.data
lats = ds.lat.data
lons = np.around(lons, 1) 
lats = np.around(lats, 1)

g = Geod(ellps='WGS84')
if info_map[current_ID]['name'] == 'azim_real':
    sin_real = ds.sin_real.values
    cos_real = ds.cos_real.values
    df = ds.to_dataframe()
    df['azim_real'] = df.apply(lambda row:g.inv(0.0, 0.0, row.sin_real, row.cos_real)[0], axis=1)  
    ds = df.to_xarray()
    whichone = ds.azim_real

elif info_map[current_ID]['name'] == 'azim_simu':
    sin_simu = ds.sin_simu.values
    cos_simu = ds.cos_simu.values
    df = ds.to_dataframe()
    df['azim_simu'] = df.apply(lambda row:g.inv(0.0, 0.0, row.sin_simu, row.cos_simu)[0], axis=1)  
    ds = df.to_xarray()
    whichone = ds.azim_simu
else:
    whichone = ds[info_map[current_ID]['name']]   

dataarray = xr.DataArray(whichone.data.T, coords=[lats, lons], dims=['lats','lons'])

levels = np.arange(info_map[current_ID]['range'][0], info_map[current_ID]['range'][1], info_map[current_ID]['range'][2])
cycle_dataarray, cycle_lons = add_cyclic_point(dataarray, coord=lons)
cycle_lons, cycle_lats = np.meshgrid(cycle_lons, lats)
contourf = axe.contourf(cycle_lons, cycle_lats, cycle_dataarray, levels=levels, transform=ccrs.PlateCarree(), add_colorbar=False, cmap='')
cax = fig.add_axes([axe.get_position().x1+0.03, axe.get_position().y0, 0.02, axe.get_position().y1-axe.get_position().y0-0.02]) 
cbar = fig.colorbar(contourf, cax=cax, drawedges=True, spacing='uniform', extend='neither') 
cbar.ax.tick_params(labelsize=20) 
cbar.set_ticks(levels)
plt.savefig()
plt.show()
