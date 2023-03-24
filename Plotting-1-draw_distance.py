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

Times = FontProperties(fname="Times.ttf", size=25) 
config = {
    "mathtext.fontset":'stix',
}
mpl.rcParams.update(config)
mpl.rcParams['axes.unicode_minus']=False

fig=plt.figure(figsize=(20,10))
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
            'name': 'dist_sim2r', # difference distance between real-end point  and simu-end point. 
            'range': [0, 135, 5],  # min, max, step
            'map_title': 'Distance between Argo and simulation\'s End-points ', 
            'unit': 'km',
            'isPercentage': '0',
            'fig_name': 'separation_distance'
        },# 0
        {
            'name': 'diff_dist', # difference distance between real and simu.
            'range': [-70000, 70000, 5000],  # min, max, step
            'map_title': 'Displacement differences between Argo and simulation', 
            'unit': 'm',
            'isPercentage': '0',
            'fig_name': 'diff_distance'
        }, # 1
        {
            'name': 'diff_dist_per', # difference distance percentage between real and simu.
            'range': [0, 1, 0.05],  # min, max, step
            'map_title': 'Percentage of displacement differences percent between Argo and simulation ', 
            'unit': '%',
            'isPercentage': '1',
            'fig_name': 'diff_dist_percent'
        },# 2
        # dist_real, dist_simu
        {
            'name': 'dist_real', 
            'range': [0, 180000, 5000],  # min, max, step
            'map_title': 'Argo Displacement in a cycle', 
            'unit': 'm',
            'isPercentage': '0',
            'fig_name': 'distance_real'
        },# 3
        {
            'name': 'dist_simu',
            'range': [0, 180000, 5000],  # min, max, step
            'map_title': 'Simulation Displacement in a cycle', 
            'unit': 'm',
            'isPercentage': '0',
            'fig_name': 'distance_simu'
        },# 4
    ]

ds = xr.open_dataset('file')
whichone = ds[info_map[current_ID]['name']] 

if info_map[current_ID]['isPercentage']==True:
    whichone = whichone*100
if info_map[current_ID]['unit']=='km':
    whichone = whichone/1000

lons = ds.lon.data
lats = ds.lat.data
lons = np.around(lons, 1) 
lats = np.around(lats, 1) 
dataarray = xr.DataArray(whichone.data.T, coords=[lats, lons], dims=['lats','lons'])
levels = np.arange(info_map[current_ID]['range'][0], info_map[current_ID]['range'][1], info_map[current_ID]['range'][2]) 
cycle_dataarray, cycle_lons = add_cyclic_point(dataarray, coord=lons)
cycle_lons, cycle_lats = np.meshgrid(cycle_lons, lats)
contourf = axe.contourf(cycle_lons, cycle_lats, cycle_dataarray, levels=levels, transform=ccrs.PlateCarree(), add_colorbar=False, cmap='') 
cax = fig.add_axes([axe.get_position().x1+0.03, axe.get_position().y0, 0.02, axe.get_position().y1-axe.get_position().y0-0.02 ]) 
cbar = fig.colorbar(contourf, cax=cax, drawedges=True, spacing='uniform', extend='neither') 
cbar.ax.tick_params(labelsize=20)  
cax.set_title(info_map[current_ID]['unit'], fontproperties=Times, fontsize=25, x=0.6)
plt.savefig()
plt.show()
