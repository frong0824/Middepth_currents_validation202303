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
            'name': 'diff_v', # difference velocity between real and simu.
            'range': [-10, 10.5, 0.5],  # min, max, step
            'map_title': 'Velocity differences between Argo and simulation', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'diff_velocity'
        }, # 0
        {
            'name': 'diff_vecu', # difference velocity between real and simu.
            'range': [-5, 5, 0.5],  # min, max, step
            'map_title': 'U-Velocity differences between Argo and simulation', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'diff_vecu'
        }, # 1
        {
            'name': 'diff_vecv', # difference velocity between real and simu.
            'range': [-5, 5, 0.5],  # min, max, step
            'map_title': 'V-Velocity differences between Argo and simulation', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'diff_vecv'
        }, # 2
        {
            'name': 'diff_v_per', # # difference velocity percentage between real and simu. 
            'range': [-100, 105, 5],  # min, max, step
            'map_title': 'Percentage of velocity differences percent between Argo and simulation ', 
            'unit': '%',
            'isPercentage': '1',
            'fig_name': 'diff_velocity_percent'
        }, # 3
        {
            'name': 'diff_v_per', # # difference velocity percentage between real and simu. 
            # 'range': [0, 105, 5],  
            'range': [-20, 105, 5],  # only positive
            # 'range': [-120, 50, 5],  # only negative
            'map_title': 'Percentage of velocity differences percent between Argo and simulation ', 
            'unit': '%',
            'isPercentage': '1',
            'fig_name': 'diff_velocity_percent_positive' # diff_velocity_percent_negative
        }, # 4
        # v_real, v_simu, vecu_real, vecu_simu, vecv_real, vecv_simu
        {
            'name': 'v_real', # difference velocity between real and simu.
            'range': [0, 18.5, 0.5],  # min, max, step
            'map_title': 'Argo Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'v_real'
        }, # 5
        {
            'name': 'v_simu', # difference velocity between real and simu.
            'range': [0, 18.5, 0.5],  # min, max, step
            'map_title': 'Simulation Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'v_simu'
        }, # 6
        {
            'name': 'vecu_real', # difference velocity between real and simu.
            'range': [-7, 15, 1],  # min, max, step
            'map_title': 'Argo U-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'vecu_real'
        }, # 7
        {
            'name': 'vecu_simu', # difference velocity between real and simu.
            'range': [-7, 15, 1],  # min, max, step
            'map_title': 'Simulation U-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'vecu_simu'
        }, # 8
        {
            'name': 'vecv_real', # difference velocity between real and simu.
            'range': [-10, 10, 0.5],  # min, max, step
            'map_title': 'Argo V-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'vecv_real'
        }, # 9
        {
            'name': 'vecv_simu', # difference velocity between real and simu.
            'range': [-10, 10, 0.5],  # min, max, step
            'map_title': 'Simulation V-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '0',
            'fig_name': 'vecv_simu'
        }, # 10
        {
            'name': 'vecu_per', # difference velocity between real and simu.
            'range': [-1, 1, 0.1],  # min, max, step
            'map_title': 'Simulation V-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '1',
            'fig_name': 'vecv_simu'
        }, # 11
        {
            'name': 'vecv_per', # difference velocity between real and simu.
            'range': [-1, 1, 0.1],  # min, max, step
            'map_title': 'Simulation V-Velocity', 
            'unit': 'cm s$^{-1}$',
            'isPercentage': '1',
            'fig_name': 'vecv_simu'
        }, # 12

    ]

ds = xr.open_dataset('file')
whichone = ds[info_map[current_ID]['name']] 
if info_map[current_ID]['isPercentage']== '1' :
    whichone = whichone*100 # 100%

lons = ds.lon.data
lats = ds.lat.data
lons = np.around(lons, 1) 
lats = np.around(lats, 1)

dataarray = xr.DataArray(whichone.data.T, coords=[lats, lons], dims=['lats','lons'])
levels = np.arange(info_map[current_ID]['range'][0], info_map[current_ID]['range'][1], info_map[current_ID]['range'][2])  
cycle_dataarray, cycle_lons = add_cyclic_point(dataarray, coord=lons)
cycle_lons, cycle_lats = np.meshgrid(cycle_lons, lats)
contourf = axe.contourf(cycle_lons, cycle_lats, cycle_dataarray, levels=levels, transform=ccrs.PlateCarree(), add_colorbar=False, cmap='')  
cax = fig.add_axes([axe.get_position().x1+0.03, axe.get_position().y0, 0.02, axe.get_position().y1-axe.get_position().y0-0.02]) 
cbar = fig.colorbar(contourf, cax=cax, drawedges=True, spacing='uniform', extend='neither') 
cbar.ax.tick_params(labelsize=20) 
cax.set_title(info_map[current_ID]['unit'], fontproperties=Times, fontsize=25, x=0.6)

plt.savefig()
plt.show()
