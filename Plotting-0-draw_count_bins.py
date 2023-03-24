import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.util import add_cyclic_point
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
axe.add_feature(cfeat.LAND.with_scale('50m'), color='grey')
axe.set_extent([-179.999, 179.999, -80, 80], crs=ccrs.PlateCarree())
axe.set_xticks(np.arange(-180.0, 180.0, 60), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(-60, 80, 30), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())

axe.tick_params(labelcolor='k',length=8)
labels = axe.get_xticklabels() + axe.get_yticklabels()
[label.set_fontproperties(Times) for label in labels]


ds = xr.open_dataset('file')
t = ds['Join_Count']
ele = 30
lons = ds.lon.data
lats = ds.lat.data
dataarray = xr.DataArray(t.data.T, coords=[lats, lons], dims=['lat','lon'])
dataarray = dataarray.where(dataarray > ele) 

levels = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
mesh = dataarray.plot.pcolormesh(ax=axe, levels=levels, transform=ccrs.PlateCarree(), add_colorbar=False, cmap='', hatch="/", zorder=1)
mpl.rcParams['hatch.linewidth'] = 0.1
cax = fig.add_axes([axe.get_position().x1+0.03, axe.get_position().y0, 0.02, axe.get_position().y1-axe.get_position().y0-0.02]) 
cbar = fig.colorbar(mesh, cax=cax, drawedges=True, spacing='uniform', extend='both')
cbar.ax.tick_params(labelsize=20) 
cbar.set_ticks(levels)

plt.savefig() 
plt.show()

