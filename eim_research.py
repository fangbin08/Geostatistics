import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as cb
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import glob
import os
import gdal
import rasterio
from rasterio.transform import Affine
from rasterio.crs import CRS
import fiona
import scipy.ndimage
# Ignore runtime warning
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

path_gis_data = '/Users/binfang/Documents/SMAP_Project/data/gis_data'
path_hvm = '/Users/binfang/Documents/Household_variable_mapping'
path_eim = '/Users/binfang/Documents/EIM_project'

########################################################################################################################
# 1. Process tickborne data
# Lyme
df_tickborne_lyme = pd.read_excel(path_eim + '/Engineering in medicine - tickborne illness data.xlsx', index_col=0,
                             sheet_name='Lyme_disease')
df_tickborne_lyme = df_tickborne_lyme.dropna(axis='index', how='any')
df_tickborne_lyme['PostalCode'] = df_tickborne_lyme['PostalCode'].astype(str)
# Remove digits after bar symbol
df_tickborne_lyme['PostalCode'] = df_tickborne_lyme['PostalCode'].apply(lambda x: x.split('-')[0])
df_tickborne_lyme['PostalCode'] = df_tickborne_lyme['PostalCode'].apply(lambda x: x[0:5])
# df_tickborne_lyme['PostalCode'] = df_tickborne_lyme['PostalCode'].astype(int)
df_tickborne_lyme = df_tickborne_lyme.sort_values(by='PostalCode', ascending=True)
df_tickborne_lyme_sum = pd.DataFrame(df_tickborne_lyme['PostalCode'].value_counts(ascending=True))

# Ehrlichiosis
df_tickborne_ehr = pd.read_excel(path_eim + '/Engineering in medicine - tickborne illness data.xlsx', index_col=0,
                             sheet_name='Ehrlichiosis')
df_tickborne_ehr = df_tickborne_ehr.dropna(axis='index', how='any')
df_tickborne_ehr['PostalCode'] = df_tickborne_ehr['PostalCode'].astype(str)
# Remove digits after bar symbol
df_tickborne_ehr['PostalCode'] = df_tickborne_ehr['PostalCode'].apply(lambda x: x.split('-')[0])
df_tickborne_ehr_sum = pd.DataFrame(df_tickborne_ehr['PostalCode'].value_counts(ascending=True))

# Rickettsia
df_tickborne_ric = pd.read_excel(path_eim + '/Engineering in medicine - tickborne illness data.xlsx', index_col=0,
                             sheet_name='Rickettsia')
df_tickborne_ric = df_tickborne_ric.dropna(axis='index', how='any')
df_tickborne_ric['PostalCode'] = df_tickborne_ric['PostalCode'].astype(str)
# Remove digits after bar symbol
df_tickborne_ric['PostalCode'] = df_tickborne_ric['PostalCode'].apply(lambda x: x.split('-')[0])
df_tickborne_ric_sum = pd.DataFrame(df_tickborne_ric['PostalCode'].value_counts(ascending=True))

df_tickborne_sum = df_tickborne_lyme_sum.join(df_tickborne_ehr_sum, how='left', lsuffix='_lyme', rsuffix='_ehr')
df_tickborne_sum = df_tickborne_sum.join(df_tickborne_ric_sum, how='left', rsuffix='_ric')
df_tickborne_sum = df_tickborne_sum.rename({'PostalCode': 'PostalCode_ric'}, axis=1)
df_tickborne_sum.reset_index(inplace=True)
df_tickborne_sum = df_tickborne_sum.rename({'index': 'zipcode_tickborne'}, axis=1)


# 2. Process GIS shapefile data
va_shape = fiona.open(path_gis_data + '/va_zipcode/va_zipcode.shp')
zipcode_lat = [feature['properties']['INTPTLAT10'] for feature in va_shape]
zipcode_lon = [feature['properties']['INTPTLON10'] for feature in va_shape]
zipcode = [feature['properties']['ZCTA5CE10'] for feature in va_shape]

df_va_tickborne = pd.DataFrame(list(zip(zipcode, zipcode_lat, zipcode_lon)),
                               columns=['zipcode', 'zipcode_lat', 'zipcode_lon'])
df_va_tickborne['zipcode_lat'] = df_va_tickborne['zipcode_lat'].astype(float)
df_va_tickborne['zipcode_lon'] = df_va_tickborne['zipcode_lon'].astype(float)

# Join the tickborne table to VA zipcode table
df_va_tickborne = df_va_tickborne.join(df_tickborne_sum.set_index('zipcode_tickborne'), on='zipcode')
df_va_tickborne = df_va_tickborne.fillna(0)
df_va_tickborne = df_va_tickborne.rename({'PostalCode_lyme': 'Lyme_disease',
                                          'PostalCode_ehr': 'Ehrlichiosis',
                                          'PostalCode_ric': 'Rickettsia'}, axis=1)
df_va_tickborne = df_va_tickborne.astype({'Lyme_disease': int, 'Ehrlichiosis': int, 'Rickettsia': int})

df_va_tickborne.to_csv(path_eim + '/df_va_tickborne.csv')


# 3. Draw maps
df_va_tickborne = pd.read_csv(path_eim + '/df_va_tickborne.csv', index_col=0)
va_shape = fiona.open(path_gis_data + '/va_zipcode/va_zipcode.shp')
shp_extent = list(va_shape.bounds)
shp_extent = list((shp_extent[0], shp_extent[2], shp_extent[1], shp_extent[3]))
shp_reader = shpreader.Reader(path_gis_data + '/va_zipcode/va_zipcode.shp')

cmap = cm.get_cmap('hot_r')
zipcode_values = df_va_tickborne['zipcode'].astype(str)

dise_list = ['Lyme_disease', 'Ehrlichiosis', 'Rickettsia']
df_values_all = []
sm_all = []
norm_all = []
for ipt in range(3):
    values = df_va_tickborne[dise_list[ipt]]
    df_values = values
    norm = colors.Normalize(vmin=values.min(), vmax=values.max())
    # df_values = norm(values).data
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(values)
    df_values_all.append(df_values)
    sm_all.append(sm)
    norm_all.append(norm)
    del(values, df_values, sm, norm)


# Draw subplot maps
fig = plt.figure(figsize=(5, 5), facecolor='w', edgecolor='k', dpi=300)
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95, hspace=0.2, wspace=0.2)
for ipt in range(3):
    ax = fig.add_subplot(3, 1, ipt+1, projection=ccrs.PlateCarree())
    ax.set_extent(shp_extent, ccrs.PlateCarree())
    for record, zip_poly in zip(shp_reader.records(), shp_reader.geometries()):
        color = df_values_all[ipt][int(np.where(zipcode_values == record.attributes['ZCTA5CE10'])[0])]
        ax.add_geometries([zip_poly], ccrs.PlateCarree(), facecolor=cmap(norm_all[ipt](color)), edgecolor='black', linewidth=0.1)
        ax.text(-83.63, 39.23, dise_list[ipt], fontsize=7, horizontalalignment='left',
                verticalalignment='top', weight='bold')

    cbar = plt.colorbar(sm_all[ipt], extend='both', orientation='vertical', pad=0.1)
    # cbar.set_ticks([])
    # cbar = fig.colorbar(sm, cax=ax, extend='both', orientation='horizontal', pad=0.1)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.locator_params(nbins=6)
    cbar.set_label('Count', fontsize=7, x=2, y=1.10, labelpad=-15, rotation=0)

plt.savefig(path_eim + '/results/tickborne_distribution.png')
plt.close()

# df_va_tickborne.plot('zipcode', y=['Lyme_disease', 'Ehrlichiosis', 'Rickettsia'], kind="line", figsize=(5, 5))
