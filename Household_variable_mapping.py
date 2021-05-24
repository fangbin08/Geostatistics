import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
# import openturns as ot
from scipy.interpolate import griddata
from sklearn.gaussian_process import GaussianProcessRegressor
from scipy.interpolate import Rbf
from sklearn.cluster import KMeans
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

path_gis_data = '/Users/binfang/Documents/SMAP_Project/data/gis_data'
path_hvm = '/Users/binfang/Documents/Household_variable_mapping'
df_csv_file = pd.read_csv(path_hvm + '/csv_files/KE_2014.csv', index_col=0)

# df_csv_file['assigned_lat'] = df_csv_file['assigned_lat'].apply(pd.to_numeric, errors='coerce')
# df_csv_file['assigned_long'] = df_csv_file['assigned_long'].apply(pd.to_numeric, errors='coerce')
# df_csv_file = df_csv_file.dropna(subset=['assigned_lat', 'assigned_long'])

var_floors = np.array(df_csv_file['floors'], dtype=float)
lat = np.array(df_csv_file['assigned_lat'], dtype=float)
lon = np.array(df_csv_file['assigned_long'], dtype=float)

# Group the variable by coordinates and calculate the accumulated values
coords = np.stack((lat, lon))
coords_unique, coords_ind = np.unique(coords, axis=1, return_index=True)
coords_list_group = [np.where((lat == coords_unique[0, x]) & (lon == coords_unique[1, x]))
                     for x in range(coords_unique.shape[1])]
var_floors_acc = np.array([np.sum(var_floors[coords_list_group[x]]) for x in range(len(coords_list_group))])
var_floors_acc[np.isnan(var_floors_acc)] = 0
lat_acc = coords_unique[0, :]
lon_acc = coords_unique[1, :]
var_floors_array = np.stack((lat_acc, lon_acc, var_floors_acc), axis=1)

columns = ['lat_acc', 'lon_acc', 'var_floors_acc']
df_floors_array = pd.DataFrame(var_floors_array, columns=columns)
df_floors_array.to_csv(path_hvm + '/csv_files_accum/df_floors_array.csv', index=True)


# # Plot the data with a scatter plot and a color map.
# fig = plt.figure()
# plt.scatter(lon_acc, lat_acc, c=var_floors_acc, cmap='viridis', s=1)
# plt.colorbar()
# plt.show()
#
# ########################################################################################################################
# inputDimension = 2
# coordinates = np.stack((var_floors_array[:, 1], var_floors_array[:, 0]), axis=1)
# observations = var_floors_array[:, -1].reshape(-1, 1)
# coordinates = ot.Sample(list(coordinates))
# observations = ot.Sample(observations.tolist())
# # input_train = ot.Sample(coordinates)
# # output_train = ot.Sample(observations, 1)
# basis = ot.ConstantBasisFactory(inputDimension).build()
# covariance_kernel = ot.SquaredExponential(1)
# algo = ot.KrigingAlgorithm(coordinates, observations, covariance_kernel, basis)
#
# # Fit
# algo.run()
# result = algo.getResult()
# krigingMetamodel = result.getMetaModel()
#
# # Create the 2D domain
# # kpInterval = ot.Interval([0., 0.], [1., 1.])
# kpInterval = ot.Interval([np.min(lon_acc), np.max(lat_acc)], [np.max(lon_acc), np.min(lat_acc)])
# # Define the number of interval in each direction of the box
# nx = 100
# ny = 100
# kpIndices = [nx-1, ny-1]
# kpMesher = ot.IntervalMesher(kpIndices)
# kpMeshBox = kpMesher.build(kpInterval)
#
# # Predict
# vertices = kpMeshBox.getVertices()
# predictions = krigingMetamodel(vertices)
# # Format for plot
# X = np.array(vertices[:,0]).reshape((ny,nx))
# Y = np.array(vertices[:,1]).reshape((ny,nx))
# predictions_array = np.array(predictions).reshape((ny,nx))
#
# # Plot
# fig = plt.figure()
# plt.pcolormesh(X, Y, predictions_array, shading='auto')
# plt.colorbar()
# plt.scatter(lon_acc, lat_acc, c=var_floors_acc, cmap='viridis', s=1)
# plt.show()


########################################################################################################################

# from sklearn.gaussian_process import GaussianProcessRegressor
#
# lons = np.linspace(np.min(lon_acc), np.max(lon_acc), 100)
# lats = np.linspace(np.max(lat_acc), np.min(lat_acc), 100)
# obs_x, obs_y = np.meshgrid(lons, lats)
#
# points = zip(obs_x,  obs_y)
# values = observations    # Replace with your observed data
#
# gp = GaussianProcessRegressor(kernel=None)
# gp.fit(points, values)
# XY_pairs = np.column_stack([obs_x.flatten(), obs_y.flatten()])
# predicted = gp.predict(XY_pairs).reshape(X.shape)


# df_floors = pd.read_csv('/Users/binfang/Documents/Household_variable_mapping/results/test.csv', index_col=0)
df_floors = pd.read_csv(path_hvm + '/csv_files_accum/df_floors_array.csv', index_col=0)


lon_pts = df_floors['lon_acc']
lat_pts = df_floors['lat_acc']

# Group the points by k-means clustering (True lat/lon)
n_clusters = len(lon_pts)//10
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(np.stack((lon_pts, lat_pts), axis=1))
coordinates_clusters_true = kmeans.cluster_centers_

lon_pts = (lon_pts - np.min(lon_pts))/(np.max(lon_pts) - np.min(lon_pts))
lat_pts = (lat_pts - np.min(lat_pts))/(np.max(lat_pts) - np.min(lat_pts))
coordinates = np.stack((lon_pts, lat_pts), axis=1)
observations = np.array(df_floors['var_floors_acc'], dtype=float).reshape(-1, 1)

# Group the points by k-means clustering
kmeans.fit(coordinates)
coordinates_clusters = kmeans.cluster_centers_
clusters = kmeans.predict(coordinates)
clusters_ind = [np.where(clusters == x) for x in range(n_clusters)]
observations_clusters = np.array([np.sum(observations[clusters_ind[x]]) for x in range(n_clusters)])

# plt.scatter(coordinates[:, 0], coordinates[:, 1], c=clusters, s=50, cmap='viridis')
# plt.scatter(coordinates_clusters[:, 0], coordinates_clusters[:, 1], c='black', s=200, alpha=0.5)

lon_grid = np.linspace(np.min(coordinates_clusters[:, 0]), np.max(coordinates_clusters[:, 0]), 1000)
lat_grid = np.linspace(np.max(coordinates_clusters[:, 1]), np.min(coordinates_clusters[:, 1]), 1000)
obs_x, obs_y = np.meshgrid(lon_grid, lat_grid)
obs_coordinates = np.stack((obs_x.flatten(), obs_y.flatten()), axis=1)

# Linear interpolation
grid_z0 = griddata(coordinates_clusters, observations_clusters, (obs_x, obs_y), method='linear', rescale=True)

fig = plt.figure()
plt.pcolormesh(obs_x, obs_y, grid_z0, shading='auto', cmap='Spectral_r')
plt.colorbar()
plt.scatter(coordinates_clusters[:, 0], coordinates_clusters[:, 1], c=observations_clusters, cmap='Spectral_r', edgecolors='k', s=7)
plt.show()


# Gaussian Process interpolation
gp = GaussianProcessRegressor(kernel=None, n_restarts_optimizer=9)
gp.fit(coordinates_clusters, observations_clusters)
grid_z1, sigma = gp.predict(obs_coordinates, return_std=True)
grid_z1 = np.reshape(grid_z1, (1000, 1000))
grid_z1[grid_z1 <= 0] = 0
grid_z1[grid_z1 > np.max(observations_clusters)] = np.nan

fig = plt.figure()
plt.pcolormesh(obs_x, obs_y, grid_z1, shading='auto', cmap='Spectral_r')
plt.colorbar()
plt.scatter(coordinates_clusters[:, 0], coordinates_clusters[:, 1], c=observations_clusters, cmap='Spectral_r', edgecolors='k', s=7)
plt.show()


# Radial Basis Function interpolation
rbf_adj = Rbf(coordinates_clusters[:, 0], coordinates_clusters[:, 1], observations_clusters, function='gaussian')
grid_z2 = rbf_adj(obs_x.flatten(), obs_y.flatten()).reshape(1000, 1000)
grid_z2[grid_z2 <= 0] = 0
grid_z2[grid_z2 > np.max(observations_clusters)] = np.nan

fig = plt.figure()
plt.pcolormesh(obs_x, obs_y, grid_z2, shading='auto', cmap='Spectral_r')
plt.colorbar()
plt.scatter(coordinates_clusters[:, 0], coordinates_clusters[:, 1], c=observations_clusters, cmap='Spectral_r', edgecolors='k', s=7)
plt.show()


# Georeferenced map
shape_country = ShapelyFeature(Reader(path_gis_data + '/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp').geometries(),
                                ccrs.PlateCarree(), edgecolor='black', facecolor='none')
shape_coast = ShapelyFeature(Reader(path_gis_data + '/gshhg-shp-2.3.7/WDBII_shp/i/WDBII_border_i_L1.shp').geometries(),
                                ccrs.PlateCarree(), edgecolor='black', facecolor='none')
extent = np.array([np.min(df_floors['lon_acc']), np.max(df_floors['lon_acc']),
                   np.min(df_floors['lat_acc']), np.max(df_floors['lat_acc'])])

fig = plt.figure(figsize=(6, 6), facecolor='w', edgecolor='k', dpi=150)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.add_feature(shape_country, linewidth=1)
ax.add_feature(shape_coast, linewidth=1)
img = ax.imshow(grid_z1, origin='upper', cmap='Spectral_r', transform=ccrs.PlateCarree(),
                extent=extent)
ax.scatter(coordinates_clusters_true[:, 0], coordinates_clusters_true[:, 1], c=observations_clusters, cmap='Spectral_r',
           edgecolors='k', s=7)
cbar = fig.colorbar(img, extend='both', pad=0.1)
gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.5, alpha=0.5, color='black')
gl.xlocator = mticker.MultipleLocator(base=2)
gl.ylocator = mticker.MultipleLocator(base=2)
gl.xlabel_style = {'size': 6}
gl.ylabel_style = {'size': 6}
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
plt.savefig(path_hvm + '/results/kenya.png')
plt.close()
