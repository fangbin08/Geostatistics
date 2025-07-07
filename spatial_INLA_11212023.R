###
library(tidyr)
library(tidyverse)
library(dplyr)
library(readxl)
library(dplyr)
library(sp)
library(sf)
library(inlabru)
library(INLA)
library(parallel)
library(readr)
library(magrittr)
###

path_hvm = '/Users/binfang/Documents/Household_variable_mapping/csv_files_accum/'
path_hvm_raster = '/Volumes/Elements/Datasets/Household_variable_mapping/csv_raster_data/'
path_csv = '/Users/binfang/Documents/Household_variable_mapping/csv_files/'
path_csv_cat = '/Users/binfang/Documents/Household_variable_mapping/categorical_data/'
path_pred = '/Users/binfang/Documents/Household_variable_mapping/predictions/'
path_fi = '/Users/binfang/Documents/Household_variable_mapping/results/'


# Import the data
df_ext_obs_all <- read.csv(file = paste0(path_hvm, 'df_ext_obs_all_cat.csv'), na.strings = "")
df_ext_obs_all <- df_ext_obs_all[, -grep('X', colnames(df_ext_obs_all))]


df_csv_file <- read.csv(file = paste0(path_hvm, 'df_cat_nonzero.csv'), na.strings = "")
df_csv_file <- subset(df_csv_file, select = -c(DATE))
df_raster <- read.csv(file = paste0(path_hvm, 'df_raster_all_nonnan_cat.csv'),
                      row.names = 1, na.strings = "")
df_raster$urban <- as.factor(df_raster$urban)
df_raster$climate <- as.factor(df_raster$climate)
df_raster$region <- as.factor(df_raster$region)
df_raster_time <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(df_raster))), c('time'))
df_raster_time[is.na(df_raster_time)] = 8.7908778
df_raster <- data.frame(cbind(df_raster_time, df_raster))
df_ext_obs_no_wwtp <- subset(df_ext_obs_all, select = -c(wwtp_dist, pigs2, poultry2, ruminants2))


# Process model input cluster data
target_var <- df_csv_file$Pigs
cluster <- df_csv_file$Cluster
time <- df_csv_file$Log_Time
Latitude <- df_csv_file$Latitude
Longitude <- df_csv_file$Longitude


# Generate model input data for brm model
# brm_model <- data.frame(cbind(cluster, target_var, time, Latitude, Longitude,
#                               df_ext_obs_no_wwtp)) # Type 1

# brm_model <- data.frame(cbind(cluster, target_var, time, Latitude, Longitude,
#                               df_ext_obs_no_wwtp, df_ext_obs_all$wwtp_dist))
# colnames(brm_model)[ncol(brm_model)] <- 'wwtp_dist' # Type 2

brm_model <- data.frame(cbind(cluster, target_var, time, Latitude, Longitude,
                              df_ext_obs_no_wwtp, df_ext_obs_all$pigs2))
colnames(brm_model)[ncol(brm_model)] <- 'pigs2' # Type 3

# brm_model <- data.frame(cbind(cluster, target_var, time, Latitude, Longitude,
#                               df_ext_obs_no_wwtp))

brm_model_cplt <- brm_model[complete.cases(brm_model), ]
brm_model_cplt$urban <- as.factor(brm_model_cplt$urban)
brm_model_cplt$land <- as.factor(brm_model_cplt$land)
brm_model_cplt$climate <- as.factor(brm_model_cplt$climate)
brm_model_cplt$region <- as.factor(brm_model_cplt$region)

# set.seed(1001)
# brm_model_cplt_sample <- brm_model_cplt %>% group_by(cluster) %>% slice_sample(n=1)

write.csv(brm_model_cplt, paste0(path_pred, '/model_pigs.csv'), 
          row.names = FALSE)

write.csv(df_raster, paste0(path_pred, '/df_raster.csv'), row.names = FALSE)


####################################################################################
# Build the INLA model
# read in data #
csv_data <- read.csv(file = paste0(path_pred, 'brm_model_poultry.csv'), na.strings = "")
glimpse(csv_data)
## remake categorical variables 
csv_data$climate <- as.factor(csv_data$climate)
csv_data$land <- as.factor(csv_data$land)
csv_data$region <- as.factor(csv_data$region)
csv_data$urban <- as.factor(csv_data$urban)

# can pilot using a sample of data
#set.seed(10111985)
#d<-sample_n(csv_data,20000)

# for simplicity just renamed the csv_data data frame
class_name <- unique(csv_data$target_var)[2]
d<-csv_data
N<-nrow(d)
d$Y_bin<-as.numeric(d$target_var == class_name)

#** set coordinate system and reproject into meters **#
coordinates(d)=~Longitude+Latitude
#initial coordinate system
proj4string(d)= CRS("+proj=longlat +ellps=sphere")  
#project into Mollweide with KM as units
D_proj <- spTransform(d, CRS("+proj=moll +units=km")) %>% 
          st_as_sf() 

# make dummy vars for factors  (required by inlabru)
model.matrix(~climate, d) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> clim
model.matrix(~land, d) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> land
model.matrix(~region, d) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> reg
model.matrix(~urban, d) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> urb
## append to original (use D from now on)
D<-cbind(d,clim,land,reg,urb)


############################################################################################
#** Spatial binomial using R-INLA **#
##

###########################################
#** make SPDE mesh for INLA and inlabru **#
# continents as boundary #
library(rworldmap)
data("coastsCoarse")
bound<-coastsCoarse %>% 
       spTransform(CRS("+proj=moll +units=km")) %>%
       inla.sp2segment()
plot(bound)

# pull the projected coordinates (in KM)
coo <- st_coordinates(D_proj)
max_edge <- diff(range(coo[,1]))/(180) #<-was 100
max_edge # should be much less than range
mesh <- inla.mesh.2d(boundary = bound, #<-boundary of observed data
  loc = coo, 
  max.edge = c(max_edge,50*max_edge),
  cutoff = 0.20*max_edge) 
##
mesh$n #<-projects 250K locations to 18K locations
plot(mesh)
points(coo,col='red',cex=0.1)

# set priors #
matern_spde <-
  inla.spde2.pcmatern(mesh,
                    # set priors for spatial SD and range
                    prior.sigma = c(3.0,0.10), #P(SD>3.0) = 0.1
                    prior.range = c(100,0.10) #P(range<200km) = 0.1
  )


gc()
# options(expressions = 500000)
###
fit_sp<-bru(Y_bin ~ Intercept(1) + time +
           # dummy vars for all categorical variables #
           climate2 + climate3 + climate4 + climate5 + 
           land1 + land2 + land3 + land4 + land5 + land6  +
           region2 + region3 + region4 + region5 + region6 +
           urban1 + urban2 + urban3 +
           # fixed effects #
           accessibility + aridity + cropland +  density + elevation + 
             evi + footprint + gdp + growing + hdi + irrigation + 
             lst_range + night_light + pasture + pet + water_dist +
           # spatial effect via SPDE #
           spat(main=coo, model = matern_spde),
         # data and likelihood #
          data=D, family="binomial")
summary(fit_sp)

saveRDS(fit_sp, paste0(path_hvm, 'inla_model_cat_',  class_name,  '.rds'), compress = TRUE)


# fit_sp <- readRDS(paste0(path_hvm, 'inla_model_cat_pigs_Pigs_Yes.rds'))

# # correlation between sites
# spde.posterior(fit_sp, "spat", 
#                what = "matern.correlation") %>%
# plot() + labs(x="Distance between sites (km)",
#               y="Correlation of spatial effect") +
#           theme_bw(base_size = 13)
# 
# # exmaine the (latent) spatial field 
# gproj <- inla.mesh.projector(mesh,dims = c(1200, 1200)) #was 1000 x 1000
# sp.mean <- inla.mesh.project(gproj, 
#                              fit_sp$summary.random$spat$mean)
# sp.sd <- inla.mesh.project(gproj, 
#                            fit_sp$summary.random$spat$sd)
# library(gridExtra)
# library(lattice)
# gridExtra::grid.arrange(levelplot(sp.mean, scales=list(draw=F), xlab='', ylab='', 
#                                   main='Spatial effect mean',col.regions = heat.colors(16)),
#                         levelplot(sp.sd, scal=list(draw=F), xla='', yla='', 
#                                   main='Spatial effect sd' ,col.regions = heat.colors(16)), 
#                         nrow=2)


############################################################################################################
# Generate predictions
# for simplicity just renamed the floors data frame

df_raster <- read.csv(file = paste0(path_hvm, 'df_raster_all_nonnan_cat.csv'),
                      row.names = 1, na.strings = "")
  
df_raster$urban <- as.factor(df_raster$urban)
df_raster$land <- as.factor(df_raster$land)
df_raster$climate <- as.factor(df_raster$climate)
df_raster$region <- as.factor(df_raster$region)
df_raster_time <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(df_raster))), c('time'))
df_raster_time[is.na(df_raster_time)] = 8.7908778
df_raster <- data.frame(cbind(df_raster_time, df_raster))

d_test<-df_raster
# d_test <- subset(d_test, select = -c(wwtp_dist, poultry2, ruminants2))
# d_test <- subset(d_test, select = -c(wwtp_dist, pigs2, poultry2, ruminants2))
# N_test<-nrow(d_test)
# d_test$Y_bin<-as.numeric(d_test$target_var == 'Flr_Fin')

# make dummy vars for factors  (required by inlabru)
model.matrix(~climate, d_test) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> clim
model.matrix(~land, d_test) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> land
model.matrix(~region, d_test) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> reg
model.matrix(~urban, d_test) %>% # create dummy variable for month
  as.data.frame() %>%
  dplyr::select(-1) %>% # drop intercept
  bind_cols() -> urb
## append to original (use D from now on)
D_test<-cbind(d_test,clim,land,reg,urb)
# D_test <- D_test[1:10000,]

#** set coordinate system and reproject into meters **#
coordinates(D_test)=~Longitude+Latitude
#initial coordinate system
proj4string(D_test)= CRS("+proj=longlat +ellps=sphere")  
#project into Mollweide with KM as units
D_proj_test <- spTransform(D_test, CRS("+proj=moll +units=km")) %>% 
  st_as_sf() 

# record prediction data coordinates to evaluate spatial effect
pred_coord_X<-st_coordinates(D_proj_test)[,1]
pred_coord_Y<-st_coordinates(D_proj_test)[,2]



# Split the data to be predicted into 10 parts
row_number <- nrow(D_proj_test)
quotient <- row_number %/% 10
mod_number <- row_number %/% 10 * 10
remainder <- row_number %% 10
df_seq <- seq(0, mod_number, by = quotient)
df_seq[length(df_seq)] <- row_number # Length of the last interval 

for (x in 1:10)
  {D_proj_test_sub <- D_proj_test[(df_seq[x]+1):df_seq[x+1],]
   y_pred <- predict(object = fit_sp, 
                     newdata = D_proj_test_sub, formula = ~ Intercept + time +
                      # dummy vars for all categorical variables #
                      climate2 + climate3 + climate4 + climate5 +
                      land1 + land2 + land3 + land4 + land5 + land6 +
                      region2 + region3 + region4 + region5 + region6 +
                      urban1 + urban2 + urban3 +
                      # fixed effects #
                      accessibility + aridity + cropland +  density + elevation + 
                      evi + footprint + gdp + growing + hdi + irrigation + 
                      lst_range + night_light +pasture + pet + water_dist +
                      # Note that we evaluate the spatial field at prediction coords
                      spat_eval(cbind(pred_coord_X[(df_seq[x]+1):df_seq[x+1]],
                                      pred_coord_Y[(df_seq[x]+1):df_seq[x+1]]))
                    ) %>%
                    mutate(pred_prob = plogis(q0.5),
                           pred_lower = plogis(q0.025),
                           pred_upper = plogis(q0.975))
   write.csv(y_pred, paste0(path_pred, 'y_pred_cat_part_', 
                            sprintf("%02d", x), '.csv'), row.names = FALSE)
  rm(D_proj_test_sub, y_pred)
  print(x)}

rm(df_raster, D_proj_test)



# head(y_pred)


# y_pred <- predict(fit_sp, NULL, newdata = D_test, type='response')
# y_pred_prob <- plogis(y_pred$Predictor$median)
# pred_coord <- coordinates(D_test[1:100,])
















# ## take random 20% of data
# #use 80% of dataset as training set and 20% as holdout
# samp <- sample(c(TRUE, FALSE), nrow(D),
#                replace=TRUE, prob=c(0.8,0.2))
# # data
# train <- D[samp, ]
# test <- D[!samp, ]
# # coordinates
# coo_train <- coo[samp,]
# coo_test <- coo[!samp,]
# ## mesh, priors, estimate using train
# # build mesh
# mesh_train <- inla.mesh.2d(loc = coo_train,
#                            max.edge = c(max_edge, 10*max_edge), # buffer to avoid edge effects
#                            cutoff = 0.20*max_edge)
# # priors #
# matern_spde_train <-
#   inla.spde2.pcmatern(mesh_train,
#                       # set priors for spatial SD and range
#                       prior.sigma = c(3, 0.10), #P(SD>3) = 0.1
#                       prior.range = c(100, 0.10) #P(range<100km) = 0.1
#   )
# # estimate #
# ###
# fit_train<-bru(Y_bin ~ Intercept(1) +
#                  # non-linear smooth time trend #
#                  trend(time, model="rw1") +
#                  # dummy vars for all categorical variables #
#                  climate2 + climate3 + climate4 + climate5 +
#                  region2 + region3 + region4 + region5 + region6 +
#                  urban1 + urban2 + urban3 +
#                  # linear fixed effects #
#                  accessibility + aridity + cropland +
#                  density + elevation + evi + footprint + growing +
#                  irrigation + lst_range + nighttime_light +
#                  pasture + pet + water_dist +
#                  # spatial effect via SPDE #
#                  spat(main=coo_train, model = matern_spde_train),
#                # data and likelihood #
#                data=train,family="binomial")
# ## predict the hold-out
# pxl<-fm_pixels(mesh_train,dims = c(50,50))
# test_pred<-predict(fit_train, pxl, newdata=test, type="response")
# head(test_pred)

###################################################################################
# Generate feature importance

fit_sp <- readRDS(paste0(path_hvm, 'inla_model_cat_poultry_Ptry_Yes.rds'))

fit_sp_table <- fit_sp$summary.fixed
fit_sp_table$value <- abs(fit_sp_table$mean / fit_sp_table$sd)
fit_sp_table <- fit_sp_table['value']
fit_sp_table <- fit_sp_table[!rownames(fit_sp_table) %in% 'Intercept',, drop = FALSE]
# sorted_row <- fit_sp_table[order(-fit_sp_table$value), ]
# fit_sp_table['value'] <- sorted_row


write.csv(fit_sp_table, 
          paste0(path_fi, 'results_250108/feature_importance_new/df_feature_imp_ruminant.csv'), 
          row.names = TRUE)




# Prepare data for SHAP estimation
# Extract fixed effect means
fixed_effects <- fit_sp$summary.fixed$mean
fixed_effects_df <- data.frame(Variable = fit_sp$names.fixed, Coefficient = fixed_effects)

# Extract predictions
predictions <- fit_sp$summary.fitted.values$mean
predictions_df <- data.frame(Predictions = predictions)

dataset <- D@data

# Save to CSV
setwd('/Users/binfang/Documents/Household_variable_mapping/results/results_250108/SHAP')
write.csv(fixed_effects_df, "fixed_effects_poultry.csv", row.names = FALSE)
write.csv(predictions_df, "predictions_poultry.csv", row.names = FALSE)
write.csv(dataset, "dataset_poultry.csv", row.names = FALSE)
