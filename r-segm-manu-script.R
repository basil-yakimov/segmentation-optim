
#### ---- libraries ----

require(raster)
require(cleangeo)
require(dplyr)
require(lidR)
require(ForestTools)
require(DEoptim)

#### ---- data ----

crowns <- raster::shapefile('data/crowns.shp')
rawCHM <- raster::raster('data/CHM.tif')
NDVI <- raster::shapefile('data/NDVI-vector.shp')

#### ---- pre-processing ----

IoU <- function(crowns, segments, segm_method){
  
  crowns@data$crID <- 1:nrow(crowns@data)
  crowns@data$crA <- raster::area(crowns)
  
  segments.cl <- cleangeo::clgeo_CleanByPolygonation.SpatialPolygons(segments, verbose = F)
  segments.cldf <- sp::SpatialPolygonsDataFrame(segments.cl, data = segments@data, match.ID = F)
  segments <- segments.cldf
  
  segments@data$segID <- 1:nrow(segments@data)
  segments@data$segA <- raster::area(segments)
  
  inter <- raster::intersect(crowns, segments)
  
  df <- inter@data |> dplyr::select(crID, segID, crA, segA) |> 
    dplyr::mutate(id = paste(crID, segID, sep = "/"),
                  A = raster::area(inter)) |>
    dplyr::group_by(id) |>
    dplyr::summarise(dplyr::across(crID:segA, dplyr::first),
                     A = sum(A)) |>
    dplyr::mutate(IoU = A/(crA + segA - A)*100)
  
  TP <- sum(df$IoU > 50)
  
  METRICS <- data.frame(segm_method = segm_method,
                        TP = TP,
                        FP = length(segments) - TP,
                        FN = length(crowns) - TP)
  
  METRICS$precision <- (METRICS$TP / (METRICS$TP + METRICS$FP))
  METRICS$recall <- (METRICS$TP / (METRICS$TP + METRICS$FN))
  METRICS$F_score <- ((2 * METRICS$precision * METRICS$recall) / (METRICS$precision + METRICS$recall))
  
  if(is.na(METRICS$F_score)){METRICS$F_score <- 0}
  
  return(METRICS)
}

# ---- –ALGORITHMS FUNC.– ----

save_dir <- 'D:\\NIKITA\\SCience\\Article\\public-script\\iter\\'
save_dir_gens <- paste0(save_dir, 'GENs\\'); dir.create(path = save_dir_gens)


## ----- Watershed -----

fn_ws <- function(par){
  
  name <- paste0(round(par[1], digits = 5), "#", round(par[2], digits = 5), "#", round(par[3], digits = 5), "#", round(par[4], digits = 5))
  AP <- data.frame(sigma = par[1], th_tree = par[2], tol = par[3], ext = par[4])
  
  start_time <- Sys.time()
  
  date <- as.character(start_time)
  date <- strsplit(date, " ")
  date <- date[[1]][1]
  
  par[4] <- floor(par[4])
  
  if(par[1] == 0){CHM <- rawCHM} else 
  {
    CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  }
  
  smoothedCHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  raster::inMemory(CHM)
  CHM <- raster::readAll(CHM)
  
  polygons <- lidR::watershed(CHM, th_tree = par[2], tol = par[3], ext = par[4])()
  
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "Watershed")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  OUT <- 1 - stat$F_score
  
  nameslist <- list.files(path = save_dir, pattern = "rda", all.files = T, full.names = F); iteration <- length(nameslist) + 1
  INFO <- list(param = AP, stat = stat, time = elapsed_time, F_score = stat$F_score)
  save(INFO, file = paste0(save_dir, "ws-", iteration, '-', date, '#', name, '_', round(stat$F_score, digits = 5), '.rda'))
  
  print(INFO)
  
  return(OUT)
}

par_def <- c(0.67, 0, 1, 1)
1 - fn_ws(par_def)

## ----- MCWS -----

# par <- c(4.9317, 0.1202, 1.7371, 7.42, 5.0904)

fn_mcws <- function(par){
  
  name <- paste0(round(par[1], digits = 5), "#", round(par[2], digits = 5), "#", round(par[3], digits = 5), "#", round(par[4], digits = 5))
  AP <- data.frame(sigma = par[1], th_tree = par[2], tol = par[3], ext = par[4])
  
  start_time <- Sys.time()
  
  date <- as.character(start_time)
  date <- strsplit(date, " ")
  date <- date[[1]][1]
  
  if(par[1] == 0){CHM <- rawCHM} else 
  {
    CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  }
  
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  raster::inMemory(CHM)
  CHM <- raster::readAll(CHM)
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons <- ForestTools::mcws(treetops = treetops, CHM = CHM, format = "polygons", minHeight = par[5])
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "MCWS")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  OUT <- 1 - stat$F_score
  
  nameslist <- list.files(path = save_dir, pattern = "rda", all.files = T, full.names = F); iteration <- length(nameslist) + 1
  INFO <- list(param = AP, stat = stat, time = elapsed_time, F_score = stat$F_score)
  save(INFO, file = paste0(save_dir, "mcws-", iteration, '-', date, '#', name, '_', round(stat$F_score, digits = 5), '.rda'))
  
  print(INFO)
  
  return(OUT)
}

par_def <- c(0.67, 0.2414, 3.192, 2.51, 0)
1 - fn_mcws(par_def)

## ----- Dalponte -----

# par <- c(9.9353, 0.1099, 0.9641, 7.458, 5.0784, 0.0573, 0.5532)

fn_d <- function(par){
  
  name <- paste0(round(par[1], digits = 5), "#", round(par[2], digits = 5), "#", round(par[3], digits = 5), "#", round(par[4], digits = 5))
  AP <- data.frame(sigma = par[1], th_tree = par[2], tol = par[3], ext = par[4])
  
  start_time <- Sys.time()
  
  date <- as.character(start_time)
  date <- strsplit(date, " ")
  date <- date[[1]][1]
  
  if(par[1] == 0){CHM <- rawCHM} else 
  {
    CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  }
  
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  raster::inMemory(CHM)
  CHM <- raster::readAll(CHM)
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons  <- lidR::dalponte2016(CHM, treetops, th_tree = par[5], th_seed = par[6], 
                                      th_cr = par[7], max_cr = 308)()
  
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "Dalponte")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  OUT <- 1 - stat$F_score
  
  nameslist <- list.files(path = save_dir, pattern = "rda", all.files = T, full.names = F); iteration <- length(nameslist) + 1
  INFO <- list(param = AP, stat = stat, time = elapsed_time, F_score = stat$F_score)
  save(INFO, file = paste0(save_dir, "d-", iteration, '-', date, '#', name, '_', round(stat$F_score, digits = 5), '.rda'))
  
  print(INFO)
  
  return(OUT)
}

par <- c(0.67, 0.2414, 3.192, 2.51, 0, 0.45, 0.55)
1 - fn_d(par)

## ----- Silva -----

# par <- c(10.128, 0.0357, 2.1206, 7.6843, 0.4549, 0.4279)

fn_s <- function(par){
  
  name <- paste0(round(par[1], digits = 5), "#", round(par[2], digits = 5), "#", round(par[3], digits = 5), "#", round(par[4], digits = 5))
  AP <- data.frame(sigma = par[1], th_tree = par[2], tol = par[3], ext = par[4])
  
  start_time <- Sys.time()
  
  date <- as.character(start_time)
  date <- strsplit(date, " ")
  date <- date[[1]][1]
  
  if(par[1] == 0){CHM <- rawCHM} else 
  {
    CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  }
  
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  raster::inMemory(CHM)
  CHM <- raster::readAll(CHM)
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons <- lidR::silva2016(CHM, treetops, max_cr_factor = par[5], exclusion = par[6])()
  
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "Silva")
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  OUT <- 1 - stat$F_score
  
  nameslist <- list.files(path = save_dir, pattern = "rda", all.files = T, full.names = F); iteration <- length(nameslist) + 1
  INFO <- list(param = AP, stat = stat, time = elapsed_time, F_score = stat$F_score)
  save(INFO, file = paste0(save_dir, "s-", iteration, '-', date, '#', name, '_', round(stat$F_score, digits = 5), '.rda'))
  
  print(INFO)
  
  return(OUT)
}

par <- c(0.67, 0.2414, 3.192, 2.51, 0.6, 0.3)
1 - fn_s(par)

# ---- –DE OPTIMIZATION– ----

save_dir <- 'D:\\NIKITA\\SCience\\Article\\public-script\\iter\\de-';

# Watershed

par_lower <-  c(0,  0,  0,  1)
par_upper <-  c(15, 15, 5, 50)

numCores <- detectCores() - 6

cl <- makeCluster(numCores)

DE.WS <-  DEoptim(fn_ws, par_lower, par_upper,
                  control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                            cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                            parVar = c('crowns', 'rawCHM', 'NDVI', 'IoU', 'save_dir')))

save(DE.WS, file = paste0(save_dir_gens, 'DE.WS_', 0,".rda"))

for (ii in 1:29){
  
  numCores <- detectCores() - 6
  
  cl <- makeCluster(numCores)
  
  DE.WS <-  DEoptim(fn_ws, par_lower, par_upper,
                    control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                              cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                              parVar = c('crowns', 'rawCHM','NDVI', 'IoU', 'save_dir'),
                                              storepopfrom = 0, initialpop = DE.WS$member$pop))
  
  save(DE.WS, file = paste0(save_dir_gens, 'DE.WS_', ii,".rda"))}



# MCWS

par_lower <- c(0,  0,    0.001, 0,  0)
par_upper <- c(15, 0.3,  5,     15, 15)

numCores <- detectCores() - 6

cl <- makeCluster(numCores)

DE.MCWS <-  DEoptim(fn_mcws, par_lower, par_upper,
                  control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                            cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                            parVar = c('crowns', 'rawCHM', 'NDVI', 'IoU', 'save_dir')))

save(DE.MCWS, file = paste0(save_dir_gens, 'DE.MCWS_', 0,".rda"))

for (ii in 1:29){
  
  numCores <- detectCores() - 6
  
  cl <- makeCluster(numCores)
  
  DE.MCWS <-  DEoptim(fn_mcws, par_lower, par_upper,
                    control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                              cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                              parVar = c('crowns', 'rawCHM','NDVI', 'IoU', 'save_dir'),
                                              storepopfrom = 0, initialpop = DE.MCWS$member$pop))
  
  save(DE.MCWS, file = paste0(save_dir_gens, 'DE.MCWS_', ii,".rda"))}




# Dalponte

par_lower <- c(0,  0,   0.001, 0,  0,  0, 0)
par_upper <- c(15, 0.3, 5,     15, 15, 1, 1)

numCores <- detectCores() - 6

cl <- makeCluster(numCores)

DE.D <-  DEoptim(fn_d, par_lower, par_upper,
                    control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                              cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                              parVar = c('crowns', 'rawCHM', 'NDVI', 'IoU', 'save_dir')))

save(DE.D, file = paste0(save_dir_gens, 'DE.D_', 0,".rda"))

for (ii in 1:29){
  
  numCores <- detectCores() - 6
  
  cl <- makeCluster(numCores)
  
  DE.D <-  DEoptim(fn_d, par_lower, par_upper,
                      control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                                cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                                parVar = c('crowns', 'rawCHM','NDVI', 'IoU', 'save_dir'),
                                                storepopfrom = 0, initialpop = DE.D$member$pop))
  
  save(DE.D, file = paste0(save_dir_gens, 'DE.D_', ii,".rda"))}


# Silva

lower <- c(0,  0,   0.001,  0,  0.15, 0)
upper <- c(15, 0.3, 5,      15,    2, 1)

numCores <- detectCores() - 6

cl <- makeCluster(numCores)

DE.S <-  DEoptim(fn_s, par_lower, par_upper,
                 control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                           cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                           parVar = c('crowns', 'rawCHM', 'NDVI', 'IoU', 'save_dir')))

save(DE.S, file = paste0(save_dir_gens, 'DE.S_', 0,".rda"))

for (ii in 1:29){
  
  numCores <- detectCores() - 6
  
  cl <- makeCluster(numCores)
  
  DE.S <-  DEoptim(fn_s, par_lower, par_upper,
                   control = DEoptim.control(strategy = 2, NP = 20, itermax = 1, parallelType = 'parallel',
                                             cluster = cl, packages = c("raster", "sp", "cleangeo", "spatialEco", "lidR", "ForestTools", "tidyverse", "dplyr"), 
                                             parVar = c('crowns', 'rawCHM','NDVI', 'IoU', 'save_dir'),
                                             storepopfrom = 0, initialpop = DE.S$member$pop))
  
  save(DE.S, file = paste0(save_dir_gens, 'DE.S_', ii,".rda"))}


# ---- –RS OPTIMIZATION– ----

save_dir <- 'D:\\NIKITA\\SCience\\Article\\public-script\\iter\\rs-';

# Watershed

par_lower <-  c(0,  0,  0,  1)
par_upper <-  c(15, 15, 5, 50)

#


for(i in 1:iter_max){
  
  rs_sigma <- runif(1, par_lower[1], par_upper[1])
  rs_th_tree <- runif(1, par_lower[2], par_upper[2])
  rs_tol <- runif(1, par_lower[3], par_upper[3])
  rs_ext <- runif(1, par_lower[4], par_upper[4])
  
  par_list[[i]] <- c(rs_sigma, rs_th_tree, rs_tol, rs_ext)
}

##

cl <- makeCluster(10) #cl <- makeCluster(3)
registerDoParallel(cl)

foreach(par = par_list) %dopar% fn_ws(par)

stopCluster(cl)


# MCWS

par_lower <- c(0,  0,    0.001, 0,  0)
par_upper <- c(15, 0.3,  5,     15, 15)

for(i in 1:iter_max){
  
  rs_sigma <- runif(1, par_lower[1], par_upper[1])
  rs_reg_slope <- runif(1, par_lower[2], par_upper[2])
  rs_reg_int <- runif(1, par_lower[3], par_upper[3])
  rs_hmin_tt <- runif(1, par_lower[4], par_upper[4])
  rs_hmin_cr <- runif(1, par_lower[5], par_upper[5])
  
  par_list[[i]] <- c(rs_sigma, rs_reg_slope, rs_reg_int, rs_hmin_tt, rs_hmin_cr)
}

##

cl <- makeCluster(10) #cl <- makeCluster(3)
registerDoParallel(cl)

foreach(par = par_list) %dopar% fn_mcws(par)

stopCluster(cl)


# Dalponte

par_lower <- c(0,  0,   0.001, 0,  0,  0, 0)
par_upper <- c(15, 0.3, 5,     15, 15, 1, 1)

# sigma, reg.slope, reg.int, hmin, th_tree, th_seed, th_cr

for(i in 1:iter_max){
  
  rs_sigma <- runif(1, par_lower[1], par_upper[1])
  rs_reg_slope <- runif(1, par_lower[2], par_upper[2])
  rs_reg_int <- runif(1, par_lower[3], par_upper[3])
  rs_hmin <- runif(1, par_lower[4], par_upper[4])
  rs_th_tree <- runif(1, par_lower[5], par_upper[5])
  rs_th_seed <- runif(1, par_lower[6], par_upper[6])
  rs_th_cr <- runif(1, par_lower[7], par_upper[7])
  
  par_list[[i]] <- c(rs_sigma, rs_reg_slope, rs_reg_int, rs_hmin, rs_th_tree, rs_th_seed, rs_th_cr)
}

##

cl <- makeCluster(10) #cl <- makeCluster(3)
registerDoParallel(cl)

foreach(par = par_list) %dopar% fn_d(par)

stopCluster(cl)


# Silva

par_lower <- c(0,  0,   0.001,  0,  0.15, 0)
par_upper <- c(15, 0.3, 5,      15,    2, 1)

# AP <- data.frame(sigma = par[1], reg.slope = par[2], reg.int = par[3], hmin = par[4], max_cr_factor = par[5], exclusion = par[6])

for(i in 1:iter_max){
  
  rs_sigma <- runif(1, par_lower[1], par_upper[1])
  rs_reg_slope <- runif(1, par_lower[2], par_upper[2])
  rs_reg_int <- runif(1, par_lower[3], par_upper[3])
  rs_hmin <- runif(1, par_lower[4], par_upper[4])
  rs_max_cr_factor <- runif(1, par_lower[5], par_upper[5])
  rs_exclusion <- runif(1, par_lower[6], par_upper[6])
  
  par_list[[i]] <- c(rs_sigma, rs_reg_slope, rs_reg_int, rs_hmin,  rs_max_cr_factor, rs_exclusion)
}

##

cl <- makeCluster(10) #cl <- makeCluster(3)
registerDoParallel(cl)

foreach(par = par_list) %dopar% fn_s(par)

stopCluster(cl)



# ---- –ABLATION– ----


## --- ablation function ---

ablation <- function(optimal_params, objective_function, baseline_values = NULL) {
  # optimal_params: найденный оптимальный вектор параметров
  # objective_function: ваша функция accuracy/целевая функция
  # baseline_values: базовые значения для "обнуления" параметров
  
  n_params <- length(optimal_params)
  
  # # Если базовые значения не заданы, используем средние по диапазону или нули
  # if (is.null(baseline_values)) {
  #   baseline_values <- rep(0, n_params)  # или другие разумные значения
  # }
  
  # Вычисляем базовое значение F(θ) с оптимальными параметрами
  base_performance <- objective_function(optimal_params)
  
  # Вектор для хранения дельт
  delta_F <- numeric(n_params)
  
  # Абляция для каждого параметра
  for (i in 1:n_params) {
    # Создаем модифицированный вектор параметров
    modified_params <- optimal_params
    modified_params[i] <- baseline_values[i]  # "обнуляем" i-й параметр
    
    # Вычисляем производительность без i-го параметра
    performance_without_i <- objective_function(modified_params)
    
    # Вычисляем дельту
    delta_F[i] <- base_performance - performance_without_i
  }
  
  # Создаем dataframe с результатами
  results <- data.frame(
    parameter = 1:n_params,
    baseline_value = baseline_values,
    optimal_value = optimal_params,
    base_performance = base_performance,
    performance_without_param = base_performance - delta_F,
    delta_F = delta_F,
    abs_delta_F = abs(delta_F)
  )
  
  return(results)
}

# ---- Watershed ----

par_optim.ws <- c(9.6069, 5.2519, 0.02, 11)
par_base.ws <- c(0.67, 0, 1, 1)

segm_ws <- function(par){
  
  CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  polygons <- lidR::watershed(CHM, th_tree = par[2], tol = par[3], ext = par[4])()
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "watershed")
  
  return(stat$F_score)
}


# ---- MCWS ----

par_optim.mcws <- c(4.9317, 0.1202, 1.7371, 7.42, 5.0904)
par_base.mcws <- c(0.67, 0.2414, 3.192, 2.51, 0)

segm_mcws <- function(par){
  
  CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons <- ForestTools::mcws(treetops = treetops, CHM = CHM, format = "polygons",
                                minHeight = par[5])
  
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "MCWS")
  
  return(stat$F_score)
}


# ---- Dalponte ----

par_optim.d <- c(9.9353, 0.1099, 0.9641, 7.458, 5.0784, 0.0573, 0.5532)
par_base.d <- c(0.67, 0.2414, 3.192, 2.51, 0, 0.45, 0.55)

segm_d <- function(par){
  
  CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons <- lidR::dalponte2016(CHM, treetops, th_tree = par[5], th_seed = par[6], 
                                 th_cr = par[7], max_cr = 308)()
  
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "Dalponte")
  
  return(stat$F_score)
}


# ---- Silva ----

par_optim.s <- c(10.128, 0.0357, 2.1206, 7.6843, 0.4549, 0.4279)
par_base.s <- c(0.67, 0.2414, 3.192, 2.51, 0.6, 0.3)

segm_s <- function(par){
  
  CHM <- spatialEco::raster.gaussian.smooth(terra::rast(rawCHM), n = 21, s = par[1]) |> raster::raster() # сглаживание 
  CHM <- raster::mask(CHM, NDVI)
  CHM <- CHM - CHM@data@min
  
  lin <- function(x) {(x)*par[2] + par[3]}
  treetops <- lidR::locate_trees(CHM, lidR::lmf(ws = lin, hmin = par[4]))
  
  polygons <- lidR::silva2016(CHM, treetops, max_cr_factor = par[5], exclusion = par[6])()
  
  polygons <- terra::rast(polygons)
  polygons <- terra::as.polygons(polygons, extent = F, dissolve = T)
  polygons <- as(polygons, "Spatial")
  
  stat <- IoU(crowns = crowns, segments = polygons, segm_method = "Dalponte")
  
  return(stat$F_score)
}

# ---- study and visualization ----

ab_res_ws <- ablation(par_optim.ws, segm_ws, par_base.ws)
ab_res_mcws <- ablation(par_optim.mcws, segm_mcws, par_base.mcws)
ab_res_d <- ablation(par_optim.d, segm_d, par_base.d)
ab_res_s <- ablation(par_optim.s, segm_s, par_base.s)


ab_res_ws$parameter_name <- c("sigma", "th_tree", "tol", "ext")
ab_res_mcws$parameter_name <- c("sigma", "reg_slope", "reg_int", "hmin.tt", "hmin.cr")
ab_res_d$parameter_name <- c("sigma", "reg_slope", "reg_int", "hmin", "th_tree", "th_seed", "th_cr")
ab_res_s$parameter_name <- c("sigma", "reg_slope", "reg_int", "hmin", "max_cr_factor", "exclusion")

ab_res_ws$method <- "Watershed" ; ab_res_ws$gp <- 1
ab_res_mcws$method <- "MCWS" ; ab_res_mcws$gp <- 2
ab_res_d$method <- "Dalponte" ; ab_res_d$gp <- 3
ab_res_s$method <- "Silva" ; ab_res_s$gp <- 4

list <- list(ab_res_ws, ab_res_mcws, ab_res_d, ab_res_s)

df <- bind_rows(list, .id = "gp")


plot <- ggplot(data = df, aes(x = reorder(parameter_name, delta_F), y = delta_F, fill = delta_F))+
  theme_bw(base_size = 15)+
  geom_bar(stat = "identity") +
  scale_fill_viridis_c() +
  labs(x = NULL, y = "ΔF", fill = "Impact")+
  facet_wrap(~ method, nrow = 2, scales = 'free')+
  coord_flip()
# annotate(geom = "text", x = 30, y = max(dfi_w$F_score) + 0.007, fontface = "italic",
#          label = round(max(dfi_w$F_score), digits = 5))+
# xlim(0, 31)

# plot

ggsave("ablation_graph-4.png", plot, path = "###",
       width = 18.46, height = 13, units = "cm",
       scale = 1.1) 


