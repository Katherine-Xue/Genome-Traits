
library(ggplot2)

setwd("0_MetaDiverse/MapPrediction")

load("0_MetaDiverse/GenomeTraits.RData")

colnames(df1)
df2 <- df1[,c(1,11,12,31,35,34,32,33)]
### load rasters ----

library(SLGACloud)

AllLayers

names(AllLayers) <- c("Clay","Sand","TP","pH","TN","SOC","CEC",
                           "MAP","ADM","MAT",
                           "TRA", "PVA")
###
library(dplyr)
colnames(df2)

df3 <- df2 %>%
  group_by(X_Longtitude, Y_Latitude) %>%
  summarise(
    across(c(GC, Genome_Size, geneDensity,cells,`Cells(log10)`), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

###
lonlat <- df3[,c("X_Longtitude", "Y_Latitude")]

names(AllLayers)

extr1 <- raster::extract(AllLayers,lonlat)

df4 <- cbind(df3,extr1[,-1])

save(df4, file = "Cov.train.RData")

####
library(randomForest)
library(missForest)
library(caret)

colnames(df4)

rf_model <- randomForest(
  x = df4[,8:19],     # predictors
  y = df4[,7],        # response (column 1)
  importance = TRUE,
  mtry=5,
  nodesize=5,
  ntree = 500)

rf_model

varImpPlot(rf_model)

#### prediction
colnames(df4)
colnames(df4)[7] <- "cell_log"

covarable <- colnames(df4)[8:19]
covarable

Imp <- data.frame()

tab.res <- data.frame(matrix(data=NA, nrow = 5, ncol = 5))
rownames(tab.res) <- colnames(df4)[3:7]
colnames(tab.res) <- c('ME', 'RMSE', 'r2', 'R2', 'rhoC')

for (i in 3:7){
  set.seed(123)
  pre.name <- paste0("res_",colnames(df4)[i])
  
  varname <- colnames(df4)[i]
  varname
  pre.name <- make_cv_RF(data = df4[,c(colnames(df4)[i],covarable)],
                         property = colnames(df4)[i],
                         covars =  covarable,
                         raster = ras.all, 
                         spatial.cv = F, 
                         nfolds = 10,
                         tunning = F)
  
  pred <- pre.name[[2]]$pred
  obs <- pre.name[[2]][,1]
  
  res <- eval_all(obs, pred)
  
  tab.res[varname, ] <- as.numeric(res[1, ])   # ✅ 用 k，不是 i
 
  plot(pre.name[[2]][,1], pre.name[[2]]$pred); abline(0,1,col="red")
  
  pre.name[[1]]$variable.importance[order(pre.name[[1]]$variable.importance)]
  pre.Imp <- data.frame("IncMSE"=pre.name[[1]]$variable.importance)
  pre.Imp$Imp <- round(pre.Imp[,1]/sum(pre.Imp[,1])*100,2)
  colSums(pre.Imp)
  
  pre.Imp$Var <- colnames(df4)[i]
  pre.Imp$factor <- row.names(pre.Imp)
  
  Imp <- rbind(Imp,pre.Imp)
  
  res.model <- pre.name
  res.model[[1]]
  path1 <- paste0('./','Model_',colnames(df4)[i],'_prediction.Rdata')
  save(res.model, file =path1)
}

write.csv(data.frame(tab.res),file = "tab.res.csv")
write.csv(Imp,file = "Imp.csv")

#### map prediction
library(raster)
library(terra)
library(ranger)

### resample 
r <- AllLayers[[1]]
res(AllLayers)
res(r) <- 0.01
ras <- resample(AllLayers, r)
names(ras)

ras1 <- stat(ras)
ras1
names(ras1)
class(ras1)

plot(ras1)

writeRaster(
  ras1,
  file.path(getwd(), "TERN_covariates_stack_0.01.tif"),
  overwrite = TRUE,
  wopt = list(
    gdal = c(
      "COMPRESS=LZW",
      "TILED=YES",
      "BIGTIFF=YES"
    )
  )
)

### bloklize

ras1 <- raster::stack("MapPrediction/TERN_covariates_stack_0.01.tif")
covariates.stack <- raster::stack(ras1)

bs <- raster::blockSize(covariates.stack,chunksize = 10000000, minblocks = 50)
bs

###

components <- colnames(df4)[6:7]
components

d=1

for (d in 1:length(components)){ # for all depths
  
  # which depth
  compNumb <- components[d]
  compNumb
  # load model
  x <- load(paste0('./','Model_',compNumb,'_prediction.Rdata'))
    model <- get(x)[[1]] ; rm(x)
  model
  FUN <- function(i){
    #tile_temp <- readAll(crop(covariates.stack, tiles.pol[i]))
    tile_temp <- readAll(crop(covariates.stack, extent(covariates.stack, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, ncol(covariates.stack))))
    
    #tile_temp2 <- as.data.frame(tile_temp, na.)
    #predict(object = model, data = tile_temp2)
    pred <- raster::predict(object = tile_temp, 
                            model = model, 
                            fun = function(model, ...) predict(model, ...)$predictions)
    #tile_temp2 <- as.data.frame(tile_temp2, )
    predvar05 <- raster::predict(object = tile_temp, 
                                 model = model, 
                                 fun = function(model, ...) predict(model,type = 'quantiles', quantiles =c(0.05), ...)$predictions[,1])
    predvar95 <- raster::predict(object = tile_temp, 
                                 model = model, 
                                 fun = function(model, ...) predict(model,type = 'quantiles', quantiles =c(0.95), ...)$predictions[,1])
    
    # link to path where the files are saves
    saveRDS(pred, file = paste0('./pred', compNumb, '/pred/', i, '.rds'))
    saveRDS(predvar05, file = paste0('./pred', compNumb, '/predvar05/', i, '.rds'))
    saveRDS(predvar95, file = paste0('./pred', compNumb, '/predvar95/', i, '.rds'))
    
  }
  
  # create files in the path
  dir.create(paste0('./pred', compNumb))
  
  dir.create(paste0('./pred', compNumb, '/pred'))
  dir.create(paste0('./pred', compNumb, '/predvar05'))
  dir.create(paste0('./pred', compNumb, '/predvar95'))
  
  ## not in parallel 
  library(plyr)
  llply(seq_along(1:bs$n), FUN)
  
  #####
  ##### MERGE PREDICTED MAPS
  # list all files for predicted maps
  rasters1 <- list.files(paste0('./pred', compNumb, '/pred'))
  
  # read all rasters and merge
  rast.list <- list()
  for(i in 1:length(rasters1)) { 
    tile_temp <- readRDS(paste0('./pred', compNumb, '/pred/', rasters1[i]))
    rast.list[i] <- tile_temp}
  
  # And then use do.call on the list of raster objects
  names(rast.list)[1:2] <- c('x', 'y')
  rast.list$fun <- mean
  rast.list$na.rm <- TRUE
  y_pred <- do.call(mosaic, rast.list)
  
  plot(y_pred)
  
  writeRaster(y_pred, filename =paste0("./MDS_", compNumb, "_pred.tif"), overwrite = TRUE)
  
  #####
  ##### MERGE PREDICTED MAPS
  
  rasters1 <- list.files(paste0('./pred', compNumb, '/predvar05'))
  
  rast.list <- list()
  for(i in 1:length(rasters1)) { 
    tile_temp <- readRDS(paste0('./pred', compNumb, '/predvar05/', rasters1[i]))
    rast.list[i] <- tile_temp}
  
  # And then use do.call on the list of raster objects
  names(rast.list)[1:2] <- c('x', 'y')
  rast.list$fun <- mean
  rast.list$na.rm <- TRUE
  y_predvar <- do.call(mosaic, rast.list)
  
  plot(y_predvar)
  writeRaster(y_predvar, filename =paste0("./MDS_", compNumb, "_predvar05p.tif"), overwrite = TRUE)
  
  #####
  ##### MERGE PREDICTED MAPS
  
  rasters1 <- list.files(paste0('./pred', compNumb, '/predvar95'))
  
  rast.list <- list()
  for(i in 1:length(rasters1)) { 
    tile_temp <- readRDS(paste0('./pred', compNumb, '/predvar95/', rasters1[i]))
    rast.list[i] <- tile_temp}
  
  # And then use do.call on the list of raster objects
  names(rast.list)[1:2] <- c('x', 'y')
  rast.list$fun <- mean
  rast.list$na.rm <- TRUE
  y_predvar <- do.call(mosaic, rast.list)
  
  plot(y_predvar)
  writeRaster(y_predvar, filename =paste0("./MDS_", compNumb, "_predvar95p.tif"), overwrite = TRUE)
  
  
  gc()
}

###

list_ras <- list.files("./",pattern="pred.tif$")
list_ras

pre.raster <- stack(paste0('./', list_ras))
pre.raster 

plot(pre.raster)

#### plot
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(tidyterra)
library(raster)
library(colorspace)
library(terra)
library(ggpubr)
library(viridisLite)

# pal1 <- diverge_hcl(100, h = c(246, 40), c = 96)
# pal1[ceiling(length(pal1) / 2)] <- "bisque"
# 
# pal2 <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
# pal2[ceiling(length(pal2) / 2)] <- "grey80"
# 
# pal3 <- diverge_hcl(100, h = c(130, 43), c = 98, l = c(70, 90))

pal1 <- rev(viridis(100))
pal2 <- rev(terrain.colors(100))
pal3 <- rev(colorRampPalette(rev(brewer.pal(11, "Spectral")))(100))

pals <- list(pal1,pal2,pal3,pal1,pal2)

###
ras <- rast(pre.raster)
names(ras)
ras <- ras[[c(3,5,4,2,1)]]
ras

components <- c("GC","Genome Size", "Gene Density","Cells","Cells(log10)")

g <- list()

for (i in 4:5){
g[[i]] <- ggplot()+
  geom_spatraster(data = ras[[i]],alpha=0.9) + 
 # facet_wrap(~lyr) +
  # scale_fill_distiller(palette = "Spectral",
  #                      na.value = "white")+
  scale_fill_gradientn( colours = pals[[i]],
                       # breaks=c(0,0.1,0.5,1,2,5,10,20,30,40,50), #,3,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85
                        na.value = NA, trans='sqrt',
                       # name = "Abundance(%)"
                      )+
  #scale_fill_gradientn()+
 ggtitle(components[i])+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))
}

ggarrange(g[[1]],g[[2]],g[[3]],ncol=3)

ggarrange(g[[4]],g[[5]],nrow=1)




