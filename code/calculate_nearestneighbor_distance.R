

# from 
# https://stackoverflow.com/questions/34752674/measuring-distances-between-every-patch-using-r

rm(list = ls())

library(raster)
library(rgeos)
library(rgdal)

#dir <- '/home/leecb/Documentos/Academico/artigos/ms_Ascensao_GLT_roads_ConsBiol/maps'
dir <- 'D:/regenerabilidade_RJ/Connectivity_RJ'
setwd(dir)

mapa <- readOGR('.', layer = 'ffcover01_under700m_res20m_HABMAT_pid_buff5m_shape')
#mapa <- raster('fcover01_under700m_res20m_HABMAT.tif')
mapa
mapa$value
# table(mapa[]) # ok, only 0 and 1
plot(mapa)

#mapa.shape <- rasterToPolygons(mapa, dissolve=TRUE)
NNdistances <- gDistance(mapa, byid=T)

diag(NNdistances) <- NA
NNdistances[1:3,1:10]
NNdist.perpatch <- apply(NNdistances, MARGIN=1, min, na.rm = T)
n <- nrow(NNdistances)
NN.2.dist.perpatch <- apply(NNdistances, MARGIN=1, function(x) sort(x)[2])
NN.3.dist.perpatch <- apply(NNdistances, MARGIN=1, function(x) sort(x)[3])
avg.dist.perpatch <- apply(NNdistances, MARGIN=1, mean, na.rm = T)



write.table(NNdist.perpatch, 'NNdist_perpatch.csv', sep = '\t')

head(NNdist.perpatch)

mapa$NNdist <- NNdist.perpatch
mapa$NNdist_2 <- NN.2.dist.perpatch
mapa$NNdist_3 <- NN.3.dist.perpatch
mapa$avg.dist <- avg.dist.perpatch
head(mapa@data)
writeOGR(mapa, '.', 'landscape_GLT_dist_per_patch', driver = 'ESRI Shapefile', overwrite_layer = T)

mean(mapa$NNdist)
hist(mapa$NNdist, breaks = 20, main = '')
quantile(mapa$NNdist, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
quantile(mapa$NNdist, p = c(0.025, 0.1, 0.5, 0.84, 0.975))

mean(mapa$NNdist_2)
hist(mapa$NNdist_2, breaks = 20, main = '')
quantile(mapa$NNdist_2, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
quantile(mapa$NNdist_2, p = c(0.025, 0.1, 0.5, 0.84, 0.975))

mean(mapa$NNdist_3)
hist(mapa$NNdist_3, breaks = 20, main = '')
quantile(mapa$NNdist_3, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
quantile(mapa$NNdist_3, p = c(0.025, 0.1, 0.23, 0.5, 0.84, 0.975))

mean(mapa$avg.dist)
hist(mapa$avg.dist, breaks = 20, main = '')
quantile(mapa$avg.dist, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
quantile(mapa$avg.dist, p = c(0.025, 0.1, 0.5, 0.84, 0.975))
