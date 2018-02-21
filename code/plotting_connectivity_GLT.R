####################################################
#
# Plotting results of connectivity assessment
# from LSMetrics
#
# Golden Lion Tamarin Landscape
# around Sao Joao river watershed, RJ, Brazil
#
# Bernardo Niebuhr - bernardo_brandaum at gmail.com
#
# No rights in the world are reserved
# Fell free to use, modify, and share
####################################################

# Loading libraries

# Working directories
# Data dir
datadir <- '/home/leecb/Documentos/UNESP/artigos/ms_Niebuhr_etal_2018_GLT_roads/data/calculated_metrics'

# Output dir
outdir <- '/home/leecb/Documentos/UNESP/artigos/ms_Niebuhr_etal_2018_GLT_roads/results'

# Go to data dir
setwd(datadir)

# Loading data
prefix <- "binary_forest_MLD_roads_river_mosaic"

# Patch size
patch <- read.table(paste(prefix, "_patch_AreaHA.txt", sep=""), header = T, sep = ",")

# Functional connectivity
nm <- list.files(pattern = "*clean_AreaHA.txt")

# Getting list of files and list of gap crossing values in order
a <- strsplit(nm, split = "_*m_")
b <- c()
for(i in 1:length(nm)) b <- c(b, strsplit(a[[i]][1], split="_"))
ll <- length(b[[1]])
aa <- c()
for(i in 1:length(nm)) aa <- c(aa, as.numeric(b[[i]][ll]))
aa
bb <- sort(aa) # order of gap crossing values
nm2 <- nm[order(aa)] # list of functional connectivity files in order

# Read table for each file - each gap crossing value
con <- list()
for(i in 1:length(nm2)) con[[i]] <- read.table(nm2[i], header = T, sep = ",")
names(con) <- bb

# Edges
nm <- list.files(pattern = "EDGE")
nm <- nm[grep("txt", nm)]

# Getting list of files and list of edge distance values in order
a <- strsplit(nm, split = "_*m.txt")
b <- c()
for(i in 1:length(nm)) b <- c(b, strsplit(a[[i]][1], split="EDGE"))
aa <- c()
for(i in 1:length(nm)) aa <- c(aa, as.numeric(b[[i]][2]))
aa
bb <- sort(aa) # order of edge distance values
nm2 <- nm[order(aa)] # edge distance files in order

# Read each file and get area/percentage of edge and core
edge <- matrix(, nrow = length(nm), ncol = 3)
for(i in 1:length(nm2)) {
  edge[i,] <- read.table(nm2[i], header = T, sep = ",")[,2]
}
rownames(edge) <- bb
colnames(edge) <- c("matrix", "edge", "interior")
edge

############################################################
# Plotting results

# Output dir
setwd(outdir)

pdf("plots_connectivity.pdf")

# Patch area
brk <- c(0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
n_lev <- hist(patch$HA, breaks = brk, plot = F)$counts
freq_lev <- n_lev/sum(n_lev)

#op <- par(mar = c(7, 4.5, 4, 2) + 0.1, mfrow=c(2,1))
op <- par(mar = c(8, 4.5, 4, 2) + 0.1)
barplot(freq_lev*100, las = 1, ylab = "Percentage of patches", xlab = "", 
        col = "grey", cex.lab = 1.7, cex.main = 1.5, 
        axes = FALSE, ylim = c(0, 80))
axis(2)
att <- seq(0.8, 0.8+1.2*(length(brk)-2), 1.2) #at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2, 10.4, 11.6)
lev <- c()
for(i in 1:(length(brk)-1)) lev <- append(lev, paste(brk[i],"-",brk[i+1],sep=""))
axis(1, at = att, labels = lev, las = 2)
mtext("Class of patch size (ha)", side = 1, line = 6, cex = 1.5)
par(op)

patch2 <- patch[order(patch$HA),]
patch2$cum.area <- cumsum(patch2$HA)
when.break <- sapply(brk[2:length(brk)], function(x) max(which(patch2$HA < x)))
#cum.area <- ifelse(freq_lev > 0, patch2$cum.area[when.break], 0)
cum.area <- patch2$cum.area[when.break]
cum.area.f <- cum.area[1]
for(i in 2:length(cum.area)) cum.area.f <- c(cum.area.f, cum.area[i]-cum.area[i-1])
cum.area.rel <- cum.area.f/sum(cum.area.f)

op <- par(mar = c(8, 4.5, 4, 2) + 0.1)
barplot(cum.area.f, las = 1, ylab = "Area (ha)", xlab = "",
        col = "grey", cex.lab = 1.7, cex.main = 1.5, ylim = c(0, max(cum.area.f)+5000),
        axes = FALSE)
axis(2)
axis(1, at = att, labels = lev, las = 2)
mtext("Class of patch size (ha)", side = 1, line = 6, cex = 1.5)
text(att-0.1,
     cum.area.f + 2000, paste(round(cum.area.rel, 2),"%A", sep = ""), cex = 0.7)
text(c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2, 10.4, 11.6)-0.1, 
     cum.area.f + 4000, paste(round(freq_lev, 2),"%NP", sep = ""), cex = 0.7)
par(op)

op <- par(mar = c(5, 8, 4, 2) + 0.1)
barplot(cum.area.f, las = 1, xlab = "Area (ha)", ylab = "", horiz = T,
        col = "grey", cex.lab = 1.7, cex.main = 1.5, xlim = c(0, max(cum.area.f)+9000),
        axes = FALSE)
axis(1)
axis(2, at = att, labels = lev, las = 2)
mtext("Class of patch size (ha)", side = 2, line = 6.5, cex = 1.5)
text(cum.area.f + 4300, att+0.2, 
     paste(round(cum.area.rel, 2),"%A", sep = ""), cex = 0.7)
text(cum.area.f + 6000, att-0.2, paste(round(freq_lev, 5),"%NP", sep = ""), cex = 0.7)
par(op)

# Connectivity

con
avg.cluster <- sapply(con, colMeans)[2,]
max.cluster <- sapply(con, function(x) sapply(x, max))[2,]
tot <- sum(patch[,2])
max.cluster.rel <- max.cluster/tot

plot(as.numeric(attr(avg.cluster, "names")), avg.cluster, pch=21, bg=1,
     xlab = "Functional distance (m)", ylab = "Expected cluster size (ha)")
lines(as.numeric(attr(avg.cluster, "names")), avg.cluster)

# para o maximo nao faz sentido fazer agora...
plot(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel, pch=21, bg=1,
     xlab = "Functional distance (m)", 
     ylab = "Expected cluster size\n(% of total remaining forest)")
lines(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel)

# Edges

edge <- as.data.frame(edge)
tot <- sum(edge[1,2:3])
edge$cum.area <- 100*edge$edge/tot

########################
# AQUI MUDAR O EIXO X (OS LABELS, PARA OUTROS VALORES!!!)
dists <- as.numeric(rownames(edge))
plot(dists, edge$cum.area, type = "b", log = "",
     xlab = "Edge distance (m)", ylab = "Cumulative area (%)")
plot(dists, edge$cum.area, type = "b", log = "x",
     xlab = "Edge distance (m)", ylab = "Cumulative area (%)")
segments(x0 = dists[4], y0 = edge$cum.area[1]-3, x1 = dists[4], y1 = edge$cum.area[4], lty = 2)
segments(x0 = dists[1]-5, y0 = edge$cum.area[4], x1 = dists[4], y1 = edge$cum.area[4], lty = 2)
segments(x0 = dists[9], y0 = edge$cum.area[1]-3, x1 = dists[9], y1 = edge$cum.area[9], lty = 2)
segments(x0 = dists[1]-5, y0 = edge$cum.area[9], x1 = dists[9], y1 = edge$cum.area[9], lty = 2)

edge$cum.area <- edge$cum.area/100
brk <- c(50, 100, 250, 500, 750, 1000)
when.break <- sapply(brk, function(x) max(which(dists < x)))
#cum.area <- ifelse(freq_lev > 0, patch2$cum.area[when.break], 0)
cum.area <- edge$cum.area[when.break]
cum.area.f <- cum.area[1]
for(i in 2:length(cum.area)) cum.area.f <- c(cum.area.f, cum.area[i]-cum.area[i-1])

#att <- seq(0.8, 0.8+1.2*(length(brk)-2), 1.2) #at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2, 10.4, 11.6)
#lev <- c()
#for(i in 1:(length(brk)-1)) lev <- append(lev, paste(brk[i],"-",brk[i+1],sep=""))
op <- par(mar = c(7, 4.5, 4, 2) + 0.1)
barplot(cum.area.f*100, las = 1, ylab = "% of remaining forest", xlab = "", 
        col = "grey", cex.lab = 1.7, cex.main = 1.5, 
        axes = FALSE)
axis(2)
axis(1, at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8), 
     labels = c("0-50", "50-100", "100-250", "250-500", "500-750", 
                "750-1000"), las = 2)
mtext("Edge distance (m)", side = 1, line = 5.5, cex = 1.5)
par(op)

dev.off()

############################################################
# Ploting together

pdf("ploting_connectivity2.pdf", width = 15, height = 15)
# png('plots_connectivity.png', width = 18, height = 18, units = 'cm', res = 500)

par(mfrow = c(2,2))

# Patch size
brk <- c(0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
n_lev <- hist(patch$HA, breaks = brk, plot = F)$counts
freq_lev <- n_lev/sum(n_lev)

att <- seq(0.8, 0.8+1.2*(length(brk)-2), 1.2) #at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2, 10.4, 11.6)
lev <- c()
for(i in 1:(length(brk)-1)) lev <- append(lev, paste(brk[i],"-",brk[i+1],sep=""))

patch2 <- patch[order(patch$HA),]
patch2$cum.area <- cumsum(patch2$HA)
when.break <- sapply(brk[2:length(brk)], function(x) max(which(patch2$HA < x)))
#cum.area <- ifelse(freq_lev > 0, patch2$cum.area[when.break], 0)
cum.area <- patch2$cum.area[when.break]
cum.area.f <- cum.area[1]
for(i in 2:length(cum.area)) cum.area.f <- c(cum.area.f, cum.area[i]-cum.area[i-1])
cum.area.rel <- cum.area.f/sum(cum.area.f)

op <- par(mar = c(5, 10, 4, 3) + 0.2)
barplot(cum.area.f, las = 1, xlab = "Area (ha)", ylab = "", horiz = T,
        col = "grey", cex.lab = 1.7, cex.axis = 3,
        xlim = c(0, max(cum.area.f)+11000),
        axes = FALSE)
axis(1, cex.axis = 1.8, lwd = 3)
axis(2, at = att, labels = lev, las = 2, cex.axis = 1.8, lwd = 3)
mtext("Class of patch size (ha)", side = 2, line = 8.5, at = 3.5, cex = 1.8)
text(cum.area.f + 4000, att+0.2, 
     paste(round(cum.area.rel, 2),"%A", sep = ""), cex = 1.2)
text(cum.area.f + 5700, att-0.2, paste(round(freq_lev, 5),"%NP", sep = ""), cex = 1.2)

mtext("A", side = 3, at = 0, cex = 2, )
par(op)


# Connectivity

con
avg.cluster <- sapply(con, colMeans)[2,]
max.cluster <- sapply(con, function(x) sapply(x, max))[2,]
tot <- sum(patch[,2])
max.cluster.rel <- max.cluster/tot

plot(as.numeric(attr(avg.cluster, "names")), avg.cluster, pch=19, cex = 2, lwd = 2,
     bg=1, cex.lab = 1.7,
     xlab = "Functional distance (m)", ylab = "", las=2, axes=F)
axis(1, cex.axis = 1.8, lwd = 3)
axis(2, las = 2, cex.axis = 1.8, lwd = 3)
mtext("Expected cluster size (ha)", side = 2, line = 5.5, cex=1.8)

lines(as.numeric(attr(avg.cluster, "names")), avg.cluster)

mtext("B", side = 3, at = 0, cex = 2)
# para o maximo nao faz sentido fazer agora...
# plot(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel, pch=21, bg=1,
#      xlab = "Functional distance (m)", 
#      ylab = "Expected cluster size\n(% of total remaining forest)")
# lines(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel)

# Edges

edge <- as.data.frame(edge)
tot <- sum(edge[1,2:3])
edge$cum.area <- 100*edge$edge/tot

########################
# AQUI MUDAR O EIXO X (OS LABELS, PARA OUTROS VALORES!!!)
op <- par(mar = c(5, 8, 4, 3) + 0.2)

dists <- as.numeric(rownames(edge))
#plot(dists, edge$cum.area, type = "b", log = "",
#     xlab = "Edge distance (m)", ylab = "Cumulative area (%)")
plot(dists, edge$cum.area, type = "b", log = "x", pch = 19, cex = 2, lwd = 2, 
     xlab = "Edge distance (m)", ylab = "", cex.lab = 1.7, axes = F, cex.axis = 1.5,
     ylim = c(0,100), xlim = c(25,1000))
axis(1, cex.axis = 1.8, at = c(25, 50, 100, 250, 500, 1000), lwd = 3)
axis(2, at = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100), cex.axis = 1.8, lwd = 3)
mtext("Cumulative area (%)", side = 2, line = 4.5, cex = 1.8)
segments(x0 = dists[4], y0 = edge$cum.area[1]-20, x1 = dists[4], y1 = edge$cum.area[4], lty = 2, lwd = 3)
segments(x0 = dists[1]-8, y0 = edge$cum.area[4], x1 = dists[4], y1 = edge$cum.area[4], lty = 2, lwd = 3)
segments(x0 = dists[9], y0 = edge$cum.area[1]-20, x1 = dists[9], y1 = edge$cum.area[9], lty = 2, lwd = 3)
segments(x0 = dists[1]-8, y0 = edge$cum.area[9], x1 = dists[9], y1 = edge$cum.area[9], lty = 2, lwd = 3)

mtext("C", side = 3, at = 25, cex = 2)
par(op)

dev.off()


####
# Plot as png

# pdf("ploting_connectivity2.pdf", width = 15, height = 15)
png('plots_connectivity_res150.png', width = 35, height = 35, units = 'cm', res = 150)

par(mfrow = c(2,2))

# Patch size
brk <- c(0, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
n_lev <- hist(patch$HA, breaks = brk, plot = F)$counts
freq_lev <- n_lev/sum(n_lev)

att <- seq(0.8, 0.8+1.2*(length(brk)-2), 1.2) #at = c(0.8, 2, 3.2, 4.4, 5.6, 6.8, 8, 9.2, 10.4, 11.6)
lev <- c()
for(i in 1:(length(brk)-1)) lev <- append(lev, paste(brk[i],"-",brk[i+1],sep=""))

patch2 <- patch[order(patch$HA),]
patch2$cum.area <- cumsum(patch2$HA)
when.break <- sapply(brk[2:length(brk)], function(x) max(which(patch2$HA < x)))
#cum.area <- ifelse(freq_lev > 0, patch2$cum.area[when.break], 0)
cum.area <- patch2$cum.area[when.break]
cum.area.f <- cum.area[1]
for(i in 2:length(cum.area)) cum.area.f <- c(cum.area.f, cum.area[i]-cum.area[i-1])
cum.area.rel <- cum.area.f/sum(cum.area.f)
freq_lev_100 <- ifelse(freq_lev < 0.001, '<0.1', round(freq_lev*100, 1))

op <- par(mar = c(5, 10, 4, 3) + 0.2)
barplot(cum.area.f, las = 1, xlab = "", ylab = "", horiz = T,
        col = "grey", cex.lab = 1.9, cex.axis = 3,
        xlim = c(0, max(cum.area.f)+11000),
        axes = FALSE)
axis(1, cex.axis = 1.8, lwd = 3)
axis(2, at = att, labels = lev, las = 2, cex.axis = 1.8, lwd = 3)
mtext("Class of patch size (ha)", side = 2, line = 8.5, at = 3.5, cex = 1.8)
mtext("Area (ha)", side = 1, line = 3, cex = 1.8)
text(cum.area.f + 5500, att+0.2, 
     paste(round(100*cum.area.rel, 1),"%A", sep = ""), cex = 1.7)
text(cum.area.f + 6600, att-0.3, paste(freq_lev_100, "%NP", sep = ""), cex = 1.7)

mtext("A", side = 3, at = 0, cex = 2)
par(op)


# Connectivity

con
avg.cluster <- sapply(con, colMeans)[2,]
max.cluster <- sapply(con, function(x) sapply(x, max))[2,]
tot <- sum(patch[,2])
max.cluster.rel <- max.cluster/tot

plot(as.numeric(attr(avg.cluster, "names")), avg.cluster, pch=19, cex = 2, lwd = 2,
     bg=1, cex.lab = 1.9,
     xlab = "", # xlab = "Functional distance (m)", 
     ylab = "", las=2, axes=F)
axis(1, cex.axis = 1.8, lwd = 3)
axis(2, las = 2, cex.axis = 1.8, lwd = 3)
mtext("Expected cluster size (ha)", side = 2, line = 5.5, cex=1.8)
mtext("Gap crossing distance (m)", side = 1, line = 3, cex = 1.8)

mtext("P", side = 2, line = 6, at = 14200, las = 2, cex = 1.6)
lines(as.numeric(attr(avg.cluster, "names")), avg.cluster, lwd = 3)

mtext("B", side = 3, at = 0, cex = 2)
# para o maximo nao faz sentido fazer agora...
# plot(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel, pch=21, bg=1,
#      xlab = "Functional distance (m)", 
#      ylab = "Expected cluster size\n(% of total remaining forest)")
# lines(as.numeric(attr(max.cluster.rel, "names")), 100*max.cluster.rel)

# Edges

edge <- as.data.frame(edge)
tot <- sum(edge[1,2:3])
edge$cum.area <- 100*edge$edge/tot

########################
# AQUI MUDAR O EIXO X (OS LABELS, PARA OUTROS VALORES!!!)
op <- par(mar = c(5, 8, 4, 3) + 0.2)

dists <- as.numeric(rownames(edge))
#plot(dists, edge$cum.area, type = "b", log = "",
#     xlab = "Edge distance (m)", ylab = "Cumulative area (%)")
plot(dists, edge$cum.area, type = "b", log = "x", pch = 19, cex = 2, lwd = 2, 
     xlab = "", ylab = "", cex.lab = 1.9, axes = F, cex.axis = 1.5,
     ylim = c(0,100), xlim = c(25,1000))
axis(1, cex.axis = 1.8, at = c(25, 50, 100, 250, 500, 1000), lwd = 3)
axis(2, at = c(0,20,40,60,80,100), labels = c(0,20,40,60,80,100), cex.axis = 1.8, lwd = 3)
mtext("Cumulative area (%)", side = 2, line = 4.5, cex = 1.8)
mtext("Edge distance (m)", side = 1, line = 3, cex = 1.8)
segments(x0 = dists[4], y0 = edge$cum.area[1]-20, x1 = dists[4], y1 = edge$cum.area[4], lty = 2, lwd = 3)
segments(x0 = dists[1]-8, y0 = edge$cum.area[4], x1 = dists[4], y1 = edge$cum.area[4], lty = 2, lwd = 3)
segments(x0 = dists[9], y0 = edge$cum.area[1]-20, x1 = dists[9], y1 = edge$cum.area[9], lty = 2, lwd = 3)
segments(x0 = dists[1]-8, y0 = edge$cum.area[9], x1 = dists[9], y1 = edge$cum.area[9], lty = 2, lwd = 3)

mtext("C", side = 3, at = 25, cex = 2)
par(op)

dev.off()
