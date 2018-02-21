# Analysis Genetics vs roads GLT

# Folders
datadir <- "/home/leecb/Documentos/UNESP/artigos/ms_Niebuhr_etal_2018_GLT_roads/data"

# Read data
setwd(datadir)

gen.roads <- read.table("KinshipGen_Roads.csv", header = T, sep = ";")
head(gen.roads, 20)
str(gen.roads)

# Organize data
gen.roads[,c("Management", "Roads", "Dirt", "Paved", "BR.101")] <- 
  as.data.frame(lapply(gen.roads[,c("Management", "Roads", "Dirt", "Paved", "BR.101")], 
                       FUN = function(X) { as.factor(as.character(X))} ))

dist <- read.table("Moraes_etal_data_landscape_genetics_Leontopithecus_rosalia.csv", header = T)
head(dist)
nrow(dist)

gen.roads$Euclidean_distance <- dist$Euclidean_distance

# Exploratory plots
plot(gen.roads$Fij ~ gen.roads$Management)
plot(gen.roads$Fij ~ gen.roads$Roads)
plot(gen.roads$Fij ~ gen.roads$Dirt)
plot(gen.roads$Fij ~ gen.roads$Paved)
plot(gen.roads$Fij ~ gen.roads$BR.101)
plot(gen.roads$Fij ~ gen.roads$Region)
plot(gen.roads$Fij ~ gen.roads$Euclidean_distance)

m1 <- glm(gen.roads$Fij ~ gen.roads$Management + gen.roads$Euclidean_distance)
res <- resid(m1)

m2 <- glm(gen.roads$Fij ~ gen.roads$Management + gen.roads$Euclidean_distance + gen.roads$road.type - 1)
summary(m2)
m3 <- glm(gen.roads$Fij ~ gen.roads$Management * gen.roads$road.type - 1+ gen.roads$Euclidean_distance)
summary(m3)

plot(res ~ gen.roads$Roads)
plot(res ~ gen.roads$Dirt)
plot(res ~ gen.roads$Paved)
plot(res ~ gen.roads$BR.101)
plot(res ~ gen.roads$Region)

summary(glm(res ~ gen.roads$Roads))
summary(glm(res ~ gen.roads$Dirt))
summary(glm(res ~ gen.roads$Paved))
summary(glm(res ~ gen.roads$BR.101))
summary(glm(res ~ gen.roads$Region))

# Classify road type
gen.roads$road.type <- 0
gen.roads$road.type <- ifelse(gen.roads$Dirt == 1, 1, gen.roads$road.type)
gen.roads$road.type <- ifelse(gen.roads$Paved == 1, 2, gen.roads$road.type)
gen.roads$road.type <- ifelse(gen.roads$BR.101 == 1, 3, gen.roads$road.type)
gen.roads$road.type <- as.factor(gen.roads$road.type)
plot(gen.roads$Fij ~ gen.roads$road.type)
plot(res ~ gen.roads$road.type)

summary(glm(gen.roads$Fij ~ gen.roads$road.type-1))
summary(glm(res ~ gen.roads$road.type-1))

# Test each management type separetely
# Test interaction road type + management

table(gen.roads$road.type)/sum(table(gen.roads$road.type))

# But we should consider types of roads that happen exclusively
gen.roads$only_paved <- as.factor(as.numeric(gen.roads$Paved) - as.numeric(gen.roads$BR.101))

gen.roads$road.type2 <- 0
gen.roads$road.type2 <- ifelse(gen.roads$Dirt == 1 & gen.roads$Paved == 0, 1, gen.roads$road.type2)
gen.roads$road.type2 <- ifelse(gen.roads$only_paved == 1 & gen.roads$BR.101 == 0 & gen.roads$Dirt == 0, 2, gen.roads$road.type2)
gen.roads$road.type2 <- ifelse(gen.roads$BR.101 == 1 & gen.roads$Dirt == 0, 3, gen.roads$road.type2)
gen.roads$road.type2 <- ifelse(gen.roads$only_paved == 1 & gen.roads$Dirt == 1, 4, gen.roads$road.type2)
gen.roads$road.type2 <- ifelse(gen.roads$BR.101 == 1 & gen.roads$Dirt == 1, 5, gen.roads$road.type2)
gen.roads$road.type2 <- as.factor(gen.roads$road.type2)
plot(gen.roads$Fij ~ gen.roads$road.type2)
plot(res ~ gen.roads$road.type2)

vioplot::vioplot()

summary(glm(gen.roads$Fij ~ gen.roads$road.type2-1))
summary(glm(res ~ gen.roads$road.type2-1))
aov(res ~ gen.roads$road.type2)

table(gen.roads$road.type2)
sum(table(gen.roads$road.type2))

par(mfrow=c(2,1))
mu<-2
si<-0.6
bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si)) 
uniform<-runif(2000,-4,4)
normal<-rnorm(2000,0,3)
library(vioplot)
vioplot::vioplot(bimodal,uniform,normal)
boxplot(bimodal,uniform,normal, plot = F, )

vioplot(gen.roads$Fij[gen.roads$road.type2 == 0],
        gen.roads$Fij[gen.roads$road.type2 == 1],
        gen.roads$Fij[gen.roads$road.type2 == 2],
        gen.roads$Fij[gen.roads$road.type2 == 3],
        gen.roads$Fij[gen.roads$road.type2 == 4],
        gen.roads$Fij[gen.roads$road.type2 == 5])

vioplot(res[gen.roads$road.type2 == 0],
        res[gen.roads$road.type2 == 1],
        res[gen.roads$road.type2 == 2],
        res[gen.roads$road.type2 == 3],
        res[gen.roads$road.type2 == 4],
        res[gen.roads$road.type2 == 5])
