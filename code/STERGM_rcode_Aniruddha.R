
library(statnet)
library(ndtv)
library(htmlwidgets)
library(latticeExtra)



data(samplk)
ls()


samplist <- list(samplk1,samplk2,samplk3)
sampdyn <- networkDynamic(network.list = samplist)


sampdyn

par(mfrow = c(2,2), oma=c(1,1,1,1), mar=c(4,1,1,1))
plot(network.extract(sampdyn, at = 0), main = "Time 1", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(sampdyn, at = 1), main = "Time2", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(sampdyn, at = 2), main = "Time3", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(sampdyn, main = "Collapsed", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)

data(short.stergm.sim)

par(mfrow = c(2,2), oma=c(1,1,1,1), mar=c(4,1,1,1))
plot(network.extract(short.stergm.sim, at = 0), main = "Time 1", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(short.stergm.sim, at = 5), main = "Time 5", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(short.stergm.sim, at = 15), main = "Time 20", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(short.stergm.sim, at = 24), main = "Time 30", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)



data(windsurfers)

par(mfrow = c(2,2), oma=c(1,1,1,1), mar=c(4,1,1,1))
plot(network.extract(windsurfers, at = 0), main = "Time 1", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(windsurfers, at = 10), main = "Time 5", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(windsurfers, at = 20), main = "Time 15", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)
plot(network.extract(windsurfers, at = 30), main = "Time 25", 
     displaylabels = T, label.cex = 0.6, vertex.cex = 2, pad = 0.5)

tSnaStats(sampdyn,"degree") # Changes in degree centrality




tErgmStats(sampdyn, "~ edges+triangle") # Notice the increase in triangles



library(ndtv)
data(short.stergm.sim)



render.d3movie(short.stergm.sim, 
               plot.par=list(displaylabels=T))


render.d3movie(short.stergm.sim,
               plot.par=list(displaylabels=T),
               output.mode = 'htmlWidget') # using htmlwidgets package here


proximity.timeline(short.stergm.sim,default.dist = 6,
                   mode = 'sammon',labels.at = 17,vertex.cex = 4)



stergm(my.network,                               #do not run this
       formation =  ~ edges + nodefactor('age'),
       dissolution =  ~ edges + gwesp(0, fixed=T),
       estimate =  `insert method`
)


samp.fit <- stergm(short.stergm.sim,
                   formation =  ~edges+mutual+gwesp(0,fixed = T),
                   dissolution = ~edges+gwesp(0,fixed = T),
                   estimate = "CMLE",
                   times = c(1:25)
)

summary(samp.fit)






samp.fit1 <- stergm(samplist,
                   formation =  ~mutual+cyclicalties+gwesp(0,fixed=T),
                   dissolution = ~mutual+cyclicalties+gwesp(0,fixed=T),
                   estimate = "CMLE",
                   times = c(1:3)
)
summary(samp.fit1)




###########################33

#model 1:
samp.fit <- stergm(samplist,
                   formation =  ~edges+mutual+cyclicalties+transitiveties,
                   dissolution = ~edges+mutual+cyclicalties+transitiveties,
                   estimate = "CMLE",
                   times = c(1:3)
)

summary(samp.fit)




samp.fit1 <- stergm(samplist,
                   formation =  ~edges+mutual+transitiveties,
                   dissolution = ~edges+cyclicalties+transitiveties,
                   estimate = "CMLE",
                   times = c(1:3)
)

summary(samp.fit1)


samp.fit2 <- stergm(samplist,
                    formation = ~edges + gwesp(0, fixed=T),
                    dissolution = ~edges + gwesp(0, fixed=T),
                    estimate = "CMLE",
                    times = c(1:3)
)

summary(samp.fit2)



samp.fit3 <- stergm(samplist,
                    formation = ~edges + gwesp(0, fixed=T)+mutual+transitiveties,
                    dissolution = ~edges + gwesp(0, fixed=T)+mutual+transitiveties,
                    estimate = "CMLE",
                    times = c(1:3)
)

summary(samp.fit3)


samp.fit4 <- stergm(windsurfer,
                    formation = ~edges + gwesp(0, fixed=T)+mutual+cyclicalties,
                    dissolution = ~edges + gwesp(0, fixed=T)+mutual+cyclicalties,
                    estimate = "CMLE",
                    times = c(1:3)
)

summary(samp.fit4)


samp.fit5 <- stergm(samplist,
                    formation = ~edges + gwesp(0, fixed=T)+mutual+cyclicalties+transitiveties,
                    dissolution = ~edges + gwesp(0, fixed=T)+mutual+cyclicalties+transitiveties,
                    estimate = "CMLE",
                    times = c(1:3)
)

summary(samp.fit5)

summary(samp.fit)
summary(samp.fit1)
summary(samp.fit2)
summary(samp.fit3)
summary(samp.fit4)
summary(samp.fit5)









#model 2
data(windsurfers)
tErgmStats(windsurfers, "~ edges+triangle")
tSnaStats(windsurfers,"degree")



samp.fit <- stergm(windsurfers,
                   formation =  ~edges,
                   dissolution = ~edges,
                   estimate = "CMLE",
                   times = c(1:31)
)

summary(samp.fit)



samp.fit5 <- stergm(windsurfers,
                    formation = ~edges + gwesp(0, fixed=T),
                    dissolution = ~edges + gwesp(0, fixed=T),
                    estimate = "CMLE",
                    times = c(1:31)
)

summary(samp.fit5)

summary(samp.fit)
summary(samp.fit5)




######model 3
short.stergm.sim



data(short.stergm.sim)
tErgmStats(short.stergm.sim, "~ edges+triangle")
tSnaStats(windsurfers,"degree")



samp.fit <- stergm(short.stergm.sim,
                   formation =  ~edges,
                   dissolution = ~edges,
                   estimate = "CMLE",
                   times = c(1:25)
)

summary(samp.fit)



samp.fit5 <- stergm(short.stergm.sim,
                    formation = ~edges + gwesp(0, fixed=T),
                    dissolution = ~edges ,
                    estimate = "CMLE",
                    times = c(1:25)
)

summary(samp.fit5)

summary(samp.fit)
summary(samp.fit5)





##########

data(McFarland_cls33_10_16_96)




data(McFarland_cls33_10_16_96)
tErgmStats(McFarland_cls33_10_16_96, "~ edges+triangle")
tSnaStats(McFarland_cls33_10_16_96,"degree")



samp.fit <- stergm(cls33_10_16_96,
                   formation =  ~edges,
                   dissolution = ~edges,
                   estimate = "CMLE",
                   times = c(1:49)
)

summary(samp.fit)



samp.fit5 <- stergm(cls33_10_16_96,
                    formation = ~edges + gwesp(0, fixed=T),
                    dissolution = ~edges + gwesp(0, fixed=T),
                    estimate = "CMLE",
                    times = c(1:49)
)

summary(samp.fit5)

summary(samp.fit)
summary(samp.fit5)

library(scatterplot3d)

compute.animation(sampdyn)
compute.animation(sampdyn)
par(mfrow = c(1,1))
timePrism(sampdyn,at=c(0,1,2),
          displaylabels=TRUE,planes = TRUE,
          label.cex=0.5)



