####Least Cost Path
library(raster)
library(gdistance)
#
#
#Need to transform raster to transitionlayer
r1=raster("~/SDM/Pathway/lano/lano_fall.tif")
values(r1)[values(r1) < 0] = NA
tr_fall <- transition(r1, transitionFunction=mean, directions=8)

r2=raster("~/SDM/Pathway/lano/lano_spring.tif")
values(r2)[values(r2) < 0] = NA
tr_spring <- transition(r2, transitionFunction=mean, directions=8)

summer_loc=read.csv("~/SDM/Pathway/lano/lano_summer_loc.csv",header=FALSE)
winter_loc=read.csv("~/SDM/Pathway/lano/lano_winter_loc.csv",header=FALSE)

sp_summer=SpatialPoints(summer_loc[,1:2])
sp_winter=SpatialPoints(winter_loc[,1:2])

####Creating null maps for desnity
not=raster(ext=extent(r2),res=0.5)
r1_a=resample(r1,not,method="bilinear")
r2_a=resample(r2,not,method="bilinear")

tr_fall_a <- transition(r1_a, transitionFunction=mean, directions=8)
tr_spring_a <- transition(r2_a, transitionFunction=mean, directions=8)

AtoB_density=raster(res=res(r2_a),ext=extent(r2_a))
BtoA_density=raster(res=res(r1_a),ext=extent(r1_a))

values(AtoB_density)=0
values(BtoA_density)=0


####Spring
plot(raster(tr_spring))
print(nrow(winter_loc))
print(nrow(summer_loc))
for(i in 1:nrow(summer_loc)){
  for(j in 1:nrow(winter_loc)){
    AtoB=shortestPath(tr_spring_a,sp_winter[j,],sp_summer[i,],output="SpatialLines")
    lines(AtoB,col='black')
    
    AtoB_a=shortestPath(tr_spring_a,sp_winter[j,],sp_summer[i,],output="TransitionLayer")
    AtoB_b=raster(AtoB_a)
    AtoB_b[is.na(AtoB_b[])] <- 0 
    AtoB_density=calc(stack(AtoB_density,AtoB_b),fun=sum)
    AtoB_density[is.na(AtoB_density[])] <- 0 
    

  }
  print("summer #")
  print(i)
}



#####Fall
plot(raster(tr_fall))
print(nrow(winter_loc))
print(nrow(summer_loc))

for(i in 1:nrow(winter_loc)){
  for(j in 1:nrow(summer_loc)){
    BtoA=shortestPath(tr_fall_a,sp_summer[j,],sp_winter[i,],output="SpatialLines")
    lines(BtoA,col='black',lwd=1.5)
    
    BtoA_a=shortestPath(tr_fall_a,sp_summer[j,],sp_winter[i,],output="TransitionLayer")
    BtoA_b=raster(BtoA_a)
    BtoA_b[is.na(BtoA_b[])] <- 0 
    BtoA_density=calc(stack(BtoA_density,BtoA_b),fun=sum)
    BtoA_density[is.na(BtoA_density[])] <- 0 
    

  }
  print("winter #")
  print(i)
}


