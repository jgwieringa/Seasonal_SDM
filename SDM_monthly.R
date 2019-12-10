####For Loop SDM
library(randomForest)
library(raster)
library(rgeos)
library(maptools)
library(dismo)
library(kernlab)
library(sp)
library(ecospat)
library(ENMeval)
library(maps)

#####Get gbif data ready
###Data downloaded from gbif
####sorted into months
####REad in each month and then reduce sampling bias
ep12=read.csv("~/SDM/Monthly/Dec/epfu_dec.csv")

ep12=ep12[,-1]
ep12=ep12[ep12$lon < -52 & ep12$lat > 15, ]

ep12a=gridSample(ep12,r,n=1)

###re-write csv files
write.csv(ep11a,"~/SDM/Monthly/Nov/epfu_nov.csv")


###Sets the java path so rJava runs properly
Sys.setenv(JAVA_HOME='C://Program Files/Java/jre1.8.0_211')
library(rJava)

ext=extent(-180,-51.5,15,80)

#Creating data matrix to collect data
a1=matrix(nrow=36,ncol=16)
a2=matrix(nrow=36,ncol=25)
colnames(a1)=c('Month','Speices','AUC_rf','AUC_gl','AUC_mx','AUC_bc','thr_AUC','AUC_AUCw_en','AUC_tssw_en','tss_rf','tss_gl','tss_mx','tss_bc','thr_tss','tss_AUCw_en','tss_tssw_en')
colnames(a2)=c("Month","Species","Prec","srad","tavg","tmax","tmin","vapr","wind","humam_inf","elev","NDVI","forest","gl_intercept","Prec","srad","tavg","tmax","tmin","vapr","wind","humam_inf","elev","NDVI","forest")

###Be in folder that contains the seasonal folders
folders=list.dirs(path='.')

month=c("apr","aug","dec","feb","jan","jul","jun","mar","may","nov","oct","sep")

for(j in 2:13){
  
  files1 <- list.files(path=folders[j],pattern = '.tif$', full.names = TRUE)
  pred=stack(files1)
  
  
  files <- list.files(path=folders[j],pattern = '.csv$', full.names = TRUE)
  
  species=c("labo","laci","lano")
  
  for(i in 1:length(files)){
    loc=read.csv(files[i])
    loc=loc[,-1]
    pres=extract(pred, loc)
    set.seed(0)
    backgr=randomPoints(pred, 500)
    
    ###sudo absence points
    absvals=extract(pred, backgr)
    pb=c(rep(1, nrow(pres)), rep(0, nrow(absvals)))
    sdmdata=data.frame(cbind(pb, rbind(pres, absvals)))
    
    group <- kfold(loc, 5)
    pres_train=loc[group != 1, ]
    pres_test=loc[group == 1, ]
    
    ###BioClim Model
    bc=bioclim(pred, pres_train)
    backg=randomPoints(pred, n=1000, ext=ext, extf = 1.25)
    colnames(backg) = c('Longitude', 'Latitude')
    group=kfold(backg, 5)
    backg_train=backg[group != 1, ]
    backg_test=backg[group == 1, ]
    
    ####Actual prediction from BC
    pb=predict(pred, bc, ext=ext, progress='')
    
    
    ####maxent
    ###Need Java installed
    ###Need MaxENT in the dismo folder
    jar=paste(system.file(package="dismo"), "~/MaxENT/maxent.jar", sep='')
    xm=maxent(pred, pres_train)
    px=dismo::predict(pred, xm, ext=ext, progress='')
    
    ####RAndom Forest
    colnames(backg_train)=c("lon","lat")
    train=rbind(pres_train, backg_train)
    pb_train=c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
    envtrain=extract(pred, train)
    envtrain=data.frame( cbind(pa=pb_train, envtrain) )
    envtrain=na.omit(envtrain)
    model=factor(pa)~prec+srad+tavg+tmax+tmin+vapr+wind+human_inf+elev+NDVI+forest
    rf1=randomForest(model, data=envtrain)
    pr=dismo::predict(pred, rf1, ext=ext)
    
    ####GLM
    glm1=glm(pa~prec+srad+tavg+tmax+tmin+vapr+wind+human_inf+elev+NDVI+forest, data = envtrain)
    pg=dismo::predict(pred,glm1,ext=ext)
    
    
    ###Weighted Mean
    emx=evaluate(pres_test, backg_test, xm, pred)
    erf=evaluate(pres_test, backg_test, rf1, pred)
    egl=evaluate(pres_test, backg_test, glm1, pred)
    ebc=evaluate(pres_test, backg_test, bc, pred)
    
    auc <- sapply(list(erf, egl, emx, ebc), function(x) x@auc)
    w <- (auc-0.5)^2
    models=stack(pr,pg,px,pb)
    m1 <- weighted.mean(models, w)
    
    #tss weighted mean
    colnames(backgr)=c("lon","lat")
    tss=rbind(loc,backgr)
    tss_pts=SpatialPoints(cbind(tss[,1],tss[,2]))
    tss2=cbind(tss,sdmdata$pb)
    
    tss.mx_values=extract(px,tss_pts,method='simple')
    tss.glm_values=extract(pg,tss_pts,method='simple')
    tss.rf_values=extract(pr,tss_pts,method='simple')
    tss.bc_values=extract(pb,tss_pts,method='simple')
    
    tss2=cbind(tss,sdmdata$pb,tss.mx_values,tss.bc_values,tss.glm_values,tss.rf_values)
    tss_final=na.omit(tss2)
    
    tss.mx=ecospat.max.tss(tss_final[,4],tss_final[,3])
    tss.bc=ecospat.max.tss(tss_final[,5],tss_final[,3])
    tss.glm=ecospat.max.tss(tss_final[,6],tss_final[,3])
    tss.rf=ecospat.max.tss(tss_final[,7],tss_final[,3])
    
    tss_values=c(tss.rf[[2]][1,2],tss.glm[[2]][1,2],tss.mx[[2]][1,2],tss.bc[[2]][1,2])
    tss_values2=as.numeric(tss_values[1:4])
    models=stack(pr,pg,px,pb)
    m2 <- weighted.mean(models, tss_values2)
    
    #####Evaluating ensemble models
    tss.en_auc_values=extract(m1,tss_pts,method='simple')
    tss.en_tss_values=extract(m2,tss_pts,method='simple')
    tss.en_p1=cbind(tss,sdmdata$pb,tss.en_auc_values,tss.en_tss_values)
    tss.en_p2=na.omit(tss.en_p1)
    tss.en.auc=ecospat.max.tss(tss.en_p2[,4],tss.en_p2[,3])
    tss.en.tss=ecospat.max.tss(tss.en_p2[,5],tss.en_p2[,3])
    
    
    tss_en_no_neg=tss.en_p2
    tss_en_no_neg[,3]=ifelse(tss_en_no_neg[,3]<0,0,tss_en_no_neg[,3])
    tss_en_no_neg[,4]=ifelse(tss_en_no_neg[,4]<0,0,tss_en_no_neg[,4])
    tss_en_no_neg[,5]=ifelse(tss_en_no_neg[,5]<0,0,tss_en_no_neg[,5])
    
    auc.en.auc=SDMTools::auc(tss_en_no_neg[,3],tss_en_no_neg[,4])
    auc.en.tss=SDMTools::auc(tss_en_no_neg[,3],tss_en_no_neg[,5])
    
    
    ###Generate threshold
    tr1=ecospat.mpa(m1,pres_test,perc = 0.9)
    tr2=ecospat.mpa(m2,pres_test,perc = 0.9)
    
    ###Generate p/np raster
    m1_tr=m1 > tr1
    m2_tr=m2 > tr2
    
    ###Writing four rasters
    ###Two for each time period
    ###Two continous and two binary

    ###Add data to matrix
    a1[j-1+((i-1)*12),1]=month[j-1]
    a1[j-1+((i-1)*12),2]=species[i]
    a1[j-1+((i-1)*12),3:6]=w
    a1[j-1+((i-1)*12),7]=tr1
    a1[j-1+((i-1)*12),8]=auc.en.auc-0.5
    a1[j-1+((i-1)*12),9]=auc.en.tss-0.5
    a1[j-1+((i-1)*12),10:13]=tss_values2
    a1[j-1+((i-1)*12),14]=tr2
    a1[j-1+((i-1)*12),15]=tss.en.auc[[2]][1,2]
    a1[j-1+((i-1)*12),16]=tss.en.tss[[2]][1,2]
    
    print(species[i])
    
    ###variable importance matrix
    a2[j-1+((i-1)*12),1]=month[j-1]
    a2[j-1+((i-1)*12),2]=species[i]
    a2[j-1+((i-1)*12),3:13]=importance(rf1)
    a2[j-1+((i-1)*12),14:25]=glm1$coefficients
    
    xm

  }
  
  print(month[j-1])
  
}

write.csv(a1,'~/SDM/Monthly/stats2.csv')
write.csv(a2,'~/SDM/Monthly/var_weights.csv')

#######Combing each species rasters into one stack
lano_jan_pnp=raster("lano_jan_tss_p_np.tif")
lano_feb_pnp=raster("lano_feb_tss_p_np.tif")
lano_mar_pnp=raster("lano_mar_tss_p_np.tif")
lano_apr_pnp=raster("lano_apr_tss_p_np.tif")
lano_may_pnp=raster("lano_may_tss_p_np.tif")
lano_jun_pnp=raster("lano_jun_tss_p_np.tif")
lano_jul_pnp=raster("lano_jul_tss_p_np.tif")
lano_aug_pnp=raster("lano_aug_tss_p_np.tif")
lano_sep_pnp=raster("lano_sep_tss_p_np.tif")
lano_oct_pnp=raster("lano_oct_tss_p_np.tif")
lano_nov_pnp=raster("lano_nov_tss_p_np.tif")
lano_dec_pnp=raster("lano_dec_tss_p_np.tif")

lano_pnp=stack(lano_jan_pnp,lano_feb_pnp,lano_mar_pnp,lano_apr_pnp,lano_may_pnp,lano_jun_pnp,lano_jul_pnp,lano_aug_pnp,lano_sep_pnp,lano_oct_pnp,lano_nov_pnp,lano_dec_pnp)

names(lano_pnp)=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

actual=stack(wind_thr,lano_pnp)

#######Correlation between wind density and probabilities
wind=raster("windturbine_density.tif")
loc=read.csv("windfarm_loc.csv")
tr_wind=ecospat.mpa(wind,loc,perc = 0.9)
wind_thr=wind>tr_wind

#####Create stack
cor_stack=stack(avg_thr,lano_pnp)
cor=layerStats(cor_stack,'pearson',na.rm=T)


a=raster("./Jan/wind.tif")
b=raster("./Feb/wind.tif")
c=raster("./Mar/wind.tif")
d=raster("./Apr/wind.tif")
e=raster("./May/wind.tif")
f=raster("./Jun/wind.tif")
g=raster("./Jul/wind.tif")
h=raster("./Aug/wind.tif")
i=raster("./Sep/wind.tif")
j=raster("./Oct/wind.tif")
k=raster("./Nov/wind.tif")
l=raster("./Dec/wind.tif")

avg=(a+b+c+d+e+f+g+h+i+j+k+l)/12

tr_avg=ecospat.mpa(avg,loc,perc = 0.9)
avg_thr=avg>tr_avg


