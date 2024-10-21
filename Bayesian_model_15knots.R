##before running the script
##### go to >>Session >> Set Working Directory >> To source file location

#### CLEAR MEMORY
rm(list=ls())
#### Prevent Warnings
options(warn=-1)
#### Clear the R console
cat("\f")



# 1. ## SOURCE FUNCTIONS AND READ IN LIBRARIES ----------------------------

time_elapsed=0
ptm1 = proc.time()
mainDir=getwd()
setwd(mainDir)
suppressPackageStartupMessages(source("Library.R"))
##### 1. ## SOURCE FUNCTIONS AND READ IN LIBRARIES
source("Library.R")

dir_create("Figures_cali")#DIRECTORY TO SAVE FIGURES
dir_fig=paste(mainDir,"Figures_cali",sep = "/")

dir_data=paste(mainDir,"data",sep = "/")#DATA DIRECTORY


load_packages(c("parallel","rstan","gamlss","data.table","rlist","akima","MASS","RColorBrewer","maps","SpatialExtremes","locfit","fields","seqinr","reshape2","ggplot2","dplyr"),quietly=TRUE)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

myPalette1 <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPalette2 <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPalette3 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")), space="Lab")


# 2. READ DATA ------------------------------------------------------------

setwd(dir_data)#set the data directory
#proxy cores
EP_cores=list.load("Eastern_Pacific_cores.rds")
WP_cores=list.load("Western_Pacific_cores.rds")
CP_core=read.delim("Line_Islands_core.txt", header = TRUE, sep = "\t", dec = ".")

Years=1855:2014#TOTAL COTEMPORARY TIME PERIOD
Years_cal=1855:2014#MODIFY THE LOWER LIMIT TO 1970 IF YOU WANT TO RUN THE VALIDATION 
Nk=15# NUMBER OF KNOTS. YOU CAN MODIFY FROM 3 TO 24
ind_t=1:length(Years)
ind_y_cal=which(Years %in% Years_cal)
T_c=length(ind_y_cal)
if (T_c<length(ind_t)){
  ind_y_val=ind_t[-ind_y_cal]
  T_v=length(ind_y_val)
}else{
  ind_y_val=ind_y_cal
  T_v=length(ind_y_val)
}


#READ THE KNOT LOCATIONS FILE
knots=list.load("knots.rds")[1:Nk,]
#READ THE LOCATIONS OF FITTING POINTS
index2_i=list.load("Fitting_points.rds")
#READ IN CONTEMPORARY DATA 
SST=list.load("1854-2014_MAYtoAPRavg_SST_FullField.rds")
xgrid = SST$xgrid
nx = SST$nx
ygrid = SST$ygrid
ny = SST$ny
grid <- SST$grid
index1 <- SST$index1
grid1=grid[index1,]#LOCATIONS OF THE 972 POINTS
nsta <- SST$nsta
climdata = SST$climdata
data = SST$data
data.un=SST$data.un


# 3. ## PCA ON FULL CONTEMPORARY SST FIELD --------------------------------

zs <- var(data)
zsvd <- svd(zs)
pcs <- t(t(zsvd$u) %*% t(data))
lambdas <- (zsvd$d/sum(zsvd$d))



# 4. ## SMOOTH AND PLOT PALEO RECONSTRUCTIONS -----------------------------

sst.smooth = gcv.recon(EP_cores$cores, WP_cores$cores,EPdeg=2, WPdeg=2, base=TRUE)

aa = 10 #10 ka
bb = 0 #0 ka
CP_core[,1] = CP_core[,1]/1000

zfit = locfit(CP_core[,2] ~ CP_core[,1], deg=2, kern="bisq", scale=TRUE)
yestbase = predict(zfit,seq(bb,aa,0.5))




# 5. ## CREATE A MATRIX OF CONTEMPORARY SSTS AT POINTS WHERE WE HA --------

reconpt <- rbind(EP_cores$coord,WP_cores$coord,c(0,204))

nreconpt <- dim(reconpt)[1] 
index <- matrix(seq(1:(nx*ny)),nx,ny)
reconsst <- c()
reconsst.a <- c()
zzpt <- c()
for (i in 1:nreconpt){
  lat <- reconpt[i,1]
  xlat <- which(ygrid == lat)
  
  lon <- reconpt[i,2]
  xlon <- which(xgrid == lon)
  
  zz <- index[xlon,xlat]
  zzpt <- c(zzpt,which(index1 == zz))
  zsst <- data[,which(index1 == zz)]
  zclim <- climdata[which(index1 == zz)]
  
  reconsst <- cbind(reconsst,zsst)			# SSTs in matrix form (160 x nreconpt)
  reconsst.a <- cbind(reconsst.a,zclim)		# SST climatology at each point (1 x nreconpt)
}


# 6. ## SCALE ALL THE PALEO DATA6. ## SCALE ALL THE PALEO DATA ------------

smoodat = c(sst.smooth$EP, sst.smooth$WP, list(cbind(seq(bb,aa,0.5),yestbase, yestbase)))
# sst.scale = scale.recon(smoodat, nreconpt, reconsst.mu, reconsst.sig)
sst.scale = scale.recon(smoodat, nreconpt, reconsst.a)
mat.unsc = sst.scale$unscaled			# Not scaled, just put into a matrix
# mat.sc = sst.scale$scaled				# Put into matrix and scaled using itself
mat.anom = sst.scale$anom				# Put into matrix and scaled using the climatology

mat.base = sst.scale$base
mat.baselim = sst.scale$baselim
# Pick which scaled series you want to use
mat = mat.baselim	
mat1 = mat
#recover the NAs
for(i in which(is.na(mat1))){
  mat1[i%%21,i%/%21+1] = apply(mat,1,mean,na.rm=TRUE)[i%%21]
}



# 7. ## GET PCS OF CONTEMPORARY and PROXY DATA AT PROXY LOCATIONS ---------

zs.p <- var(reconsst)
zsvd.p <- svd(zs.p)
pcs.p <- reconsst %*% zsvd.p$u
lambdas.p <- (zsvd.p$d/sum(zsvd.p$d))
npcs = 3
Z = pcs.p[,1:npcs] #predictors
##PC for paleo reconstruction
proxy_pcs.p <- mat1 %*% zsvd.p$u
proxy_Z = proxy_pcs.p[,1:npcs] #predictors
Tp=nrow(proxy_Z)



# # 8. # obtain the locations for each group ------------------------------

S = length(index2_i)
group_size = 30 #set the group size
number_groups = S/group_size

if ((S %% group_size)>0){
  message('Error: the number of fitting points (',S,") is not a multiple of the group size (",group_size,"). Modify the group size")
}else{
  set.seed(1)
  index2t=index2_i
  index_group=matrix(NA,group_size,number_groups)
  for (i in 1:number_groups) {
    index_group[,i]=sample(index2t,group_size,replace=FALSE)
    temp=which(index2t %in% index_group[,i])
    index2t=index2t[-temp]
    }
  index2=as.vector(index_group)#vector of the well distributed fitting points for each group
  
  S = length(index2)
  setwd(mainDir)
  list.save(index_group,paste("points_GroupSize_",group_size,".rds",sep=""))
}





# 9. #COMPUTE THE DISTANCE TO KNOTS POINTS ARRAYS -------------------------
#MATRIX OF DISTANCE FOR KNOTS TO EACH POINT OF THE FULL GRID (972 POINTS)
S_t=dim(data)[2]
dist_t = as.matrix(dist(grid1))
dist_knots_to_points_t = as.matrix(dist(rbind(knots, grid1)))[1:dim(knots)[1],(dim(knots)[1]+1):(dim(knots)[1]+S_t)]

#MATRIX OF DISTANCE FOR KNOTS TO EACH FITTING POINT (150 POINTS)

dist = as.matrix(dist(grid1[index2,]))
dist_knots_to_points = as.matrix(dist(rbind(knots, grid1[index2,])))[1:dim(knots)[1],(dim(knots)[1]+1):(dim(knots)[1]+S)]


# #10. PLOTS WITH THE FITTING POINTS ---------------------------------------


name="Fig01"
par(mfrow = c(1, 1), mar = c(1.5,1.4,0.55,0.5))
zfull = rep(NaN,(nx*ny))
zfull[index1] = apply(data.un,2,mean)
zmat = matrix(zfull,nrow=nx,ncol=ny)
zlim2=c(min(zfull[index1]),max(zfull[index1]))
image.plot(xgrid,ygrid,zmat,
           ylim=c(-10.3,10.3), xlim=range(min(xgrid),max(xgrid)),
           xlab="",ylab="", main="",
           cex.axis=1.1, col=rev(myPalette3(100)),
           zlim=zlim2,
           legend.lab = "SST (ºC)", axes = FALSE,smallplot=c(.915,0.935, 0.14,0.95),axis.args=list(tck = -0.5,mgp=c(0, 0.7, 0),col="black",lwd=1,cex.axis=1))
axis(2,at=c(-10,0,10),tck = -0.05,mgp=c(0, 0.5, 0))
axis(1,at=seq(100,300,50),tck = -0.05,mgp=c(0, 0.5, 0))
box()
maps::map("world2", 
          xlim=c(min(xgrid), max(xgrid)), ylim=c(min(ygrid)-5, max(ygrid)+5),
          plot = TRUE, add = TRUE
)


points(knots[,2],knots[,1],pch = 21, bg = "#41A59E", col = "black", cex = 1.5, lwd = 1.5)

cols2=c("#B1E4F9","#C1A119","#A050DC","#4CE08F","#909090")
for (i in 1:number_groups) {
  points(grid1[index_group[,i],2],grid1[index_group[,i],1],
         bg = cols2[i],col="black", cex = 1, lwd = 1,pch=22)
  
}


p1.base <- recordPlot()
invisible(dev.off())
setwd(dir_fig)

pdf(paste(name,".pdf",sep = ""), width = 10, height = 2)
p1.base
dev.off()

p1.base 


# #11. Settings  for Bayesian model in  stan ------------------------------

setwd(mainDir)

# sampler iterations
iterations = 100#3500
warmup = 70#1500
thin=(iterations-warmup)/2000
if (thin<=1){thin=1}
P = 2 #spatial predictors
no_nugget = FALSE
nchains = 3
n_beta_knots = dim(knots)[1]
max_treedepth=13
adapt_delta=0.99
rng_seed = 5000
init_rng_seed = 5000
set.seed(rng_seed)


model_name_ = 'normal_spatial_composite_beta'
model_file = "normal_spatial_composite_multivariate_3_PCs.stan"
# the real model name
model_name = sprintf('%s_%02dknots',model_name_,n_beta_knots)
message('This model is: ',model_name)


# # 12. initial guest of the mu regresion coefficients and variance -------

mle_loc=(1:S)*NA
mle_loc1=mle_loc
mle_loc2=mle_loc
mle_loc3=mle_loc
mle_scale = mle_loc

Z1=as.data.frame(Z)
colnames(Z1)=c("PC1","PC2","PC3")
for (i in 1:S) {
  mod<-gamlss(data[ind_y_cal,index2[i]]~PC1+PC2+PC3,family=NO(), data=Z1[ind_y_cal,], method=mixed(1,20))
  mle_loc[i]=mod$mu.coefficients[1]
  mle_scale[i]=exp(mod$sigma.coefficients[1])
  mle_loc1[i]=mod$mu.coefficients[2]
  mle_loc2[i]=mod$mu.coefficients[3]
  mle_loc3[i]=mod$mu.coefficients[4]
}



# # 13. create a list with the inputs for the Bayesian model in STAN --------

normalized_grid =sweep(grid1,2,c(min(grid1[,1]),min(grid1[,2])),FUN="-")
normalized_grid1 =sweep(normalized_grid,2,c(max(grid1[,1])-min(grid1[,1]),max(grid1[,2])-min(grid1[,2])),FUN="/")

normal_spatial_data = list(
  S = S, S_t=S_t, T = T_c, T_t=T_v,T_p=Tp, G = number_groups, 
  R = group_size, P = P, K = n_beta_knots, y = data[ind_y_cal,index2],
  dist = dist,dist_t=dist_t, dist_knots_to_points = t(dist_knots_to_points),dist_knots_to_points_t=t(dist_knots_to_points_t),  
  covars = as.matrix(normalized_grid1[index2,]),covars_t=as.matrix(normalized_grid1), covarsNS = Z[ind_y_cal,1:3],covarsNS_t = Z[ind_y_val,1:3],covarsNSp=proxy_Z,PCS = npcs)

initf_manual = function()list(
  loc = mle_loc,
  loc1 = mle_loc1,
  loc2 = mle_loc2,
  loc3 = mle_loc3,
  log_scale = log(mle_scale),
  log_nugget = log(rep(0.2,S)),
  
  basis_loc = matrix(0.01,P,n_beta_knots),
  basis_loc1 = matrix(0.01,P,n_beta_knots),
  basis_loc2 = matrix(0.01,P,n_beta_knots),
  basis_loc3 = matrix(0.1,P,n_beta_knots),
  basis_log_scale = matrix(0.01,P,n_beta_knots),
  
  
  intercept_loc = 0.1,
  intercept_loc1 = 0.1,
  intercept_loc2 = 0.1,
  intercept_loc3 = 0.1,
  intercept_log_scale = 0.1,
  
  sill_loc = 0.1,
  sill_loc1 = 0.1,
  sill_loc2 = 0.1,
  sill_loc3 = 0.1,
  sill_log_scale = 0.1,
  sill_log_nugget = 0.1,
  
  nugget_loc = 0.1,
  nugget_loc1 = 0.1,
  nugget_loc2 = 0.1,
  nugget_loc3 = 0.1,
  nugget_log_scale = 0.1,
  nugget_log_nugget = 0.1,
  
  range = max(dist)/4,
  range_loc = max(dist)/2,
  range_loc1 = max(dist)/2,
  range_loc2 = max(dist)/2,
  range_loc3 = max(dist)/2,
  range_log_scale = max(dist)/2,
  range_log_nugget = max(dist)/2,
  
  range_basis_loc = matrix(max(dist)/2,P,n_beta_knots),
  range_basis_loc1 = matrix(max(dist)/2,P,n_beta_knots),
  range_basis_loc2 = matrix(max(dist)/2,P,n_beta_knots),
  range_basis_loc3 = matrix(max(dist)/2,P,n_beta_knots),
  range_basis_log_scale = matrix(max(dist)/2,P,n_beta_knots)
)
gc()#to free up the memory


# # 14. #Run the STAN model -----------------------------------------------

ptm2 = proc.time()
stan_model = stan(model_file, data = normal_spatial_data, chains = 0)##to compile the C++ model 

stanfit = stan(fit = stan_model, data =normal_spatial_data,
               chains = nchains, iter=iterations, warmup=warmup,thin=thin,
               seed=rng_seed, refresh=1, init=initf_manual,
               control=list(max_treedepth=max_treedepth,adapt_delta=adapt_delta),cores=nchains)


#list.save(stanfit,"model_ff_predictors_18_knots_V7_2.rds")#UNCOMMENT TO SAVE THE RAW STAN RESULTS
sprintf("elapsed time for running STAN model: %s hours",round((proc.time() - ptm2)[3]/3600,2))
gc()#to free up the memory




# # SAVE POSTERIOR SAMPLES AND DISPLAY R HAT VALUES -----------------------

val_post=rstan::extract(stanfit)##extract the posterior into a list
###SAVE THE POSTERIOR SAMPLES
list.save(val_post,paste("predictors_",Nk,"_knots_.rds",sep=""))

message('R hat values for the fitting:')
print(quantile(summary(stanfit)$summary[,10],p=c(0.005,0.25,0.5,0.75,0.995),na.rm=TRUE))
