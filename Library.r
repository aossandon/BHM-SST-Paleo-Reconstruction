load_packages <- function(packages, quietly=TRUE){
  for(package in packages){
    installed = 
      if(quietly)
        suppressWarnings(suppressMessages(suppressPackageStartupMessages(
          require(package, character.only = TRUE, quietly = TRUE)
        )))
    else
      require(package, character.only = TRUE)
    
    if(!installed){
      install.packages(package,repos = "https://cloud.r-project.org")
      if(quietly)
        suppressWarnings(suppressMessages(suppressPackageStartupMessages(
          library(package, character.only = TRUE, quietly = TRUE)
        )))
      else
        library(package, character.only = TRUE)
    }
  }
  invisible(NULL)     
}

dir_create <- function(paths){
  # only creates the directory if it doesnt exist
  # defaults to recursive
  for(path in paths)
    if(!file.exists(path)) dir.create(path, recursive=TRUE)
  
  invisible(NULL)
}



# BINxygrid function ------------------------------------------------------

BINxygrid <- function(data, xgrid, ygrid, ntime, type, undefined, condition) {
  ny <- length(ygrid)
  nx <- length(xgrid)
  ngrid <- nx*ny
  xygrid <- matrix(0, nrow=(nx*ny), ncol=2)
  
  # Create the xygrid of lat and lon for each gridcell.
  i=0
  for (iy in 1:ny){
    for (ix in 1:nx){
      i=i+1
      xygrid[i,1] = ygrid[iy]
      xygrid[i,2] = xgrid[ix] 
    }
  }
  
  # Option for data that is defined in every single gridcell.
  if (type == "ngrid"){
    nsta = nx * ny
    data <- array(data=data, dim=c(nx,ny,ntime))
    BINdata <- matrix(NA, nrow=ntime, ncol=ngrid)
    
    for (i in 1:ntime){
      data1 <- data[,,i]
      BINdata[i,] <- data1
    }
    rm("data")
    
    # Return the results.
    results <- list(data = BINdata, grid = xygrid, nsta=nsta, xgrid=xgrid, ygrid=ygrid)
    return(results)
  }
  
  
  # Option for data that undefined in some gridcells (i.e. rainfall data just over India).
  if (type == "special"){
    data <- array(data=data, dim=c(nx,ny,ntime))
    
    index <- 1:(nx*ny)
    data1 <- data[,,1]
    
    
    if (condition == "greater"){
      index1 <- index[data1 > undefined]    
    }
    
    if (condition == "less"){
      index1 <- index[data1 < undefined]   
    }
    
    nsta <- length(index1)
    BINdata <- matrix(NA, nrow=ntime, ncol=nsta)
    for (i in 1:ntime){
      data1 <- data[,,i]
      data2 <- data1[index1]
      BINdata[i,] <- data2
    }
    
    # Return the results.
    results <- list(data = BINdata, grid = xygrid, ind=index1, nsta=nsta, xgrid=xgrid, ygrid=ygrid)
    return(results)    
  }
  
  
}


gcv.recon <- function(EP, WP, EPdeg, WPdeg, base){
  nEP = length(EP)			# number of east Pacific records
  nWP = length(WP)			# number of west Pacific records
  # ntot = nEP + nWP 
  
  EP.smooth = list()
  for (i in 1:nEP){
    N = dim(EP[[i]])[1]
    zfit = locfit(EP[[i]]$SST2 ~ EP[[i]]$age, deg=EPdeg, kern="bisq", scale=TRUE)
    
    if (round(EP[[i]]$age[N]) > 10){
      maxyr = 10
    } else {
      maxyr = round(EP[[i]]$age[N])
    }
    minyr = floor(EP[[i]]$age[1])
    
    if (base == "TRUE"){
      aa = 10
      bb = 0
    }
    yestbase = predict(zfit,seq(bb,aa,0.5))
    
    # print(yestbase)
    # print(minyr)
    # print(maxyr)
    
    xsmo = rep(NaN,21)
    xsmo[(2*minyr+1):(2*maxyr+1)] = yestbase[(2*minyr+1):(2*maxyr+1)]
    
    xEP.smooth = matrix(cbind(seq(bb,aa,0.5), yestbase, xsmo), 21, 3)
    
    #xEP.smooth = matrix(cbind(seq(minyr,maxyr,1),yestbase),length(yestbase),2)
    EP.smooth[[i]] = xEP.smooth
  }
  
  
  WP.smooth = list()
  for (i in 1:nWP){
    N = dim(WP[[i]])[1]
    zfit = locfit(WP[[i]]$SST2 ~ WP[[i]]$age, deg=WPdeg, kern="bisq", scale=TRUE)
    
    if (round(WP[[i]]$age[N]) > 10){
      maxyr = 10
    } else {
      maxyr = round(WP[[i]]$age[N])
    }
    minyr = floor(WP[[i]]$age[1])
    
    if (base == "TRUE"){
      aa = 10
      bb = 0
    }
    yestbase = predict(zfit,seq(bb,aa,0.5))
    
    xsmo = rep(NaN,21)
    xsmo[(2*minyr+1):(2*maxyr+1)] = yestbase[(2*minyr+1):(2*maxyr+1)]
    
    xWP.smooth = matrix(cbind(seq(bb,aa,0.5), yestbase, xsmo), 21, 3)
    
    #xWP.smooth = matrix(cbind(seq(minyr,maxyr,1),yest),length(yest),2)
    WP.smooth[[i]] = xWP.smooth
  }
  
  
  results = list(EP = EP.smooth, WP = WP.smooth)
  return(results)
  
}

scale.recon <- function(smoodat, sites, clim){
  sst.mat = matrix(NaN,length(seq(0:20)),(sites))
  sst.mat.a = sst.mat 		# Initialize matrix for anomalies (using climatology)
  sst.mat.base = sst.mat 		# Initialize matrix for using base 0 ka as zero
  sst.mat.baselim = sst.mat 		# Initialize matrix for using base 0 ka as zero
  
  for (i in 1:sites){
    xmu = clim[i]
    # print(i)
    xi = 2*smoodat[[i]][,1]
    # print(xi)
    #paleomu = mean(smoodat[[i]][,2], na.rm=TRUE)
    #paleosig = sd(smoodat[[i]][,2], na.rm=TRUE)
    sst.mat[xi+1,i] = smoodat[[i]][,3]
    sst.mat.a[xi+1,i] =  (smoodat[[i]][,3] - xmu )
    sst.mat.base[xi+1,i] = smoodat[[i]][,2] - smoodat[[i]][1,2]
    sst.mat.baselim[xi+1,i] =  smoodat[[i]][,3] - smoodat[[i]][1,2]
    
  }
  
  results = list(unscaled = sst.mat, anom = sst.mat.a, base = sst.mat.base, baselim=sst.mat.baselim)
  return(results)
  
}


geti <- function(x,i){
  d = length(dim(x))
  eval(parse(text=paste0('x[i',paste(rep(',',d-1),collapse=''),']')))
}


exp_cov_matrix <- function(x, sill, range, nugget){
  
  S = dim(x)[2]
  cov = matrix(NA,S,S);
  # off-diagonal elements
  for (i in 1:(S-1)) {
    for (j in (i+1):S) {
      # beta0 = partial sill = scale^2
      # beta1 = range = phi
      cov[i,j] = sill * exp(-1.0/range * x[i,j])
      cov[j,i] = cov[i,j]
    }
    cov[i,i] = sill + nugget
  }
  cov[S,S] = sill + nugget
  
  return(cov)
}





STT_quantile=function(sst_new){
  library(doParallel)
  ptm2 = proc.time()
  tmp=matrix(NA,dim(sst_new)[2],dim(sst_new)[3])
  sst_sta=list("pred_2_5"=tmp,"pred_25"=tmp,"pred_50"=tmp,"pred_75"=tmp,"pred_97_5"=tmp)
  nam=list.names(sst_sta)
  pob=c(0.025,0.25,0.5,0.75,0.975)
  numWorkers <- detectCores()-3
  registerDoParallel(numWorkers)
  res1=foreach (i=1:dim(sst_new)[3], .packages=c()) %dopar% {
    temp=matrix(NA,dim(sst_new)[2],length(pob))
    for (j in 1:length(nam)) {
      temp[,j]=apply(sst_new[,,i],2,quantile,probs=pob[j])
    }
    temp
  }
  
  temp2=array(unlist(res1), dim = c(dim(sst_new)[2],length(pob),dim(sst_new)[3]))
  stopImplicitCluster()
  
  for (j in 1:length(nam)) {
    sst_sta[[nam[j]]]=temp2[,j,]
  }
  sprintf("elapsed time: %s minutes",round((proc.time() - ptm2)[3]/60,2))
  return(sst_sta)
}




function_beta <- function(y, y_hat){
  
  beta = 1 -sum((y-y_hat)^2)/sum(y^2)
  return(beta)
  
}