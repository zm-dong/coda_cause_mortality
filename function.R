
# Functions to be used in the "master.R" file.  

# ------------ 2S-CoDa function ---------------#
CoDa.2step <- function(dx, ih,k,years,ages,ses,w.base,age.weight){
  # dx life table deaths stackked horizontial with dimensions time x age for each stackked matrix 
  # years = The years included in fitting period
  # ages = single ages    
  # ses = group names  
  n <- length(years)
  m <- length(ages)
  s <- length(ses)
  w.base1 <- w.base
  age.weight1 <- age.weight
  weight.fun <- function(x){
    return(w.base*(1-w.base)^x)    
  }
  rescale <- function(ws,sum){
    return(ws/sum)
  }
  
  x <- (1:nrow(dx))
  w <- sapply(X=x,FUN=weight.fun) 
  w.scaled <- sapply(X=w,FUN=rescale,sum=sum(w)) 
  
  #--------------- fit step 1 ---------------#
  fit.one <- CoDa_multi_para_wei_both(dx,nr.rank = 1,w.base=w.base1,age.weight.2=age.weight1)  
  R2.one <- cbind(fit.one$par$vs[1]^2/sum(fit.one$par$vs^2), fit.one$par$vs[2]^2/sum(fit.one$par$vs^2),fit.one$par$vs[3]^2/sum(fit.one$par$vs^2),fit.one$par$vs[4]^2/sum(fit.one$par$vs^2)) 

  # fit Arima model to kt
  order.one=c(1,1,1)
  fcst.one <- forecast(Arima(fit.one$kt, order=order.one ,include.drift=T, method="ML"), h=ih)
  kt.fit.one <- fcst.one$fitted
  kt.for.one <- fcst.one$mean
  ax.all <- matrix(fit.one$ax,length(ages),k)
  bx.all <- matrix(fit.one$bx,length(ages),k)
  #------------------- step 2 -------------------# 
  clr.proj.fit.one <- matrix(fit.one$kt,length(years),1) %*% t(fit.one$bx)
  clr.proj.for.one <- matrix(c(fit.one$kt, kt.for.one),length(years)+ih,1) %*% t(fit.one$bx) 
  res.1.in <- array(NA,dim=c(n,m,s))
  
  data.clr <- array(fit.one$clr.cent,dim=c(n,m,s)) 
  clr.proj.fit.one.ar <- array(clr.proj.fit.one,dim=c(n,m,s)) 
  
  for(i in 1:s){
    res.1.in[,,i] <-  data.clr[,,i] - clr.proj.fit.one.ar[,,i]  
  }
  clr.proj.for.in <- array(NA,dim=c((n+ih),m,s))
  clr.proj.fit.in <- array(NA,dim=c(n,m,s))
  kt.in.all <- matrix(NA,n+ih,s)
  bx.in.all <- matrix(NA,m,s)
  
  weight.in <- matrix(age.weight1,m,s)
  R2.in <- matrix(NA,s,4)
  for(i in 1:s){
    par.1 <- svd(res.1.in[,,i])
    R2.one <- cbind(par.1$d[1]^2/sum(par.1$d^2), par.1$d[2]^2/sum(par.1$d^2),par.1$d[3]^2/sum(par.1$d^2),par.1$d[4]^2/sum(par.1$d^2)) 
    R2.in[i,] <- R2.one
    U.in <- par.1$u
    V.in <- par.1$v
    S.in <- diag(par.1$d)
    
    bx.in <- V.in[,1]
    kt.in <- S.in[1,1]*U.in[,1]
    
    bx2.in <- V.in[,2]
    kt2.in <- S.in[2,2]*U.in[,2] 
    
    kt.in.for <- c(kt.in,forecast(auto.arima(kt.in,max.d=1),h=ih)$mean)
    kt2.in.for <- c(kt2.in,forecast(auto.arima(kt2.in,max.d=1),h=ih)$mean)
    
    clr.proj.for.in[,,i] <- matrix(kt.in.for,(n+ih),1) %*% t(bx.in) + matrix(kt2.in.for,(n+ih),1) %*% t(bx2.in) 
    clr.proj.fit.in[,,i] <- matrix(kt.in,n,1) %*% t(bx.in) + matrix(kt2.in,n,1) %*% t(bx2.in)
    kt.in.all[,i] <- kt.in.for
    bx.in.all[,i] <- bx.in
    
  }
  
  clr.proj.for.in.all <- clr.proj.for.in[,,1]
  clr.proj.fit.in.all <- clr.proj.fit.in[,,1]
  
  for(i in 2:s){
    clr.proj.for.in.all <- cbind(clr.proj.for.in.all,clr.proj.for.in[,,i])
    clr.proj.fit.in.all <- cbind(clr.proj.fit.in.all,clr.proj.fit.in[,,i])
  }
  
  #--- calculate forecast ---------------#
  clr.proj.for.one.all <- clr.proj.for.one + clr.proj.for.in.all
  clr.proj.fit.one.all <- clr.proj.fit.one + clr.proj.fit.in.all
  #projections
  
  #Inv clr
  BK.proj.fit.one <- clrInv(clr.proj.fit.one.all)
  BK.proj.for.one <- clrInv(clr.proj.for.one.all)
  
  #Add geometric mean
  proj.fit.one <- BK.proj.fit.one + fit.one$ax
  proj.for.one <- BK.proj.for.one + fit.one$ax
  proj.for.one.ar <- array((unclass(proj.for.one)*100000), dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  proj.fit.one.ar <- array((unclass(proj.fit.one)*100000), dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  dx.ar <- array(dx, dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  dx.ar1 <- array(acomp(dx), dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  #----- calculate explained varation -------------#
  VarJoint <- NULL
  VarInt <- NULL
  for(i in 1:k){VarJoint[i] = (norm(clr.proj.fit.one.ar[,,i],type='F')^2)/(norm(data.clr[,,i], type = 'F')^2)}
  for(i in 1:k){VarInt[i]  = (norm(clr.proj.fit.in[,,i],type='F')^2)/(norm(data.clr[,,i], type = 'F')^2)}
  Var.res <- 1 - VarJoint - VarInt
  exp.var <- cbind(VarJoint,VarInt,Var.res)
  
  #Life table setup
  a.one <- c(rep(2.5, length(ages))) 
  n.one <- c(rep(5,length(ages)))
  radix.one <- 100000
  
  for.mx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  for.dx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  fit.dx.one <- array(NA, dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  proj.for.one.tot <- proj.for.one.ar[,,1]
  for(i in 2:k){proj.for.one.tot <- proj.for.one.tot + proj.for.one.ar[,,i]}
  
  for.Lx.one <- LifeT.get.Lx((proj.for.one.tot), radix.one, a.one, n.one)$Lx
  
  for(i in 1:k){for.mx.one[,,i] <- proj.for.one.ar[,,i]/for.Lx.one}

  return(list(dx=dx,R2=R2.one,R2.in=R2.in,axs=ax.all,bxs=bx.all,kt=fit.one$kt,kt.for=kt.for.one,kt.in.all=kt.in.all,bx.in.all=bx.in.all, dx.forcast = proj.for.one.ar,BK.fit=clr.proj.fit.one.all,dx.fit = proj.fit.one.ar,dx.obs=dx.ar,exp.var=exp.var,for.mx.one=for.mx.one))
}

#-------------------------------------------------------#

CoDa_multi_para_wei_both <- function(dx,nr.rank,x,w.base,age.weight.2){
  
  weight.fun <- function(x){
    return(w.base*(1-w.base)^x)    
  }
  rescale <- function(ws,sum){
    return(ws/sum)
  }
  
  x <- (1:nrow(dx))
  w <- sapply(X=x,FUN=weight.fun) 
  w.scaled <- sapply(X=w,FUN=rescale,sum=sum(w)) 
  
  close.dx <- acomp(dx)
  ax <- geometricmeanCol(close.dx)
  dx.cent <-close.dx-ax
  clr.cent <- clr(dx.cent)
  
  # SVD: bx and kt
  par <- svd.triplet(clr.cent,ncp=nr.rank,row.w=rev(w.scaled),col.w=age.weight.2)
  U <- par$U
  V <- par$V
  S <- diag(par$vs)
  
  if(nr.rank == 1){
    bx<- V[,1]
    kt<- S[1,1]*U[,1]
    return(list(ax=ax,bx=bx,kt=kt,dx=dx,par=par,clr.cent=clr.cent))
  }
  else if(nr.rank == 2){
    
    bx<- V[,1]
    kt<- S[1,1]*U[,1] 
    
    bx2 <- V[,2]
    kt2 <- S[2,2]*U[,2] 
    return(list(ax=ax,bx=bx,kt=kt,bx2=bx2,kt2=kt2,dx=dx,par=par))
  }
  else{
    bx<- V[,1]
    kt<- S[1,1]*U[,1] 
    
    bx2 <- V[,2]
    kt2 <- S[2,2]*U[,2]  
    
    bx3 <- V[,3]
    kt3 <- S[3,3]*U[,3]
    return(list(ax=ax,bx=bx,kt=kt,bx2=bx2,kt2=kt2,bx3=bx3,kt3=kt3,dx=dx,par=par))  
  }
}


#-------------------------------------------------------#

LifeT.get.Lx<-function(dx, radix, a, n){
  # t times x 
  dx<-dx
  lx<-matrix(radix, nrow(dx), ncol(dx))
  Lx<-matrix(NA, nrow(dx), ncol(dx))
  Tx<-matrix(NA, nrow(dx), ncol(dx))
  for(j in 1:nrow(dx)){
    
    lx[j,2:ncol(dx)]<- radix-cumsum(dx[j,1:(ncol(dx)-1)])
    
    Lx[j,]<- (lx[j,]*n) - (dx[j,]*(n-a))
    
    for(i in 1:ncol(dx)){
      Tx[j,i]<-sum(Lx[j, i:ncol(Lx)], na.rm=T)
    }
    Tx[Tx==0]<-NA
  }
  
  
  colnames(lx)<-colnames(Lx)<-colnames(dx)
  rownames(lx)<-rownames(Lx)<-rownames(dx)
  return(list(dx=dx, lx=lx, Lx=Lx))
}

#-------------------------------------------------------#

#CT-CoDa function 
CoDa.CT <- function(dx, ih,k,years,ages,ses){
  # dx life table deaths stacked horizontal with dimensions time x age for each stacked matrix 
  # years = The years included in fitting period
  # ages = single ages    
  # ses = group names  
  n <- length(years)
  m <- length(ages)
  fit.one <- CoDa_multi_para(dx,nr.rank = 3)  
  R2.one <- cbind(fit.one$par$d[1]^2/sum(fit.one$par$d^2), fit.one$par$d[2]^2/sum(fit.one$par$d^2),fit.one$par$d[3]^2/sum(fit.one$par$d^2),fit.one$par$d[4]^2/sum(fit.one$par$d^2)) 
  # fit Arima model to kt
  order.one=c(1,1,1)
  fcst.one <- forecast(auto.arima(fit.one$kt,max.d=1),h=ih)
  kt.fit.one <- fcst.one$fitted
  kt.for.one <- fcst.one$mean
  
  order.two=c(1,0,1)
  fcst.two <- forecast(auto.arima(fit.one$kt2,max.d=1),h=ih)
  kt.fit.two <- fcst.two$fitted
  kt.for.two <- fcst.two$mean

  order.3=c(1,0,1)
  fcst.3 <- forecast(auto.arima(fit.one$kt3,max.d=1),h=ih)
  kt.fit.3 <- fcst.3$fitted
  kt.for.3 <- fcst.3$mean
  
  ax.all <- matrix(fit.one$ax,length(ages),k)
  bx.all <- matrix(fit.one$bx,length(ages),k)
  bx2.all <- matrix(fit.one$bx2,length(ages),k)
  
  #projections
  clr.proj.fit.one<- matrix(fit.one$kt,length(years),1) %*% t(fit.one$bx) + matrix(fit.one$kt2,length(years),1) %*% t(fit.one$bx2) + matrix(fit.one$kt3,length(years),1) %*% t(fit.one$bx3)
  clr.proj.for.one <- matrix(c(fit.one$kt, kt.for.one),length(years)+ih,1) %*% t(fit.one$bx) + matrix(c(fit.one$kt2, kt.for.two),length(years)+ih,1) %*% t(fit.one$bx2) + matrix(c(fit.one$kt3, kt.for.3),length(years)+ih,1) %*% t(fit.one$bx3)
  #clr.proj.for.one <- matrix(c(fit.one$kt, kt.for.one),length(years)+ih,1) %*% t(fit.one$bx) + matrix(c(fit.one$kt2, kt.for.two),length(years)+ih,1) %*% t(fit.one$bx2) 
  #Inv clr
  BK.proj.fit.one <- clrInv(clr.proj.fit.one)
  BK.proj.for.one <- clrInv(clr.proj.for.one)
  
  #Add geometric mean
  proj.fit.one <- BK.proj.fit.one + fit.one$ax
  proj.for.one <- BK.proj.for.one + fit.one$ax
  proj.for.one.ar <- array((unclass(proj.for.one)*100000), dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  #proj.for.one.ar <- array((unclass(proj.for.one)*1), dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  
  proj.fit.one.ar <- array((unclass(proj.fit.one)*100000), dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  dx.ar <- array(dx, dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  dx.ar1 <- array(acomp(dx), dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  #Life table setup
  a.one <- c(rep(2.5, length(ages))) 
  n.one <- c(rep(5,length(ages)))
  radix.one <- 100000
  
  for.mx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  for.dx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  fit.dx.one <- array(NA, dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  proj.for.one.tot <- proj.for.one.ar[,,1]
  for(i in 2:k){proj.for.one.tot <- proj.for.one.tot + proj.for.one.ar[,,i]}
  
  for.Lx.one <- LifeT.get.Lx((proj.for.one.tot), radix.one, a.one, n.one)$Lx
  
  for(i in 1:k){
    for.mx.one[,,i] <- proj.for.one.ar[,,i]/for.Lx.one
  }
  return(list(dx=dx,R2=R2.one,axs=ax.all,bxs=bx.all,kt=fit.one$kt,kt.for=kt.for.one,kt2=fit.one$kt2,kt.for2=kt.for.two, dx.forcast = proj.for.one.ar,bx2.all=bx2.all,dx.fit=proj.fit.one.ar,dx.obs=dx.ar,sv=fit.one$par$d,BK.obs=fit.one$clr.cent,BK.fit=clr.proj.fit.one,for.mx.one=for.mx.one,proj.for.one.tot=proj.for.one.tot))
}

#-------------------------------------------------------#
CoDa_multi_para <- function(dx,nr.rank){
  close.dx <- acomp(dx)
  ax <- geometricmeanCol(close.dx)
  dx.cent <-close.dx-ax
  clr.cent <- clr(dx.cent)
  
  # SVD: bx and kt
  par<- svd(clr.cent, nu=nr.rank, nv=nr.rank)
  U <- par$u
  V <- par$v
  S <- diag(par$d)
  
  if(nr.rank == 1){
    
    bx<- V[,1]
    kt<- S[1,1]*U[,1]
    
    return(list(ax=ax,bx=bx,kt=kt,dx=dx,par=par))
  }
  #else{
  if(nr.rank == 2){
    
    bx<- V[,1]
    kt<- S[1,1]*U[,1] 
    
    bx2 <- V[,2]
    kt2 <- S[2,2]*U[,2] 
    
    return(list(ax=ax,bx=bx,kt=kt,bx2=bx2,kt2=kt2,dx=dx,par=par)) 
  }
  else{
    bx<- V[,1]
    kt<- S[1,1]*U[,1] 
    bx2 <- V[,2]
    kt2 <- S[2,2]*U[,2]
    bx3 <- V[,3]
    kt3 <- S[3,3]*U[,3]
  }
  return(list(ax=ax,bx=bx,kt=kt,bx2=bx2,kt2=kt2,bx3=bx3,kt3=kt3,dx=dx,par=par,clr.cent=clr.cent)) 
}

#----------------------------------------------------#
#VECM-CoDA function 
CoDa.vecm <- function(dx, ih,k,years,ages,ses,ir,kt.for){
  # dx life table deaths stacked horizontal with dimensions time x age for each stackked matrix 
  # years = The years included in fitting period
  # ages = single ages    
  # ses = group names  
  n <- length(years)
  m <- length(ages)
  s <- length(ses)
  
  close.dx <- acomp(dx)
  ax <- geometricmeanCol(close.dx)
  dx.cent <-close.dx-ax
  clr.cent <- clr(dx.cent)
  clr.cent.ar.2 <- array(clr.cent,dim=c(n,m,k))
  
  clr.proj.for.one1 <- array(NA,dim=c((n+ih),m,length(ses)))
  bx.all <- matrix(NA,m,s)
  kt.all <- matrix(NA,n,s)
  bx.all.2 <- matrix(NA,m,s)
  kt.all.2 <- matrix(NA,n,s)
  R2 <- matrix(NA,s,3)
  
  for(i in 1:k){
    par<- svd(clr.cent.ar.2[,,i], nu=2, nv=2)
    U <- par$u
    V <- par$v
    S <- diag(par$d)
    
    bx<- V[,1]
    kt<- S[1,1]*U[,1]
    
    bx2 <- V[,2]
    kt2 <- S[2,2]*U[,2]
    
    kt.all[,i] <- kt 
    bx.all[,i] <- bx
    
    kt.all.2[,i] <- kt2 
    bx.all.2[,i] <- bx2
    R2[i,] <- cbind(par$d[1]^2/sum(par$d^2), par$d[2]^2/sum(par$d^2),par$d[3]^2/sum(par$d^2))
  }  
  
  #R2.one <- cbind(fit.one$par$d[1]^2/sum(fit.one$par$d^2), fit.one$par$d[2]^2/sum(fit.one$par$d^2),fit.one$par$d[3]^2/sum(fit.one$par$d^2),fit.one$par$d[4]^2/sum(fit.one$par$d^2)) 
  ax.all <- matrix(ax,length(ages),k)
  colnames(kt.all) <- as.character(1:s)
  colnames(kt.all.2) <- as.character(1:s)
  
  jo.test <- ca.jo(kt.all,type="trace",ecdet="trend",K=2,spec="transitory")
  print(summary(jo.test))
  
  #projections
  vecvar <- vec2var(jo.test, r = ir) 
  for.kt.all <- predict(vecvar,n.ahead=ih)
  
  clr.proj.for.one <- array(NA,dim=c((n+ih),m,s))
  clr.proj.fit.one <- array(NA,dim=c((n),m,s))
  kt.all.for <- matrix(NA,(n+ih),s)
  kt.all.for.2 <- matrix(NA,(n+ih),s)
  
  for(i in 1:s){
    kt.all.for[,i] <- c(kt.all[,i],for.kt.all$fcst[[i]][,1])
    kt.all.for.2[,i] <- c(kt.all.2[,i],forecast(auto.arima(kt.all.2[,i],max.d=0),h=ih)$mean)
    
    clr.proj.for.one[,,i] <- matrix(kt.all.for[,i],(n+ih),1) %*% t(bx.all[,i]) + matrix(kt.all.for.2[,i],(n+ih),1) %*% t(bx.all.2[,i])
    clr.proj.fit.one[,,i] <- matrix(kt.all[,i],n,1) %*% t(bx.all[,i]) + matrix(kt.all.2[,i],n,1) %*% t(bx.all.2[,i])
  }
  
  clr.proj.for.one.1 <- clr.proj.for.one[,,1]
  clr.proj.fit.one.1 <- clr.proj.fit.one[,,1]
  
  for(i in 2:k){
    clr.proj.for.one.1 <- cbind(clr.proj.for.one.1,clr.proj.for.one[,,i])
    clr.proj.fit.one.1 <- cbind(clr.proj.fit.one.1,clr.proj.fit.one[,,i])
  }
  #Inv clr
  BK.proj.for.one <- clrInv(clr.proj.for.one.1)
  BK.proj.fit.one <- clrInv(clr.proj.fit.one.1)
  
  
  #Add geometric mean
  proj.for.one <- BK.proj.for.one + ax
  proj.fit.one <- BK.proj.fit.one + ax  
  
  proj.for.one.ar <- array((unclass(proj.for.one)*100000), dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  proj.fit.one.ar <- array((unclass(proj.fit.one)*100000), dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))

  dx.ar <- array(dx, dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  dx.ar1 <- array(acomp(dx), dim=c(n,m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  #Life table setup
  a.one <- c(rep(2.5, length(ages))) 
  n.one <- c(rep(5,length(ages)))
  radix.one <- 100000
  
  for.mx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  for.dx.one <- array(NA, dim=c((n+ih),m,k),dimnames=list(as.character(years[1]:(years[length(years)]+ih)),as.character(ages),ses))
  fit.dx.one <- array(NA, dim=c((n),m,k),dimnames=list(as.character(years[1]:(years[length(years)])),as.character(ages),ses))
  
  proj.for.one.tot <- proj.for.one.ar[,,1]
  for(i in 2:k){proj.for.one.tot <- proj.for.one.tot + proj.for.one.ar[,,i]}
  
  for.Lx.one <- LifeT.get.Lx((proj.for.one.tot), radix.one, a.one, n.one)$Lx
  
  for(i in 1:k){for.mx.one[,,i] <- proj.for.one.ar[,,i]/for.Lx.one}
  return(list(dx=dx,axs=ax.all,bxs=bx.all,kt=kt.all,kt.all.for=kt.all.for,R2=R2,kt.all.for.2=kt.all.for.2,bx.all.2=bx.all.2,dx.fit=proj.fit.one.ar,BK.fit=clr.proj.fit.one ,dx.forcast = proj.for.one.ar,dx.obs=dx.ar,jo.test=jo.test,for.mx.one=for.mx.one))
}


#--------------------------------------------------------#
# LC function 


# Lee-Carter model
# mx =  matrix of mx, by year (rows) and age (columns)
# t= number of years to be forecasted

COD.LC<- function(mx, t,k,ages,years){
  
  library(forecast)
  n <- length(years)
  m <- length(ages)
  d.i.for.lc <- q.i <- LC.mx.for.ar <- array(NA,dim=c((n+t),m,k))
  kt.all <- matrix(NA,(n+t),k)
  bx.all <- matrix(NA,m,k)
  R2 <- matrix(NA,1,k)
  for(i in 1:k){
    lc.temp <- LC.no.jump(mx[,,i], t=t)
    LC.mx.for.ar[,,i] <- lc.temp$forefit.mx
    kt.all[,i] <- lc.temp$kt
    bx.all[,i] <- lc.temp$bx
    R2[,i] <- lc.temp$R2
  }
  
  LC.mx.for.all <- apply(LC.mx.for.ar[,,],MARGIN=c(1,2),FUN=sum)
  
  for(i in 1:k){q.i[,,i] <-  (LC.mx.for.ar[,,i] * 5)/(1+(5-2.5)*LC.mx.for.all)}
  
  lx.all <- Lx.all <- qi.all <- matrix(NA,(n+t),m)
  
  for(p in 1:(n+t)){
    qi.all[p,] <-  Lifetablebrige.mx(LC.mx.for.all[p,] ,age=ages,nax=2.5)$qx
    Lx.all[p,] <-  Lifetablebrige.mx(LC.mx.for.all[p,] ,age=ages,nax=2.5)$Lx
  }
  
  for(i in 1:k){
    d.i.for.lc[,,i] <- LC.mx.for.ar[,,i] * Lx.all 
  }
  
  return(list(dx.forcast=d.i.for.lc,kt.all=kt.all,bx.all=bx.all,for.mx.one=LC.mx.for.ar,R2=R2))   
}

#--------------------------------------------------------#


# Lee-Carter model
# mx =  matrix of mx, by year (rows) and age (columns)
# t= number of year to be forecasted

LC.no.jump<- function(mx, t){
  
  library(forecast)
  
  ln.mx<- log(mx)
  ax<- colMeans(ln.mx, na.rm=T)
  lnmx.cent<-sweep(ln.mx, 2, ax, FUN="-")
  
  # SVD: bx and kt
  par<- svd(lnmx.cent)
  U <- par$u
  V <- par$v
  S <- diag(par$d)
  R2 <- par$d[1]^2/sum(par$d^2)
  bx<- V[,1]/ sum(V[,1])
  kt<- S[1,1]*U[,1]* sum(V[,1])
  
  # Pick an ARIMA model (order=c())
  kt.fit<-forecast(Arima(kt, order=c(0,1,0), include.drift=T), h=t)$fitted
  kt.for<-forecast(Arima(kt, order=c(0,1,0), include.drift=T), h=t)$mean
  
  variability<- cumsum((par$d)^2/sum((par$d[1:length(par$d)])^2))
  
  #projections
  logmx.proj<- matrix(c(kt.fit, kt.for),(nrow(mx)+t),1) %*% t(bx)
  lnmx.proj2<- sweep(logmx.proj, 2, ax, FUN="+")
  proj<- exp(lnmx.proj2)
  
  #jump-off
  #pert<- mx[nrow(mx),]-proj[nrow(mx),]
  
  #forecast.mx<- proj[nrow(mx):nrow(proj),]
  #for(i in 1:nrow(forecast.mx)){
  #for(j in 1:ncol(forecast.mx)){
  #    forecast.mx[i,j]<- proj[(i+(nrow(mx)-1)),j]+pert[j]
  #  }
  #}
  #forefit.mx<- rbind(proj[1:(nrow(mx)-1),], forecast.mx)
  output<-list(mx,logmx.proj,bx=bx,variability,kt=c(kt,kt.for),forefit.mx= proj,R2=R2)
  #names(output)<- c("mx", "BxKt", "VarExp","kt","for.mx")
  return(output)
}


#------------------------------------------------------------#
Lifetablebrige.mx=function(mx,age,nax){
  # vector with mx
  # age is start age of interval (age = c(0,1,seq(5,110,5))) 
  # only for data with ages higher than zero   
  n <- c(diff(age), 9999)  
  qx <- (n * mx)/(1 + (n - nax) * mx) 
  qx <- c(qx[-(length(qx))], 1)
  qx<- ifelse(qx>1,0.998,qx)
  nage <- length(age)
  npx <- 1 - qx
  l0 = 100000
  lx <- round(cumprod(c(l0, npx)))
  ndx <- -diff(lx)
  lxpn <- lx[-1]
  nLx <- n * lxpn + ndx * nax
  #nLx <- c(nLx[-(length(qx))], lx[length(qx)]*nax[length(qx)])
  nLx <- c(nLx[-(length(qx))], lx[length(qx)]/mx[length(qx)])
  
  Tx <- c(rev(cumsum(rev(nLx[-length(nLx)]))),0)
  Tx <- c(Tx [-(length(qx))], nLx[length(qx)])
  lx <- lx[1:length(age)]
  ex <- Tx/lx
  list(mx=mx,nax=nax,qx=qx,lx=lx,ndx=ndx,Lx=nLx,Tx=Tx,ex=ex)  
}


# Copying over the RobRSVD Function

RobRSVD <-
  function(data, irobust = F, huberk = 1.345, 
           iinituv = F, inits, initu, initv, 
           niter = 1000, tol = 1e-5, istablize = T, 
           uspar = 0, vspar = 0, iugcv = F, ivgcv = F,
           usparmax = 1e4, usparmin = 1e-10, nuspar = 14,iloguspar = T, 
           vsparmax = 1e4, vsparmin = 1e-10, nvspar = 14, ilogvspar = T){
    
    # Robust Functional Singular Value Decomposition
    
    ugcvmat = as.numeric()
    vgcvmat = as.numeric()
    
    # Check if the paramenters are appropriate. #
    
    size_data = c(dim(data))
    m = size_data[1]
    n = size_data[2]
    if (istablize){
      # The original data is centered by the median of the whole panel/matrix
      myscale = 1.4785*median(abs(c(data) - median(data)))
      localdata = data / myscale
    }else{
      myscale = 1
      localdata = data
    }
    
    if (!iinituv){
      # Classical SVD is applied on the localdata
      # Only the first u vector and v vector are retained 
      svdold = svd(localdata, 1, 1)
      uold = svdold$u
      vold = svdold$v
      sold = svdold$d[1]
    }else{
      uold = initu
      sold = inits
      vold = initv
      if (istablize)
        sold = sold / myscale
    }
    uold = sold * uold
    
    # Use the initial estimate of u and v vectors to reconstruct the data matrix
    Appold = uold %*% t(vold) 
    # Residuals (matrix) after the initial estimate
    Rmat = localdata - Appold
    Rvec = c(Rmat) # stack the Rmat columns
    mysigma = median(abs(Rvec))/0.675
    iter = 1
    localdiff = 9999
    diffvec = as.numeric()
    uspar_current = uspar
    vspar_current = vspar
    
    ugcvscore = as.numeric()
    vgcvscore = as.numeric()
    
    while (localdiff>tol & iter<niter){
      if (irobust){
        Wmat = huberWeightLS(Rmat/mysigma, huberk)
      }else{
        Wmat = matrix(1, m, n)
      }
      # update u first 
      if (!iugcv){
        uterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
          (2*mysigma^2) * 
          (c(t(vold) %*% (diag(n) + vspar*ssmatls(n)$y) %*% vold) * (diag(m) + uspar*ssmatls(m)$y) - diag(sum(vold^2), m))
        uterm2 = (Wmat*localdata) %*% vold
        unew = solve(uterm1) %*% uterm2
      }else{ #choosing the penality parameter by cross-validation according to GCV criteria.
        if (nuspar<0)
          stop("number of smoothing parameter can not be negative")
        else{
          if (iloguspar){
            usparvec = 10^seq(log10(usparmin), log10(usparmax), length.out=nuspar)
          }else{
            usparvec = seq(usparmin, usparmax, length.out=nuspar)
          }
          ugcvvec = as.numeric()
          ugcvmat = as.numeric()
          for (iter_uspar in 1:nuspar){
            u_nsrobust = solve(diag(colSums(diag(c(vold^2)) %*% t(Wmat)))) %*% (Wmat*localdata) %*% vold
            usterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
              (2*mysigma^2) *
              (c(t(vold) %*% (diag(n) + vspar_current*ssmatls(n)$y) %*% vold) * (diag(m)+usparvec[iter_uspar]*ssmatls(m)$y)-diag(sum(vold^2), m))
            usterm2 = (Wmat*localdata) %*% vold
            u_srobust = solve(usterm1) %*% usterm2
            smooth_u = solve(usterm1) %*% diag(colSums(diag(c(vold^2)) %*% t(Wmat)))
            gcv_ut = m * sum((u_nsrobust-u_srobust)^2) / (m - sum(diag(smooth_u)))^2
            ugcvvec = c(ugcvvec, gcv_ut)
            ugcvmat = cbind(ugcvmat, u_srobust/(sqrt(sum(u_srobust^2))))
          }
          uspar_current = usparvec[which.min(ugcvvec)]
          ugcvscore = cbind(usparvec, ugcvvec)
          
          uterm1 = diag(colSums(diag(c(vold^2)) %*% t(Wmat))) + 
            (2*mysigma^2) *
            (c(t(vold) %*% (diag(n)+vspar_current*ssmatls(n)$y) %*% vold) * (diag(m)+uspar_current*ssmatls(m)$y) - diag(sum(vold^2), m))
          uterm2 = (Wmat*localdata) %*% vold
          unew = solve(uterm1) %*% uterm2 
        }#nuspar else
      }
      
      # update v secondly
      if (!ivgcv){
        vterm1 =  diag(colSums(diag(c(unew^2))%*%Wmat)) + 
          (2*mysigma^2) * 
          (c(t(unew)%*%(diag(m)+uspar*ssmatls(m)$y)%*%unew) * (diag(n)+vspar*ssmatls(n)$y) - diag(sum(unew^2), n))
        vterm2 = t(Wmat*localdata) %*% unew
        vnew = solve(vterm1) %*% vterm2
      }else{
        if (nvspar<0)
          stop("number of smoothing parameter can not be negative")
        else{
          if (ilogvspar){
            vsparvec = 10^seq(log10(vsparmin), log10(vsparmax), length.out=nvspar)
          }else{
            vsparvec = seq(vsparmin, vsparmax, length.out=nvspar)
          }
          vgcvvec = as.numeric()
          vgcvmat = as.numeric()
          for (iter_vspar in 1:nvspar){
            v_nsrobust = solve(diag(colSums(diag(c(unew^2))%*%Wmat))) %*% t(Wmat*localdata) %*% unew
            vsterm1 = diag(colSums(diag(c(unew^2))%*%Wmat)) + 
              (2*mysigma^2) * 
              (c(t(unew)%*%(diag(m)+uspar_current*ssmatls(m)$y)%*%unew) * (diag(n)+vsparvec[iter_vspar]*ssmatls(n)$y) - diag(sum(unew^2), n))
            vsterm2 = t(Wmat*localdata) %*% unew
            v_srobust = solve(vsterm1) %*% vsterm2
            smooth_v = solve(vsterm1) %*% diag(colSums(diag(c(unew^2))%*%Wmat))
            gcv_vt = n*sum((v_nsrobust-v_srobust)^2) / (n-sum(diag(smooth_v)))^2
            vgcvvec = c(vgcvvec, gcv_vt)
            vgcvmat = cbind(vgcvmat, v_srobust/sqrt(sum(v_srobust^2)))
          }
          vspar_current = vsparvec[which.min(vgcvvec)]
          
          vgcvscore = cbind(vsparvec, vgcvvec)
          
          vterm1 = diag(colSums(diag(c(unew^2))%*%Wmat)) + 
            (2*mysigma^2) * 
            (c(t(unew)%*%(diag(m)+uspar_current*ssmatls(m)$y)%*%unew) * (diag(n)+vspar_current*ssmatls(n)$y) - diag(sum(unew^2), n))
          vterm2 = t(Wmat*localdata) %*% unew
          vnew = solve(vterm1) %*% vterm2
        }# nvspar>=0
      }# ivgcv=T
      
      Appnew = unew %*% t(vnew)
      Rmat = localdata - Appnew
      
      localdiff = max(abs(Appnew-Appold))
      Appold = Appnew
      uold = sqrt(sum(vnew^2)) * unew
      vold = vnew / sqrt(sum(vnew^2)) # v is normalized to length 1
      
      iter = iter + 1
      diffvec = c(diffvec, localdiff)
    }# while
    
    v = vold
    s = myscale * sqrt(sum(uold^2))
    u = uold / sqrt(sum(uold^2))
    
    if (iugcv)
      uspar = uspar_current
    if (ivgcv)
      vspar = vspar_current
    
    diagout = list(ugcvscore=ugcvscore, vgcvscore=vgcvscore, ugcvmat=ugcvmat, vgcvmat=vgcvmat)
    return(list(s=s, v=v, u=u, diagout=diagout))
  }

