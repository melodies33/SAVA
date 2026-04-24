library(truncnorm)
library(tidyverse)
library(foreach)
library(doParallel)
library(Rcpp)
library(onlineFDR)
library(SAVA)
library(BSDA)

#### Counterexample 1 ####

n.time = 100 # number of total times
q = 0.1 # FDR level
n.rep = 1000 # repeated number
w0 = q
k = 25
{
  mu.v = rep(0, n.time)
  mu.v[1:20] = 2.5
  mu.v[21:50] = 0.01
  mu.v[51:60] = 2.5
  mu.v[61:70] = 0.001
  mu.v[71:80] = 2.5
  mu.v[81:n.time] = 0.01
  FSR = rep(NA, n.time)
  TSR = rep(NA, n.time)
  mFSR = rep(NA, n.time)
  cl <- makeCluster(3)
  registerDoParallel(cl)
  set.seed(667)
  ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_counterexample1(n.time, mu.v, q)
                    fsrvec = re$FSP
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    fsnumber = re$fsnumber
                    decidenumber = re$decidenumber
                    c(fsrvec, tprtvec, fsnumber, decidenumber)
                  }  
  stopImplicitCluster()
  stopCluster(cl)
  result = rowMeans(ree)
  FSR = result[1:n.time]
  TSR = result[(n.time + 1):(2*n.time)]
  mFSR = result[(2*n.time+1):(3*n.time)]/result[(3*n.time+1):(4*n.time)]
  
  FSR.sava = rep(NA, n.time)
  TSR.sava = rep(NA, n.time)
  mFSR.sava = rep(NA, n.time)
  cl <- makeCluster(3)
  registerDoParallel(cl)
  set.seed(66)
  ree = foreach(m = 1:n.rep, .combine = cbind,
                .packages = c('truncnorm','SAVA')) %dopar%{
                  re = sava_ctexam1_compare(n.time, mu.v, q, w0, k)
                  fsrvec = re$FSP
                  tprtvec = re$TP / sapply(re$realsig, max, 1)
                  fsnumber = re$fsnumber
                  decidenumber = re$decidenumber
                  c(fsrvec, tprtvec, fsnumber, decidenumber)
                }  
  stopImplicitCluster()
  stopCluster(cl)
  result = rowMeans(ree)
  FSR.sava = result[1:n.time]
  TSR.sava = result[(n.time + 1):(2*n.time)]
  mFSR.sava = result[(2*n.time+1):(3*n.time)]/result[(3*n.time+1):(4*n.time)]

  datare = tibble(time = rep(1:n.time, 2),
                  FSR = c(FSR, FSR.sava),
                  TSR = c(TSR, TSR.sava),
                  mFSR = c(mFSR, mFSR.sava),
                  method = rep(c('Method 1', 'SAVA'), each = n.time))
  
  write_csv(datare, 'sava_counterexample_1.csv')
  
}

#### Counterexample 2 ####
n.time = 100 # number of total times
prob.start = 1# probability of another hypothesis starting to be observed
q = 0.1 # FDR level
n.rep = 1000
plrtio = 0.5
w0 = q
k = 25
{
  mu = 1
  FSR = rep(NA, n.time)
  TSR = rep(NA, n.time)
  mFSR = rep(NA, n.time)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(667)
  ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    
                    re = sava_counterexample2(n.time, prob.start, plrtio, q, mu, w0, k)
                    fsrvec = re$FSP
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    fsnumber = re$fsnumber
                    decidenumber = re$decidenumber
                    c(fsrvec, tprtvec, fsnumber, decidenumber)
                  }
  
  stopImplicitCluster()
  stopCluster(cl)
  result = rowMeans(ree)
  FSR = result[1:n.time]
  TSR = result[(n.time + 1):(2*n.time)]
  mFSR = result[(2*n.time + 1):(3*n.time)]/result[(3*n.time + 1):(4*n.time)]
  
  FSR.sava = rep(NA, n.time)
  TSR.sava = rep(NA, n.time)
  mFSR.sava = rep(NA, n.time)
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(667)
  ree = foreach(m = 1:n.rep, .combine = cbind,
                .packages = c('truncnorm','SAVA')) %dopar%{
                  re = sava_ctexam2_compare(n.time, prob.start, plrtio, q, mu, w0, k)
                  fsrvec = re$FSP
                  tprtvec = re$TP / sapply(re$realsig, max, 1)
                  fsnumber = re$fsnumber
                  decidenumber = re$decidenumber
                  c(fsrvec, tprtvec, fsnumber, decidenumber)
                }  
  stopImplicitCluster()
  stopCluster(cl)
  result = rowMeans(ree)
  FSR.sava = result[1:n.time]
  TSR.sava = result[(n.time + 1):(2*n.time)]
  mFSR.sava = result[(2*n.time + 1):(3*n.time)]/result[(3*n.time + 1):(4*n.time)]
  
  datare = tibble(time = rep(1:n.time, 2),
                  FSR = c(FSR, FSR.sava),
                  TSR = c(TSR, TSR.sava),
                  mFSR = c(mFSR, mFSR.sava),
                  method = rep(c('Method 2', 'SAVA'), each = n.time))
  
  write_csv(datare, 'sava_counterexample_2.csv')
  
}

##### Performances under different k in Gaussian case ######
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # FDR level
kvec = seq(2,100,2)
n.rep = 1000# repeated number
gamma = 0
{
  mu = 0.1
  plrtio = 0.5
  w0 = q/2
  FSRmat = matrix(NA, n.time, length(kvec))
  TSRmat = matrix(NA, n.time, length(kvec))
  FSPhatmat = matrix(NA, n.time, length(kvec))
  cl <- makeCluster(9)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(kvec)){
    k = kvec[i]
    ree = foreach(j = 1:n.rep, .combine = cbind,
                  .packages = c('dplyr','SAVA')) %dopar%{
                    re = sava_gaussian(n.time, prob.start, plrtio, mu, q, gamma, k, w0)
                    fsrvec = re$FSP
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    fsphat = re$FSPhat
                    c(fsrvec, tprtvec, fsphat)
                  }
    result = rowMeans(ree)
    FSRmat[,i] = result[1:n.time]
    TSRmat[,i] = result[(n.time + 1):(2*n.time)]
    FSPhatmat[,i] = result[(2*n.time + 1):(3*n.time)]
    print(i/length(kvec))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(kvec)),
                  FSR = as.vector(FSRmat),
                  TSR = as.vector(TSRmat),
                  FSPhat = as.vector(FSPhatmat),
                  k = as.factor(rep(kvec, each = n.time)))
  
  write_csv(datare, 'sava_oracle_case_different_k.csv')
  
}

##### Comparison in Gaussian case #####
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu.v = seq(0.05,0.2,0.05)
n.rep = 10000# repeated number
gamma = 0
k = 25
w0 = q
outputlength = n.time*prob.start
#### SAVA ######
# change pi+
{
  mu = 0.1
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('dplyr','SAVA')) %dopar%{
                    re = sava_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = plrtio, mu = mu, gamma = 0, k = k, w0 = w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
  
  write_csv(datare, 'sava_gaussian_pi.csv')
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('dplyr','SAVA')) %dopar%{
      re = sava_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, gamma = 0, k = k, w0 = w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mu.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
  
  write_csv(datare, 'sava_gaussian_mu.csv')
}

#### Lordpp ####
# change pi+
{
  mu = 0.1
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('dplyr','SAVA')) %dopar%{
                    re = lordpp_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = plrtio, mu = mu, gamma = gamma)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
  
  write_csv(datare, 'lordpp_gaussian_pi.csv')
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('dplyr','SAVA')) %dopar%{
      re = lordpp_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, gamma = gamma)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
  
  write_csv(datare, 'lordpp_gaussian_mu.csv')
}
#### Saffron ####
# change pi+
{
  mu = 0.1
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('dplyr','SAVA')) %dopar%{
                    re = saffron_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = plrtio, mu = mu, gamma = gamma)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
  
  write_csv(datare, 'saffron_gaussian_pi.csv')
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('dplyr','SAVA')) %dopar%{
      re = saffron_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, gamma = gamma)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
 
  write_csv(datare, 'saffron_gaussian_mu.csv')
}
#### Addis #####
# change pi+
{
  mu = 0.1
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('dplyr','SAVA')) %dopar%{
                    re = addis_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = plrtio, mu = mu, gamma = gamma)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
 
  write_csv(datare, 'addis_gaussian_pi.csv')
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('dplyr','SAVA')) %dopar%{
      re = addis_gaussian(n.time = n.time, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, gamma = gamma)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
  
  write_csv(datare, 'addis_gaussian_mu.csv')
}

##### Comparison in truncated Guassian distribution ####
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
bound = 4
mu.v = c(0.8,1,1.2,1.4)
sigma = 1
k = 25
w0 = q
n.rep = 1000# repeated number
outputlength = n.time*prob.start
#### SAVA ######
# change pi+
{
  mu = 1
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q, ratio.plus = plrtio, mu = mu, sigma = sigma, gamma = 0, k = k, w0 = w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
 
  write_csv(datare, 'sava_general_pi.csv')
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, sigma = sigma, gamma = 0, k = k, w0 = w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mu.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
 
  write_csv(datare, 'sava_general_mu.csv')
}

#### Lordpp ####
# change pi
{
  mu = 1
  # prob.start = 0.03
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q, ratio.plus = plrtio,mu = mu,sigma = sigma, gamma = 0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
  
  write_csv(datare, 'lordpp_general_pi.csv')
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q, ratio.plus = 0.5,mu = mu,sigma = sigma, gamma = 0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
 
  write_csv(datare, 'lordpp_general_mu.csv')
}

#### Saffron ####
# change pi
{
  mu = 1
  # prob.start = 0.03
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q,ratio.plus = plrtio, mu = mu, sigma = sigma, gamma = 0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
 
  write_csv(datare, 'saffron_general_pi.csv')
  
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_general(n.time = n.time, bound = bound, prob.start = prob.start, q = q,ratio.plus = 0.5, mu = mu, sigma = sigma, gamma = 0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
  write_csv(datare, 'saffron_general_mu.csv')
  
}

#### Addis ####
# change pi
{
  mu = 1
  # prob.start = 0.03
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_general(n.time = n.time,bound =  bound, prob.start = prob.start, q = q,ratio.plus = plrtio, mu = mu, sigma = sigma, gamma = 0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = n.time)))
  write_csv(datare, 'addis_general_pi.csv')
  
  
}
# change mu
{
  FSR.mu = matrix(0, outputlength, length(mu.v))
  TSR.mu = matrix(0, outputlength, length(mu.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mu.v)){
    mu = mu.v[i]
    print(i/length(mu.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_general(n.time = n.time,bound =  bound, prob.start = prob.start, q = q,ratio.plus = 0.5, mu = mu, sigma = sigma, gamma = 0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:n.time, length(mu.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mu.v, each = n.time)))
  
  write_csv(datare, 'addis_general_mu.csv')
  
}



##### Comparison under Gamma model ####
# true mu+ = mu*mumulti and mu- = mu/mumulti for Gamma mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu = 4
mumulti.v = c(1.2,1.3,1.4,1.5)
k = 25
w0 = q
gammak = 2
n.rep = 1000# repeated number
outputlength = n.time*prob.start

# change pi+
{
  mumulti = 1.3
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_gamma(gammak, n.time, prob.start, plrtio, mu, q, mumulti, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_gamma_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mumulti.v))
  TSR.mu = matrix(0, outputlength, length(mumulti.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mumulti.v)){
    mumulti = mumulti.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_gamma(gammak, n.time, prob.start, ratio.plus, mu, q, mumulti, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mumulti.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mumulti.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  multiply = as.factor(rep(mumulti.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = multiply))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = multiply))
  write_csv(datare, 'sava_gamma_multiply.csv')
}
#### Lordpp ####
# change pi+
{
  mumulti = 1.3
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_gamma(gammak, n.time, prob.start, q, plrtio, mu, mumulti)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_gamma_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mumulti.v))
  TSR.mu = matrix(0, outputlength, length(mumulti.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mumulti.v)){
    mumulti = mumulti.v[i]
    print(i/length(mumulti.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_gamma(gammak, n.time, prob.start, q, ratio.plus, mu, mumulti)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mumulti.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mumulti.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_gamma_mumulti.csv')
}

#### Saffron ####
# change pi+
{
  mumulti = 1.3
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_gamma(gammak , n.time, prob.start, q, plrtio, mu, mumulti)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_gamma_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mumulti.v))
  TSR.mu = matrix(0, outputlength, length(mumulti.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mumulti.v)){
    mumulti = mumulti.v[i]
    print(i/length(mumulti.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_gamma(gammak , n.time, prob.start, q, ratio.plus, mu, mumulti)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mumulti.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mumulti.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_gamma_mumulti.csv')
  
}

#### Addis ####
# change pi+
{
  mumulti = 1.3
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(k = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_gamma(gammak, n.time, prob.start, q, plrtio, mu, mumulti)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_gamma_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mumulti.v))
  TSR.mu = matrix(0, outputlength, length(mumulti.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mumulti.v)){
    mumulti = mumulti.v[i]
    print(i/length(mumulti.v))
    ree = foreach(k = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_gamma(gammak, n.time, prob.start, q, ratio.plus, mu, mumulti)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mumulti.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mumulti.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_gamma_mumulti.csv')
  
}



##### Comparison under beta model ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0.5
mudelta.v = c(0.12,0.13,0.14,0.15)
k = 25
beta0 = 2
w0 = q
n.rep = 1000# repeated number
outputlength = n.time*prob.start
#### SAVA #####
# change pi+
{
  mudelta = 0.14
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_beta(n.time, prob.start, plrtio, mu0, mudelta, beta0, q, k , w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_beta_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_beta(n.time, prob.start, ratio.plus, mu0, mudelta, beta0, q, k , w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mudelta.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  delta = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = delta))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = delta))
  write_csv(datare, 'sava_beta_delta.csv')
}

#### Lordpp ####
# change pi+
{
  mudelta = 0.14
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_beta_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_beta_mu.csv')
}

#### Saffron ####
# change pi+
{
  mudelta = 0.14
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_beta_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_beta(n.time, prob.start, q, ratio.plus, mu0, mudelta, beta0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_beta_mu.csv')
  
}

#### Addis ####
# change pi+
{
  mudelta = 0.14
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_beta_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_beta(n.time, prob.start, q, ratio.plus, mu0, mudelta, beta0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_beta_mu.csv')
  
}
##### Comparison under beta model with coin betting e-process ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0.5
mudelta.v = c(0.04,0.05,0.06,0.07)
k = 25
beta0 = 2
w0 = q
n.rep = 1000# repeated number
outputlength = n.time*prob.start
##### SAVA #####
# change pi+
{
  mudelta = 0.05
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_beta_coin(n.time, prob.start, plrtio, mu0, mudelta, beta0, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_beta_pi_coin.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_beta_coin(n.time, prob.start, ratio.plus, mu0, mudelta, beta0, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mudelta.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  delta = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = delta))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = delta))
  write_csv(datare, 'sava_beta_delta_coin.csv')
}
#### Lordpp ####
# change pi+
{
  mudelta = 0.05
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_beta_coin_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_beta_coin_mu.csv')
}

#### Saffron ####
# change pi+
{
  mudelta = 0.05
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_beta_coin_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_beta(n.time, prob.start, q, ratio.plus, mu0, mudelta, beta0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_beta_coin_mu.csv')
  
}

#### Addis ####
# change pi+
{
  mudelta = 0.05
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_beta(n.time, prob.start, q, plrtio, mu0, mudelta, beta0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_beta_coin_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_beta(n.time, prob.start, q, ratio.plus, mu0, mudelta, beta0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_beta_coin_mu.csv')
  
}


##### Comparison with subgaussian e-processes ####
#### Gaussian model ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0
mudelta.v = c(0.4,0.45, 0.5,0.55)
k = 25
w0 = q
subsigma = 1
n.rep = 1000# repeated number
outputlength = n.time*prob.start
##### SAVA #####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_gauss_sub(subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_subgau_gauss_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_gauss_sub(subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mudelta.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  delta = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = delta))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = delta))
  write_csv(datare, 'sava_subgau_gauss_delta.csv')
}
#### Lordpp ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_gauss_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_sub_gauss_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_gauss_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_sub_gauss_mu.csv')
}

#### Saffron ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_gauss_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_sub_gauss_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_gauss_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_sub_gauss_mu.csv')
  
}

#### Addis ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_gauss_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_sub_gauss_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_gauss_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_sub_gauss_mu.csv')
  
}

#### Uniform model ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0
mudelta.v = c(0.45,0.5,0.55,0.6)
k = 25
w0 = q
subsigma = 1
n.rep = 1000# repeated number
outputlength = n.time*prob.start
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_unif_sub(subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_subgau_unif_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_unif_sub(subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mudelta.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  delta = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = delta))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = delta))
  write_csv(datare, 'sava_subgau_unif_delta.csv')
}
#### Lordpp ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_unif_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_sub_unif_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_unif_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_sub_unif_mu.csv')
}

#### Saffron ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_unif_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_sub_unif_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_unif_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_sub_unif_mu.csv')
  
}

#### Addis ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_unif_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_sub_unif_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_unif_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_sub_unif_mu.csv')
  
}
#### Comparison under drifted Bernoulli model ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0
mudelta.v = c(0.5,0.6,0.7,0.8)
k = 25
w0 = q
berp = 0.3
subsigma = 1
n.rep = 1000# repeated number
outputlength = n.time*prob.start
#### SAVA ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_ber_sub(berp, subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_subgau_ber_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_ber_sub(berp, subsigma = 1, n.time, prob.start, plrtio, mu0, mudelta, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(mudelta.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  delta = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = delta))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = delta))
  write_csv(datare, 'sava_subgau_ber_delta.csv')
}
#### Lordpp ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_ber_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_sub_ber_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_ber_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'lordpp_sub_ber_mu.csv')
}

#### Saffron ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_ber_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_sub_ber_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_ber_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'saffron_sub_ber_mu.csv')
  
}

#### Addis ####
# change pi+
{
  mudelta = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_ber_sub(subsigma, n.time, prob.start, q, plrtio, mu0, mudelta)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_sub_ber_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(mudelta.v))
  TSR.mu = matrix(0, outputlength, length(mudelta.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(mudelta.v)){
    mudelta = mudelta.v[i]
    print(i/length(mudelta.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_ber_sub(subsigma, n.time, prob.start, q, ratio.plus, mu0, mudelta)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(mudelta.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  mu = as.factor(rep(mudelta.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = mu))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = mu))
  write_csv(datare, 'addis_sub_ber_mu.csv')
  
}


##### SAVA under misspecified uniform model: change a #####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
a0 = 1
sigma0 = a0/sqrt(3)
mu0 = 0
a.v = c(0.5,1,1.25,1.5,2,4)
k = 25
w0 = q
subsigma = sigma0
n.rep = 1000# repeated number
outputlength = n.time*prob.start
ratio.plus = 0.5

# the case when mu_delta = 0.5
{
  mudelta = 0.5
  FSR.mu = matrix(0, outputlength, length(a.v))
  TSR.mu = matrix(0, outputlength, length(a.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(a.v)){
    aa = a.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_unif_sub_a(subsigma = 1, n.time, prob.start, ratio.plus, mu0, mudelta, aa, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(a.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(a.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  a = as.factor(rep(a.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = a))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = a))
  write_csv(datare, 'sava_subgau_unif_a.csv')
}
# the case when mu_delta = 1
{
  mudelta = 1
  FSR.mu = matrix(0, outputlength, length(a.v))
  TSR.mu = matrix(0, outputlength, length(a.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(a.v)){
    aa = a.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_unif_sub_a(subsigma = 1, n.time, prob.start, ratio.plus, mu0, mudelta, aa, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(a.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(a.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  a = as.factor(rep(a.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = a))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = a))
  write_csv(datare, 'sava_subgau_unif_a_mu1.csv')
}

#### SAVA under truncated Gaussian: change K ######
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
bound.v = c(4,6,8,10,12,20)
bound0 = 4
sigma = 1
k = 25
w0 = q
n.rep = 1000# repeated number
outputlength = n.time*prob.start
# the case when mu =1
{
  mu = 1
  FSR.mu = matrix(0, outputlength, length(bound.v))
  TSR.mu = matrix(0, outputlength, length(bound.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(bound.v)){
    bound = bound.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_truncgauss_K(n.time = n.time, bound = bound, bound0 = bound0, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, sigma = sigma, gamma = 0, k = k, w0 = w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(bound.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(bound.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  K = as.factor(rep(bound.v/2, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = K))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = K))
  write_csv(datare, 'sava_truncgauss_K.csv')
}

# the case when mu = 1.5

{
  mu = 1.5
  FSR.mu = matrix(0, outputlength, length(bound.v))
  TSR.mu = matrix(0, outputlength, length(bound.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(bound.v)){
    bound = bound.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_truncgauss_K(n.time = n.time, bound = bound, bound0 = bound0, prob.start = prob.start, q = q, ratio.plus = 0.5, mu = mu, sigma = sigma, gamma = 0, k = k, w0 = w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(bound.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(bound.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  K = as.factor(rep(bound.v/2, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = K))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = K))
  write_csv(datare, 'sava_truncgauss_K_mu0.5.csv')
}




##### comparison with dependent data streams ####

#### Gaussian model ####
# true mu+ = mu0 + mudelta and mu- = mu0 - mudelta for sample mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu0 = 0
k = 25
w0 = q
subsigma = 1
rho.v = c(0.2,0.3,0.4,0.5)
n.rep = 1000# repeated number
outputlength = n.time*prob.start
#### SAVA ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_depend_gauss_sub(rho, subsigma, n.time, prob.start, plrtio, mu0, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_depend_subgau_gauss_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_depend_gauss_sub(rho, subsigma, n.time, prob.start, ratio.plus, mu0, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(rho.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'sava_depend_subgau_gauss_rho.csv')
}
#### Lordpp ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_dependgauss_sub(rho, subsigma, n.time, prob.start, q, plrtio, mu0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_depend_sub_gauss_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_dependgauss_sub(rho, subsigma, n.time, prob.start, q, ratio.plus, mu0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'lordpp_depend_sub_gauss_mu.csv')
}

#### Saffron ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_dependgauss_sub(rho, subsigma, n.time, prob.start, q, plrtio, mu0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_depend_sub_gauss_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_dependgauss_sub(rho, subsigma, n.time, prob.start, q, ratio.plus, mu0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'saffron_depend_sub_gauss_mu.csv')
  
}

#### Addis ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_dependgauss_sub(rho, subsigma, n.time, prob.start, q, plrtio, mu0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_depend_sub_gauss_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(4)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_dependgauss_sub(rho, subsigma, n.time, prob.start, q, ratio.plus, mu0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'addis_depend_sub_gauss_mu.csv')
  
}

##### Gamma model ####
# true mu+ = mu*mumulti and mu- = mu/mumulti for Gamma mean
n.time = 3000 # number of total times
prob.start = 1/3# probability of another hypothesis starting to be observed
q = 0.05 # significant level
ratioplus = seq(0.2,0.8,0.2)
mu = 4
rho.v = c(0.2,0.4,0.6,0.8)
k = 25
w0 = q
gammak = 2
n.rep = 1000# repeated number
outputlength = n.time*prob.start
#### SAVA ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6611)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = sava_depend_gamma(rho, gammak, n.time, prob.start, plrtio, mu, q, k, w0)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'sava_depend_gamma_pi.csv')
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(668)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = sava_depend_gamma(rho, gammak, n.time, prob.start, ratio.plus, mu, q, k, w0)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
    print(i/length(rho.v))
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  multiply = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = multiply))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = multiply))
  write_csv(datare, 'sava_depend_gamma_multiply.csv')
}


#### Lordpp ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(66)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_dependgamma(rho, gammak, n.time, prob.start, q, plrtio, mu)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'lordpp_depend_gamma_pi.csv')
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = lordpp_dependgamma(rho, gammak, n.time, prob.start, q, ratio.plus, mu)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'lordpp_depend_gamma_mu.csv')
}

#### Saffron ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = saffron_dependgamma(rho, gammak, n.time, prob.start, q, plrtio, mu)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'saffron_depend_gamma_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = saffron_dependgamma(rho, gammak, n.time, prob.start, q, ratio.plus, mu)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'saffron_depend_gamma_mu.csv')
  
}

#### Addis ####
# change pi+
{
  rho = 0.5
  FSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  TSR.ratioplus = matrix(NA, outputlength, length(ratioplus))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(64)
  for(i in 1:length(ratioplus)){
    plrtio = ratioplus[i]
    print(i/length(ratioplus))
    ree = foreach(kk = 1:n.rep, .combine = cbind,
                  .packages = c('truncnorm','SAVA')) %dopar%{
                    re = addis_dependgamma(rho, gammak, n.time, prob.start, q, plrtio, mu)
                    fsrvec = re$FSP
                    if (length(fsrvec) >= outputlength){
                      fsrvec = fsrvec[1:outputlength]
                    }else{
                      fsrvec = c(fsrvec, 
                                 rep(NA,(outputlength - length(fsrvec))))
                    }
                    tprtvec = re$TP / sapply(re$realsig, max, 1)
                    if (length(tprtvec) >= outputlength){
                      tprtvec = tprtvec[1:outputlength]
                    }else{
                      tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
                    }
                    c(fsrvec, tprtvec)
                  }
    result = rowMeans(ree, na.rm = T)
    FSR.ratioplus[,i] = result[1:outputlength]
    TSR.ratioplus[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(ratioplus)),
                  FSR = as.vector(FSR.ratioplus),
                  TSR = as.vector(TSR.ratioplus),
                  pi = as.factor(rep(ratioplus, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = pi))
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = pi))
  write_csv(datare, 'addis_depend_gamma_pi.csv')
  
  
}
# change mu
{
  ratio.plus = 0.5
  FSR.mu = matrix(0, outputlength, length(rho.v))
  TSR.mu = matrix(0, outputlength, length(rho.v))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  set.seed(6689)
  for(i in 1:length(rho.v)){
    rho = rho.v[i]
    print(i/length(rho.v))
    ree = foreach(kk = 1:n.rep, .combine = cbind, .packages = c('truncnorm','SAVA')) %dopar%{
      re = addis_dependgamma(rho, gammak, n.time, prob.start, q, ratio.plus, mu)
      fsrvec = re$FSP
      if (length(fsrvec) >= outputlength){
        fsrvec = fsrvec[1:outputlength]
      }else{
        fsrvec = c(fsrvec, 
                   rep(NA,(outputlength - length(fsrvec))))
      }
      tprtvec = re$TP / sapply(re$realsig, max, 1)
      if (length(tprtvec) >= outputlength){
        tprtvec = tprtvec[1:outputlength]
      }else{
        tprtvec = c(tprtvec, rep(NA, outputlength - length(tprtvec)))
      }
      c(fsrvec, tprtvec)
    }
    result = rowMeans(ree, na.rm = T)
    FSR.mu[,i] = result[1:outputlength]
    TSR.mu[,i] = result[-(1:outputlength)]
  }
  stopImplicitCluster()
  stopCluster(cl)
  datare = tibble(time = rep(1:outputlength, length(rho.v)),
                  FSR = as.vector(FSR.mu),
                  TSR = as.vector(TSR.mu),
                  rho = as.factor(rep(rho.v, each = outputlength)))
  ggplot(datare)+
    geom_line(aes(x = time, y = FSR, color = rho))
  
  ggplot(datare)+
    geom_line(aes(x = time, y = TSR, color = rho))
  write_csv(datare, 'addis_depend_gamma_mu.csv')
  
}
