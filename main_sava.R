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
  cl <- makeCluster(10)
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

##### Comparison in general distribution ####
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
  cl <- makeCluster(10)
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
  cl <- makeCluster(10)
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


