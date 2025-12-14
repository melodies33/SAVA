library(truncnorm)

levelgenek = function(n){
  re = rep(0,n)
  for (i in 1:(n-1)){
    re[i] = 0.5^i
  }
  re[n] = 0.5^(n-1)
  return(re)
}

testlevelperm = function(indcur, deltaesti, stoptime, t, startindcur, k, w0, q){
  levelpool = levelgenek(k)
  if (indcur <= k){
    orderind = (t - startindcur)%%k + 1
    level = w0*levelpool[orderind]
  }else{
    level = 0
  }
  selectindex = which(deltaesti != 0)
  selectindex = selectindex[which(selectindex < indcur)]
  selectstoptime = stoptime[selectindex]
  ll = length(selectstoptime)
  if (ll > 0){
    ordersel = (t - selectstoptime[1])%%k + 1
    level = level + (q - w0)*levelpool[ordersel]*(indcur <= selectindex[1] + k)
    if (ll > 1){
      for (i in 2:length(selectstoptime)){
        ordersel = (t - selectstoptime[i])%%k + 1
        level = level + q*levelpool[ordersel]*(indcur <= selectindex[i] + k)
      }
    }
  }
  return(level)
}

# function for estimate the proportion of directions
propesti = function(issparvec, tau){
  nums = length(issparvec)
  res = sum(issparvec)/(1-tau)/nums
  res = min(res, 1)
  return(res)
}

falsetotalnum = function(realvec, estivec){
  numberplus = sum((realvec - 1)*estivec < 0)
  numberminus = sum((realvec + 1)*estivec < 0)
  numbers = numberplus + numberminus
  return(numbers)
}

falseplusnum = function(realvec, estivec){
  numbers = sum((realvec - 1) * estivec < 0)
  return(numbers)
}


falseminusnum = function(realvec, estivec){
  numbers = sum((1 + realvec)*estivec < 0)
}
trueplusnum = function(realvec, estivec){
  plusind = which(realvec > 0)
  numbers = sum((realvec * estivec)[plusind] == 1)
  return(numbers)
}

trueminusnum = function(realvec, estivec){
  minusind = which(realvec < 0)
  numbers = sum((realvec * estivec)[minusind] == 1)
  return(numbers)
}
truetotalnum = function(realvec, estivec){
  plusind = which(realvec > 0)
  minusind = which(realvec < 0)
  numplus = sum((realvec * estivec)[plusind] == 1)
  numminus = sum((realvec * estivec)[minusind] == 1)
  return(numplus + numminus)
}
###### functions for counterexamples #####

sava_counterexample1 = function(n.time = 500, # number of total times
                                   mu.v,
                                   q = 0.1 # FDR control
){
  starttime = seq(1, n.time)
  numstream = n.time
  realtotal = seq(1,n.time,1)
  thetareal = rep(NA, n.time)
  thetareal[which(mu.v > 0)] = 1
  thetareal[which(mu.v < 0)] = -1
  pvalueplusseq = runif(numstream)
  pvalueminusseq = runif(numstream)
  for (i in 1:numstream){
    if (mu.v[i] > 0){
      pvalueplusseq[i] = 1 - pnorm(q = rnorm(1, mu.v[i], 1), mean = 0, sd = 1)
    }else{
      pvalueminusseq[i] = pnorm(q = rnorm(1, mu.v[i], 1), mean = 0, sd = 1)
    }
  }
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, n.time)
  decidetotal = rep(0, n.time)
  truetotal = rep(0, n.time)
  testlevelvec = rep(0, n.time)
  wealth = q
  selectnum = 0
  for (t in 1:n.time){
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[t] = falsetotal[t-1]
      decidetotal[t] = decidetotal[t-1]
      truetotal[t] = truetotal[t-1]
    }else{
      for (a in 1:length(active)){
        indcur = active[a]
        threshold = testlevelvec[indcur]
        deltaplus = pvalueplusseq[indcur] <= threshold
        deltaminus = pvalueminusseq[indcur] <= threshold
        isplus = pvalueplusseq[indcur] < pvalueminusseq[indcur]
        deltatotal = deltaplus + deltaminus
        if (deltatotal != 0){
          selectnum = selectnum + 1
          if (selectnum > 1){
            wealth = wealth + q
          }
          if (deltatotal == 1){
            deltaesti[indcur] = deltaplus - deltaminus
          }else if (deltatotal == 2){
            if (isplus){
              deltaesti[indcur] = 1
            }else{
              deltaesti[indcur] = -1
            }
          }
          stoptime[indcur] = t
        }
      }
      falsetotal[t] = falsetotalnum(thetareal, deltaesti)
      decidetotal[t] = sum(abs(deltaesti))
      truetotal[t] = truetotalnum(thetareal, deltaesti)
      
      for(ii in 1:t){
        if (abs(deltaesti)[ii] == 0){
          avp = min(pvalueminusseq[ii],pvalueplusseq[ii])
          if (avp <= wealth){
            testlevelvec[ii] = avp
            wealth = wealth - avp
          }
        }
      }
      
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal,
              fsnumber = falsetotal,
              decidenumber = sapply(decidetotal, max, 1)))
}

sava_ctexam1_compare = function(n.time = 500,
                                mu.v,
                                q = 0.1,
                                w0, k){
  # zetasum = zeta(1.4, shift = 1)
  starttime = seq(1, n.time)
  numstream = n.time
  realtotal = seq(1,n.time,1)
  thetareal = rep(NA, n.time)
  thetareal[which(mu.v > 0)] = 1
  thetareal[which(mu.v < 0)] = -1
  pvalueplusseq = runif(numstream)
  pvalueminusseq = runif(numstream)
  for (i in 1:numstream){
    if (mu.v[i] > 0){
      pvalueplusseq[i] = 1 - pnorm(q = rnorm(1, mu.v[i], 1), mean = 0, sd = 1)
    }else{
      pvalueminusseq[i] = pnorm(q = rnorm(1, mu.v[i], 1), mean = 0, sd = 1)
    }
  }
  selectlarger1 = 0
  thresholdvec = rep(0, numstream)
  thresholdvec[1:k] = w0/k
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, n.time)
  decidetotal = rep(0, n.time)
  truetotal = rep(0, n.time)
  rejvec = integer(0)
  for (t in 1:n.time){
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[t] = falsetotal[t-1]
      decidetotal[t] = decidetotal[t-1]
      truetotal[t] = truetotal[t-1]
    }else{
      for (a in 1:length(active)){
        indcur = active[a]
        threshold = thresholdvec[indcur]
        deltaplus = pvalueplusseq[indcur] < threshold
        deltaminus = pvalueminusseq[indcur] < threshold
        isplus = pvalueplusseq[indcur] < pvalueminusseq[indcur]
        deltatotal = deltaplus + deltaminus
        if (deltatotal != 0){
          if (deltatotal == 1){
            deltaesti[indcur] = deltaplus - deltaminus
          }else if (deltatotal == 2){
            if (isplus){
              deltaesti[indcur] = 1
            }else{
              deltaesti[indcur] = -1
            }
          }
          stoptime[indcur] = t
          selectlarger1 = selectlarger1 + 1
          endindex = min(indcur+k, numstream)
          thresholdvec[(indcur + 1):endindex] = thresholdvec[(indcur + 1):endindex] + (selectlarger1 == 1)*(q - w0)/k + (selectlarger1 > 1)*q/k
        }
      }
      falsetotal[t] = falsetotalnum(thetareal, deltaesti)
      decidetotal[t] = sum(abs(deltaesti))
      truetotal[t] = truetotalnum(thetareal, deltaesti)
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal,
              fsnumber = falsetotal,
              decidenumber = sapply(decidetotal, max, 1)))
}

sava_counterexample2 = function(n.time = 500,# number of total times
                                prob.start = 1,
                                ratio.plus = 0.5,# probability of another hypothesis starting to be observed
                                q = 0.1,# FDR control
                                mu = 2, w0, k){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)
  numplus = length(which(thetareal == 1))
  numminus = length(which(thetareal == -1))
  selectlarger1 = 0
  pvalueplusseq = runif(numstream)
  pvalueminusseq = runif(numstream)
  
  if (numplus > 0){
    pvalueplusseq[which(thetareal == 1)] = 1 - sapply(rnorm(length(which(thetareal == 1)),mu,0.5),pnorm, mean = 0, sd = 1)
  }
  if (numminus > 0){
    pvalueminusseq[which(thetareal == -1)] = sapply(rnorm(length(which(thetareal == -1)),-mu,0.5),pnorm, mean = 0, sd = 1)
  }
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, n.time)
  decidetotal = rep(0, n.time)
  truetotal = rep(0, n.time)
  for (t in 1:n.time){
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[t] = falsetotal[t-1]
      decidetotal[t] = decidetotal[t-1]
      truetotal[t] = truetotal[t-1]
    }else{
      for (a in 1:length(active)){
        indcur = active[a]
        threshold = testlevelperm(indcur, deltaesti, stoptime, t, starttime[indcur], k, w0, q)
        deltaplus = pvalueplusseq[indcur] < threshold
        deltaminus = pvalueminusseq[indcur] < threshold
        isplus = pvalueplusseq[indcur] < pvalueminusseq[indcur]
        deltatotal = deltaplus + deltaminus
        if (deltatotal != 0){
          if (deltatotal == 1){
            deltaesti[indcur] = deltaplus - deltaminus
          }else if (deltatotal == 2){
            if (isplus){
              deltaesti[indcur] = 1
            }else{
              deltaesti[indcur] = -1
            }
          }
          stoptime[indcur] = t
        }
      }
      falsetotal[t] = falsetotalnum(thetareal, deltaesti)
      decidetotal[t] = sum(abs(deltaesti))
      truetotal[t] = truetotalnum(thetareal, deltaesti)
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal,
              fsnumber = falsetotal,
              decidenumber = sapply(decidetotal, max, 1)))
}

sava_ctexam2_compare = function(n.time = 500, 
                              prob.start = 1,
                              ratio.plus = 0.5,
                              q = 0.1, 
                              mu = 2, w0, k){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)
  numplus = length(which(thetareal == 1))
  numminus = length(which(thetareal == -1))
  selectlarger1 = 0
  thresholdvec = rep(0, numstream)
  thresholdvec[1:k] = w0*levelgenek(k)
  pvalueplusseq = runif(numstream)
  pvalueminusseq = runif(numstream)
  
  if (numplus > 0){
    pvalueplusseq[which(thetareal == 1)] = 1 - sapply(rnorm(length(which(thetareal == 1)),mu,0.5),pnorm, mean = 0, sd = 1)
  }
  if (numminus > 0){
    pvalueminusseq[which(thetareal == -1)] = sapply(rnorm(length(which(thetareal == -1)),-mu,0.5),pnorm, mean = 0, sd = 1)
  }
  
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, n.time)
  decidetotal = rep(0, n.time)
  truetotal = rep(0, n.time)
  selectvec = integer(0)
  for (t in 1:n.time){
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[t] = falsetotal[t-1]
      decidetotal[t] = decidetotal[t-1]
      truetotal[t] = truetotal[t-1]
    }else{
      samptime = t - starttime[active] + 1
      for (a in 1:length(active)){
        indcur = active[a]
        timecur = samptime[a]
        threshold = thresholdvec[indcur]
        deltaplus = pvalueplusseq[indcur] < threshold
        deltaminus = pvalueminusseq[indcur] < threshold
        isplus = pvalueplusseq[indcur] < pvalueminusseq[indcur]
        deltatotal = deltaplus + deltaminus
        if (deltatotal != 0){
          if (deltatotal == 1){
            deltaesti[indcur] = deltaplus - deltaminus
          }else if (deltatotal == 2){
            if (isplus){
              deltaesti[indcur] = 1
            }else{
              deltaesti[indcur] = -1
            }
          }
          stoptime[indcur] = t
          selectlarger1 = selectlarger1 + 1
          endindex = min(indcur+k, numstream)
          thresholdvec[(indcur + 1):endindex] = thresholdvec[(indcur + 1):endindex] + (selectlarger1 == 1)*(q - w0)*levelgenek(k)[1:(endindex - indcur)] + (selectlarger1 > 1)*q*levelgenek(k)[1:(endindex - indcur)]
        }
      }
      falsetotal[t] = falsetotalnum(thetareal, deltaesti)
      decidetotal[t] = sum(abs(deltaesti))
      truetotal[t] = truetotalnum(thetareal, deltaesti)
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal,
              fsnumber = falsetotal,
              decidenumber = sapply(decidetotal, max, 1)))
}


##### algorithms in general cases ######

sava_general = function(n.time = 500, bound = 4, prob.start = 1, q = 0.1, ratio.plus = 0.5, mu = 0.2, sigma = 1, gamma = 0, k = 20, w0 = 0.05){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  samplemat = matrix(NA, numstream, n.time)
  selectvec = integer(0)
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  # generate samples
  for (i in 1:numstream){
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    samplemat[i, ] = rtruncnorm(n.time, -bound/2, bound/2, mureal, sigma)
  }
 
  selectlarger1 = 0
  thresholdvec = rep(0, numstream)
  thresholdvec[1:k] = w0/k
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (ti in 1:num_decision){
    t = decisiontime[ti]
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[ti] = falsetotal[ti-1]
      decidetotal[ti] = decidetotal[ti-1]
      truetotal[ti] = truetotal[ti-1]
    }else{
      samptime = t - starttime[active] + 1
      for (a in 1:length(active)){
        indcur = active[a]
        timecur = samptime[a]
        threshold = thresholdvec[indcur]
        cs.below = falphaplus(threshold, samplemat[indcur, 1:timecur], bound, timecur, gamma)
        cs.above = falphaminus(threshold, samplemat[indcur, 1:timecur], bound, timecur, gamma)
        
        deltaplus = cs.below >= 0
        deltaminus = cs.above < 0
        deltatotal = deltaplus - deltaminus
        if (deltatotal != 0){
          deltaesti[indcur] = deltatotal
          stoptime[indcur] = t
          selectlarger1 = selectlarger1 + 1
          endindex = min(indcur+k, numstream)
          thresholdvec[(indcur + 1):endindex] = thresholdvec[(indcur + 1):endindex] + (selectlarger1 == 1)*(q - w0)/k + (selectlarger1 > 1)*q/k
        }
      }
      falsetotal[ti] = falsetotalnum(thetareal, deltaesti)
      decidetotal[ti] = sum(abs(deltaesti))
      truetotal[ti] = truetotalnum(thetareal, deltaesti)
    }
  }
  
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}


lordpp_general = function(n.time = 500,  bound = 4,  prob.start = 1,  q = 0.1, ratio.plus = 0.5, mu = 0.2, sigma = 5, gamma = 0){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    sampnum = edtime - sttime + 1
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rtruncnorm(sampnum, -bound/2, bound/2, mureal, sigma)
    thresholdposi = threldpp(q,q/10, i, selectposi)
    thresholdnega = threldpp(q,q/10, i, selectnega)
    p.posi = wilcox.test(sampvec, mu = 0, alternative = 'greater')$p.value
    p.nega = wilcox.test(sampvec, mu = 0, alternative = 'less')$p.value
    deltaplus = p.posi < thresholdposi
    deltaminus = p.nega < thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
    }
    
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
    
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}

saffron_general = function(n.time = 3000, bound = 4, prob.start = 0.05,  q = 0.1,  ratio.plus = 0.5,  mu = 2, sigma = 1, gamma = 0.5){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  cplus = rep(0, numstream)
  cminus = rep(0, numstream)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    sampnum = edtime - sttime + 1
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rtruncnorm(sampnum, -bound/2, bound/2, mureal, sigma)
    
    p.posi = wilcox.test(sampvec, mu = 0, alternative = 'greater')$p.value
    p.nega = wilcox.test(sampvec, mu = 0, alternative = 'less')$p.value
    
    if (p.posi <= 0.5){
      cplus[i] = 1
    }
    if (p.nega <= 0.5){
      cminus[i] = 1
    }
    thresholdposi = thresaff(q, q/5, i, selectposi, cplus, 0.5)
    thresholdnega = thresaff(q,q/5, i, selectnega, cminus, 0.5)
    
    deltaplus = p.posi <= thresholdposi
    deltaminus = p.nega <= thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
    }
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
    
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}

addis_general = function(n.time = 500, bound = 4, prob.start = .3, q = 0.1, ratio.plus = 0.5, mu = 0.2, sigma = 5, gamma = 0){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  kstarposi = integer(0)
  kstarnega = integer(0)
  splus = rep(0, numstream)
  sminus = rep(0, numstream)
  cplus = rep(0, numstream)
  cminus = rep(0, numstream)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    sampnum = edtime - sttime + 1
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rtruncnorm(sampnum, -bound/2, bound/2, mureal, sigma)
    p.posi = wilcox.test(sampvec, mu = 0, alternative = 'greater')$p.value
    p.nega = wilcox.test(sampvec, mu = 0, alternative = 'less')$p.value
    
    if (p.posi <= 0.5){
      splus[i] = 1
    }
    if (p.nega <= 0.5){
      sminus[i] = 1
    }
    if (p.posi <= 0.25){
      cplus[i] = 1
    }
    if (p.nega <= 0.25){
      cminus[i] = 1
    }
    
    sindexposi = sum(splus)
    sindexnega = sum(sminus)
    
    thresholdposi = threaddis(q, q/2, i, selectposi, cplus, sindexposi, kstarposi, lamb = 0.25, tau = 0.5)
    thresholdnega = threaddis(q, q/2, i, selectnega, cminus, sindexnega, kstarnega, lamb = 0.25, tau = 0.5)
    
    deltaplus = p.posi <= thresholdposi
    deltaminus = p.nega <= thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
      kstarposi = append(kstarposi, sum(splus))
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
      kstarnega = append(kstarnega, sum(sminus))
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
     
    }
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}

##### Algorithms in Gaussian case #####
sava_gaussian = function(n.time = 3000, prob.start = 0.05, ratio.plus = 0.5, mu = 2, q = 0.1, gamma = 0, k = 10, w0){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  samplemat = matrix(NA, numstream, n.time)
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  # generate samples
  for (i in 1:numstream){
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    samplemat[i, ] = rnorm(n.time, mureal, 1)
  }
  sprtmat = matrix(NA, numstream, n.time)
  for (i in 1:numstream){
    sprtmat[i, ] = cummean(samplemat[i, ])
  }
  selectlarger1 = 0
  thresholdvec = rep(0, numstream)
  thresholdvec[1:k] = w0/k
  deltaesti = rep(0, numstream)
  stoptime = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (ti in 1:num_decision){
    t = decisiontime[ti]
    active = which((starttime <= t) & (stoptime == 0))
    if (length(active) == 0){
      falsetotal[ti] = falsetotal[ti-1]
      decidetotal[ti] = decidetotal[ti-1]
      truetotal[ti] = truetotal[ti-1]
    }else{
      samptime = t - starttime[active] + 1
      for (a in 1:length(active)){
        indcur = active[a]
        timecur = samptime[a]
        threshold = thresholdvec[indcur]
        lamposi = oracsprt(mu, -mu, sprtmat[indcur, timecur], timecur)
        lamnega = oracsprt(-mu, mu, sprtmat[indcur, timecur], timecur)
        deltaplus = (lamposi >= 1/threshold)
        deltaminus = (lamnega >= 1/threshold)
        deltatotal = deltaplus - deltaminus
        if (deltatotal != 0){
          deltaesti[indcur] = deltatotal
          stoptime[indcur] = t
          selectlarger1 = selectlarger1 + 1
          endindex = min(indcur+k, numstream)
          thresholdvec[(indcur + 1):endindex] = thresholdvec[(indcur + 1):endindex] + (selectlarger1 == 1)*(q - w0)/k + (selectlarger1 > 1)*q/k
        }
      }
      falsetotal[ti] = falsetotalnum(thetareal, deltaesti)
      decidetotal[ti] = sum(abs(deltaesti))
      truetotal[ti] = truetotalnum(thetareal, deltaesti)
    }
  }
  # for (t in 1:n.time){
  #   active = which((starttime <= t) & (stoptime == 0))
  #   if (length(active) == 0){
  #     falsetotal[t] = falsetotal[t-1]
  #     decidetotal[t] = decidetotal[t-1]
  #     truetotal[t] = truetotal[t-1]
  #   }else{
  #     samptime = t - starttime[active] + 1
  #     for (a in 1:length(active)){
  #       indcur = active[a]
  #       timecur = samptime[a]
  #       threshold = thresholdvec[indcur]
  #       lamposi = oracsprt(mu, -mu, sprtmat[indcur, timecur], timecur)
  #       lamnega = oracsprt(-mu, mu, sprtmat[indcur, timecur], timecur)
  #       deltaplus = (lamposi >= 1/threshold)
  #       deltaminus = (lamnega >= 1/threshold)
  #       deltatotal = deltaplus - deltaminus
  #       if (deltatotal != 0){
  #         deltaesti[indcur] = deltatotal
  #         stoptime[indcur] = t
  #         selectlarger1 = selectlarger1 + 1
  #         endindex = min(indcur+k, numstream)
  #         thresholdvec[(indcur + 1):endindex] = thresholdvec[(indcur + 1):endindex] + (selectlarger1 == 1)*(q - w0)/k + (selectlarger1 > 1)*q/k
  #       }
  #     }
  #     falsetotal[t] = falsetotalnum(thetareal, deltaesti)
  #     decidetotal[t] = sum(abs(deltaesti))
  #     truetotal[t] = truetotalnum(thetareal, deltaesti)
  #   }
  # }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}


lordpp_gaussian = function(n.time, prob.start, q, ratio.plus, mu, gamma){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    # if (i < numstream){
    #   times = seq(sttime, starttime[i + 1] - 1,1)
    #   edtime = starttime[i + 1] - 1
    # }else{
    #   times = seq(sttime, n.time)
    #   edtime = n.time
    # }
    sampnum = edtime - sttime + 1
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rnorm(sampnum, mureal, 1)
    sampmean = mean(sampvec)
    uvalue = sqrt(sampnum)*sampmean
    p.posi = pnorm(uvalue, 0, 1, lower.tail = F)
    p.nega = pnorm(uvalue, 0, 1)
    thresholdposi = threldpp(q,q/5, i, selectposi)
    thresholdnega = threldpp(q,q/5, i, selectnega)
    deltaplus = p.posi <= thresholdposi
    deltaminus = p.nega <= thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
    }
    # falsetotal[times] = falsetotal[sttime]
    # decidetotal[times] = decidetotal[sttime]
    # truetotal[times] = truetotal[sttime]
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}

saffron_gaussian = function(n.time, prob.start, q, ratio.plus, mu, gamma){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus) + 1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  cplus = rep(0, numstream)
  cminus = rep(0, numstream)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    # if (i < numstream){
    #   times = seq(sttime, starttime[i + 1] - 1,1)
    #   edtime = starttime[i + 1] - 1
    # }else{
    #   times = seq(sttime, n.time)
    #   edtime = n.time
    # }
    sampnum = edtime - sttime + 1
    # if (i < numstream){
    #   times = seq(sttime, starttime[i + 1] - 1,1)
    #   edtime = starttime[i + 1] - 1
    # }else{
    #   times = seq(sttime, n.time)
    #   edtime = n.time
    # }
    # sampnum = edtime - sttime + 1
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rnorm(sampnum, mureal, 1)
    sampmean = mean(sampvec)
    uvalue = sqrt(sampnum)*sampmean
    p.posi = pnorm(uvalue, 0, 1, lower.tail = F)
    p.nega = pnorm(uvalue, 0, 1)
    if (p.posi <= 0.5){
      cplus[i] = 1
    }
    if (p.nega <= 0.5){
      cminus[i] = 1
    }
    thresholdposi = thresaff(q, q/5, i, selectposi, cplus, 0.5)
    thresholdnega = thresaff(q,q/5, i, selectnega, cminus, 0.5)
    deltaplus = p.posi <= thresholdposi
    deltaminus = p.nega <= thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
    }
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
    
    
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}

addis_gaussian = function(n.time, prob.start, q, ratio.plus, mu, gamma){
  newtime = runif(n.time - 1) < prob.start
  starttime = c(1, which(newtime == 1) + 1)
  decisiontime = starttime[-1] - 1
  if (decisiontime[length(decisiontime)] < n.time){
    decisiontime = c(decisiontime, n.time)
  }
  num_decision = length(decisiontime)
  numstream = length(starttime)
  thetareal = -2*as.numeric(runif(numstream) > ratio.plus)+1
  
  realtotal = rep(0, n.time)
  realtotal[starttime[which(thetareal != 0)]] = 1
  realtotal = cumsum(realtotal)[decisiontime]
  selectposi = integer(0)
  selectnega = integer(0)
  kstarposi = integer(0)
  kstarnega = integer(0)
  splus = rep(0, numstream)
  sminus = rep(0, numstream)
  cplus = rep(0, numstream)
  cminus = rep(0, numstream)
  deltaesti = rep(0, numstream)
  falsetotal = rep(0, num_decision)
  decidetotal = rep(0, num_decision)
  truetotal = rep(0, num_decision)
  for (i in 1:numstream){
    sttime = starttime[i]
    edtime = decisiontime[i]
    sampnum = edtime - sttime + 1
    # sttime = starttime[i]
    # if (i < numstream){
    #   times = seq(sttime, starttime[i + 1] - 1,1)
    #   edtime = starttime[i + 1]
    # }else{
    #   times = seq(sttime, n.time)
    #   edtime = n.time
    # }
    # sampnum = edtime - sttime 
    if (thetareal[i] == 1){
      mureal = mu
    }else{
      mureal = - mu
    }
    sampvec = rnorm(sampnum, mureal, 1)
    sampmean = mean(sampvec)
    uvalue = sqrt(sampnum)*sampmean
    p.posi = pnorm(uvalue, 0, 1, lower.tail = F)
    p.nega = pnorm(uvalue, 0, 1)
    if (p.posi <= 0.5){
      splus[i] = 1
    }
    if (p.nega <= 0.5){
      sminus[i] = 1
    }
    
    if (p.posi <= 0.25){
      cplus[i] = 1
    }
    if (p.nega <= 0.25){
      cminus[i] = 1
    }
    sindexposi = sum(splus)
    sindexnega = sum(sminus)
    thresholdposi = threaddis(q, q/5, i, selectposi, cplus, sindexposi, kstarposi, lamb = 0.25, tau = 0.5)
    thresholdnega = threaddis(q, q/5, i, selectnega, cminus, sindexnega, kstarnega, lamb = 0.25, tau = 0.5)
   
    deltaplus = p.posi <= thresholdposi
    deltaminus = p.nega <= thresholdnega
    if (deltaplus > 0){
      selectposi = append(selectposi, i)
      kstarposi = append(kstarposi, sum(splus))
    }
    if (deltaminus > 0){
      selectnega = append(selectnega, i)
      kstarnega = append(kstarnega, sum(sminus))
    }
    deltatotal = deltaplus - deltaminus
    if (deltatotal != 0){
      deltaesti[i] = deltatotal
    }
    # falsetotal[times] = falsetotal[sttime]
    # decidetotal[times] = decidetotal[sttime]
    # truetotal[times] = truetotal[sttime]
    if (i > 1){
      falsetotal[i] = falsetotal[i-1] + falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = decidetotal[i-1] + abs(deltaesti[i])
      truetotal[i] = truetotal[i-1] + truetotalnum(thetareal[i], deltaesti[i])
    }else{
      falsetotal[i] = falsetotalnum(thetareal[i], deltaesti[i])
      decidetotal[i] = abs(deltaesti[i])
      truetotal[i] = truetotalnum(thetareal[i], deltaesti[i])
    }
  }
  fsp = falsetotal / sapply(decidetotal, max, 1)
  return(list(FSP = fsp,
              TP = truetotal,
              realsig = realtotal))
}
