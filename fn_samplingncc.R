
sampling.ncc.maxm = function(cohort, m, tau){
  
  nfull = nrow(cohort)
  cohort$ind.fail[cohort$eventime>tau] <- FALSE
  cohort_old = cohort_new = cohort
  
  ph2id_new = ph2id_old = matrix(NA, nrow = sum(cohort_old$ind.fail), ncol = m + 1, 
                                 dimnames = list(1:sum(cohort_old$ind.fail), c("case", paste("ctrl",1:m, sep=""))))
  
  case_id = which(cohort_old$ind.fail)
  cohort_old$nrisk = NA 
  nriskM = matrix(NA, nrow = nfull, ncol = m); colnames(nriskM) = paste("nrisk.",1:m,sep="")

  
  for(i in 1:nfull){
    atrisk_old = with(cohort_old, which(eventime >= eventime[i]))
    atrisk_new = setdiff(atrisk_old, as.numeric(na.omit(c(ph2id_new))))
    
    cohort_old$nrisk[i] = length(atrisk_old)
    for(mk in 1:m){
      tmpff = setdiff(atrisk_old, as.numeric(na.omit(c(ph2id_new[,1:(mk+1)]))))
      nriskM[i,mk] = length(tmpff)
    }
    
    if(i %in% case_id){
      j = which(case_id %in% i)
      ph2id_old[j, 1] = ph2id_new[j, 1] = i #
      ph2id_old[j,-1] = ph2id_new[j,-1] = sample(setdiff(atrisk_old, i), m)
      resamp = (ph2id_new[j,-1] %in% unique(as.numeric(na.omit(ph2id_new[-j,]))))
      repeat {
        if(any(resamp)) {
          ph2id_new[j, 1+which(resamp)] = sample(setdiff(atrisk_new, i), sum(resamp))
          resamptmp = (ph2id_new[j,1+which(resamp)] %in% unique(c(as.numeric(na.omit(ph2id_new[-j,])), #
                                                                  as.numeric(ph2id_new[j,-c(1+which(resamp))]))))
          resamp[which(resamp)] = resamptmp
        }
        else break
      }
    }
  }
  
  cohort_old$ind.ph2 = cohort_new$ind.ph2 = F
  cohort_old$ind.ph2[unique(c(ph2id_old))] = T
  cohort_new$ind.ph2[unique(c(ph2id_new))] = T
  
  #### inclusion probability / sampling weight ####
  inclpM = matrix(NA, nrow = nrow(cohort_old), ncol = m)
  colnames(inclpM) = paste("inclp.",1:m,sep="")
  for(mk in 1:m){
    inclpM[,mk] = with(cohort_old, mapply(function(i){
      1-prod(1-mk/(nrisk[which((eventime < eventime[i]) & (ind.fail))]-1))
    }, 1:nfull))
  }
  wgtM = apply(inclpM, 2, function(x){ifelse(cohort_old$ind.fail, 1, 1/x)})
  colnames(wgtM) = paste("wgt.",1:m,sep="")
  cohort_old = data.frame(cohort_old, inclpM, wgtM)
  
  inclpM = matrix(NA, nrow = nrow(cohort_new), ncol = m)
  colnames(inclpM) = paste("inclp.",1:m,sep="")
  for(mk in 1:m){
    inclpM[,mk] = with(cohort_new, mapply(function(i){
      1-prod(1-mk/(nriskM[which((eventime < eventime[i]) & (ind.fail)), mk]-1))
    }, 1:nfull))
  }
  wgtM = apply(inclpM, 2, function(x){ifelse(cohort_new$ind.fail, 1, 1/x)})
  colnames(wgtM) = paste("wgt.",1:m,sep="")
  cohort_new = data.frame(cohort_new, inclpM, wgtM)
  
  #### Data arrangement including duplicates ####
  ncc_old = cohort_old[c(ph2id_old),]
  ncc_old$matid = rep(1:sum(cohort_old$ind.fail), m + 1)
  ncc_old$matidm = rep(0:m, each = sum(cohort_old$ind.fail))
  cohort_old$matid = cohort_old$matidm = NA
  cohort_old = rbind(cohort_old[-unique(c(ph2id_old)),], ncc_old)
  cohort_old = cohort_old[order(as.numeric(rownames(cohort_old))),]
  cohort_old$ind.fail[which(as.numeric(rownames(cohort_old))%%1 != 0)] = F

  cohort_new = cbind(cohort_new, nriskM)
  ncc_new = cohort_new[c(ph2id_new),]
  ncc_new$matid = rep(1:sum(cohort_new$ind.fail), m + 1)
  ncc_new$matidm = rep(0:m, each = sum(cohort_new$ind.fail))
  cohort_new$matid = cohort_new$matidm = NA
  cohort_new = rbind(cohort_new[-unique(c(ph2id_new)),], ncc_new)
  cohort_new = cohort_new[order(as.numeric(rownames(cohort_new))),]

  return(list(old = cohort_old, new = cohort_new))
}