# Estimate inverse probability weighted baseline hazard function
bhaz.wgt = function(coxph_object, wgt = NULL){
  
  coxdetail_object = coxph.detail(coxph_object, riskmat = T)
  
  Z = model.matrix(coxph_object)
  Zbeta = Z %*% coxph_object$coefficients
  
  risk_mat = coxdetail_object$riskmat
  if(is.null(wgt)) wgt = coxph_object$weight
  YwexpbZ = sweep(risk_mat, 1, wgt*exp(Zbeta), `*`)
  
  dNt = matrix(0, nrow(risk_mat), ncol(risk_mat))
  evnt = apply(risk_mat, 1, function(x) { ifelse(all(x==0), 0, max(which(x == 1))) } )
  dNt[cbind(1:nrow(risk_mat), evnt)] = 1
  dNt = sweep(dNt, 1, coxph_object$y[,ncol(coxph_object$y)], `*`)
  
  Gt = colSums(YwexpbZ)
  H0t = colSums(dNt)/Gt
  
  return(H0t)
}

# Estimate Langholz-Borgan weighted baseline hazard function
bhaz.wgt.LB = function(coxph_object, dat, m){
  coxdetail_object = coxph.detail(coxph_object)
  
  Z = coxdetail_object$x
  failtime = coxdetail_object$y[,2]
  Zbeta = Z %*% coxph_object$coefficients
  
  caseid = which(1:nrow(Z) %% (m+1) == 1)
  if(!is.null(dat$nrisk)){
    wgt = dat[rownames(Z[caseid,]),]$nrisk/(m+1) 
  }
  if(is.null(dat$nrisk)){
    wgt = as.numeric(unlist(subset(dat[rownames(Z[caseid,]),], select = get(paste("nrisk.",1,sep="")))))/(m+1)
  }
  bhaz_vec = 1/rowSums(sweep(matrix(exp(Zbeta), ncol = m+1, byrow = T), 1, wgt, `*`))
  bhaz_vec = bhaz_vec[order(failtime[caseid])]
  names(bhaz_vec) = sort(failtime[caseid])
  
  return(bhaz_vec)
}

# Summarize Results from ncc.design()
summary.res = function(res){
  MEAN = apply(res, c(1,2), mean, na.rm = T)
  VAR = apply(res, c(1,2), var, na.rm = T)
  SD = apply(res, c(1,2), sd, na.rm = T)
  RE = VAR[,1]/VAR*100
  BIAS = abs((MEAN-MEAN[,1]))#/MEAN[,1])*100
  DV = VAR - VAR[,1]
  return(list(MEAN = MEAN, BIAS = BIAS, VAR = VAR, RE = RE, SD = SD, DV = DV))
}

