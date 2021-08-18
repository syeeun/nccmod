
ncc.design = function(nsimul = DUMX, output1 = "DUMY1", seed = DUMZ, nfull, maxm = 5, ...){
  
  res = array(NA, dim=c(5,3,5,nsimul),
              dimnames = list(1:maxm,
                              c("beta1", "beta2", "Lambda0"),
                              c("full", "standard/partial", "modified/partial", "standard/pseudo", "modified/pseudo"),
                              1:nsimul))
  
  
  set.seed(seed)
  for(i in 1:nsimul){
    tryCatch({
      cohort = simul.cohort(n = nfull, max.t = 10, ...)
      cohort.ncc = sampling.ncc.maxm(cohort, m = maxm, tau = 10)
      for(k in 1:maxm){
        cohort.ncc$old.unified = subset(cohort.ncc$old, subset = (as.numeric(rownames(cohort.ncc$old)) %% 1 == 0))
        cohort.ncc$new.unified = subset(cohort.ncc$new, subset = (as.numeric(rownames(cohort.ncc$new)) %% 1 == 0))
        
        fit.full = coxph(Surv(eventime, ind.fail) ~ z1 + z2, data = cohort.ncc$old.unified)
        bhaz.full = bhaz.wgt(fit.full, wgt = 1)
        fit.LB.old = coxph(Surv(eventime, ind.fail) ~ z1 + z2 + strata(matid), data = cohort.ncc$old, subset = (ind.ph2 & (matidm <= k)))
        fit.LB.new = coxph(Surv(eventime, ind.fail) ~ z1 + z2 + strata(matid), data = cohort.ncc$new, subset = (ind.ph2 & (matidm <= k)))
        fit.Sam.old = coxph(Surv(eventime, ind.fail) ~ z1 + z2, data = cohort.ncc$old.unified, subset = (ind.ph2 & (matidm <= k)), weights = get(paste("wgt.",k,sep="")))
        fit.Sam.new = coxph(Surv(eventime, ind.fail) ~ z1 + z2, data = cohort.ncc$new.unified, subset = (ind.ph2 & (matidm <= k)), weights = get(paste("wgt.",k,sep="")))
        bhaz.LB.old = bhaz.wgt.LB(fit.LB.old, dat = cohort.ncc$old, m = k)
        bhaz.LB.new = bhaz.wgt.LB(fit.LB.new, dat = cohort.ncc$old, m = k)
        bhaz.Sam.old = bhaz.wgt(fit.Sam.old)
        bhaz.Sam.new = bhaz.wgt(fit.Sam.new)
        
        res[k,,,i] = cbind(full = c(fit.full$coefficients, sum(bhaz.full)),
                           LB.old = c(fit.LB.old$coefficients, sum(bhaz.LB.old)),
                           LB.new = c(fit.LB.new$coefficients, sum(bhaz.LB.new)),
                           Sam.old = c(fit.Sam.old$coefficients, sum(bhaz.Sam.old)),
                           Sam.new = c(fit.Sam.new$coefficients, sum(bhaz.Sam.new)))
        
        # cat("================ simul:", i, "; m=", mseq[k] ,"================", "\n")
      }
      # cat("================ simul:", i,"\n")
      cat(",", i)
    }, error = function(e) e)}
  
  cat("\n")
  saveRDS(res, file = output1)
}

