simul.cohort = function(n, max.t = 10, Smax.t = 0.95, 
                        alpha = 1, relz1 = 1.5, relz2 = 2.5, 
                        Smax.censor = 0.90, 
                        #max.enter = 2, #Uncomment for left truncation
                        mu = c(z0 = 0, z1 = 0, z2 = 0),
                        sig = c(sd0 = 1, sd1 = 1, sd2 = 1),
                        rho = c(rho01 = .05, rho02 = .24, rho12 = .25)){
  
  ##########################################
  logit = function(p) {log(p/(1-p))}
  inv.logit = function(x) {1/(1+exp(-x))}
  ##########################################
  
  Sig = matrix(NA, nrow = 3, ncol = 3)
  diag(Sig) = sig^2
  Sig[upper.tri(Sig)] = Sig[lower.tri(Sig)] = rho*c(sig[1]*sig[2], sig[1]*sig[3], sig[2]*sig[3])
  Z = data.frame(mvrnorm(n = n, mu = mu, Sigma = Sig))
  
  # coefficients
  beta0 = log(-log(Smax.t)/(max.t)^alpha); beta1 = log(relz1); beta2 = log(relz2)
  lam = with(Z, exp((beta0 + beta1*z1 + beta2*z2)/alpha))
  
  # survival data
  survtime = mapply(rweibull, n = 1, shape = alpha, scale = 1/lam)
  # entrtime = runif(n, 0, max.enter) #Uncomment for left truncation
  deadtime = rexp(n, rate = -log(Smax.censor)/max.t) 
  censtime = max.t #- entrtime #Uncomment for left truncation
  eventime = mapply(min, survtime, censtime, deadtime) 
  ind.fail = (survtime <= mapply(min, censtime, deadtime)) 
  
  
  DATA = data.frame(cbind(Z, eventime, ind.fail))
  DATA = DATA[order(DATA$eventime),]
  
  rownames(DATA) = 1:n
  return(DATA)
}
