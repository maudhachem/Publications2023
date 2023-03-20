# Filename: LK2speciesDistrDelayODEsover.R
# Author: Nick Beeton, CSIRO Hobart, March 2023
# Description: This file contains the code to solve equations (22)-(27)
# and displays the solutions N1(t) and N2(t) shown in Figure 2(a)-(c)
# from the following reference
# M El-Hachem and NJ Beeton (2023) 
# Coexistence in two-species competition with delayed maturation
T = 60

library(deSolve)
for (example in c("a","b","c"))
{
  m = c(0.3, 0.3)
  mu = c(0.3, 0.3)
  r = c(1,1)
  a = c(1,1)
  {
    if (example != "c") b = c(0.1,0.1) # a, b
    else b = c(3,3)# c
  }
  
  {
    if (example == "a") Delta = c(1.5922684937891114,12.738147950312891) # a
    else Delta = c(1.5922684937891114,3.184536987578223) # b, c
  }
  p = c(20,20)
  
  {
    if (example == "a") beta = p / c(1,8) # a
    else beta = p / c(1,2) # b, c
  }
  
  Ninit = matrix(c(0.5, 0.5), 2, 1)
  
  f = function(t, state, parameters){
    with(as.list(c(state, parameters)), { 
      N1 = state[1]
      N2 = state[2]
      x1 = state[2 + (1:p1)]
      x2 = state[2 + p1 + (1:p2)]
      f1 = state[3 + p1 + p2]
      f2 = state[4 + p1 + p2]
      y1 = state[4 + p1 + p2 + (1:p1)]
      y2 = state[4 + 2*p1 + p2 + (1:p2)]
      
      l = 0:(p1-1)
      expsum1 = sum(exp(l*log(beta1*Delta1) - lfactorial(l)))
      l = 0:(p2-1)
      expsum2 = sum(exp(l*log(beta2*Delta2) - lfactorial(l)))
      
      denom1 = 1 - exp(-beta1*Delta1) * expsum1
      denom2 = 1 - exp(-beta2*Delta2) * expsum2
      
      dN1dt = r1 * (beta1 / (beta1 + mu1))^p1 * x1[p1] * exp(-Delta1 * (a1 * f1 + b2 * f2)) - m1*N1
      dN2dt = r2 * (beta2 / (beta2 + mu2))^p2 * x2[p2] * exp(-Delta2 * (a2 * f2 + b1 * f1)) - m2*N2
      
      dx1dt = numeric(p1)
      dx1dt[1] = (mu1 + beta1) * (N1 - x1[1]) /denom1
      if (p1 > 1)
        for (k in 2:p1)
          dx1dt[k] = (mu1 + beta1) * (x1[k-1] - x1[k]) / denom1
      
      dx2dt = numeric(p2)
      dx2dt[1] = (mu2 + beta2) * (N2 - x2[1]) / denom2
      if (p2 > 1)
        for (k in 2:p2)
          dx2dt[k] = (mu2 + beta2) * (x2[k-1] - x2[k]) / denom2
      
      df1dt = (N1 - y1[p1]) / Delta1
      df2dt = (N2 - y2[p2]) / Delta2
      
      dy1dt = numeric(p1)
      dy1dt[1] = beta1 * (N1 - y1[1]) /denom1
      if (p1 > 1)
        for (k in 2:p1)
          dy1dt[k] = beta1 * (y1[k-1] - y1[k]) / denom1
      
      dy2dt = numeric(p2)
      dy2dt[1] = beta2 * (N2 - y2[1]) /denom2
      if (p2 > 1)
        for (k in 2:p2)
          dy2dt[k] = beta2 * (y2[k-1] - y2[k]) / denom2
      
      list(c(dN1dt, dN2dt, dx1dt, dx2dt, df1dt, df2dt, dy1dt, dy2dt))
    })
  }
  
  parameters = c(m1 = m[1], m2 = m[2],
                 mu1 = mu[1], mu2 = mu[2],
                 r1 = r[1], r2 = r[2],
                 a1 = a[1], a2 = a[2],
                 b1 = b[1], b2 = b[2],
                 Delta1 = Delta[1], p1 = p[1],
                 Delta2 = Delta[2], p2 = p[2],
                 beta1 = beta[1], beta2 = beta[2])
  
  N = Ninit
  state = c(N1 = N[1], N2 = N[2],
            x1 = rep(0,parameters[["p1"]]),
            x2 = rep(0,parameters[["p2"]]),
            f1 = 0, f2 = 0,
            y1 = rep(0,parameters[["p1"]]),
            y2 = rep(0,parameters[["p2"]]))
  
  times = seq(0.01, T, 0.01)
  
  out = lsoda(y = state, times = times, func = f, parms = parameters)
  out = as.data.frame(out)
  # Julia order: N1 N2 f1 f2 x1 y1 x2 y2
  # R order: N1 N2 x1 x2 f1 f2 y1 y2
  # ii = c(1:2, 5:24, 45:64, 3:4, 25:44, 65:84)
  # state = ii/84
  # t = 0
  # out = f(t, state, parameters)[[1]]
  # out[order(ii)]
  
  plot(out$time, out$N1, 
       xlim = c(0,T), ylim = c(0,1), 
       col = "limegreen", type = 'l', lwd = 2, 
       xlab = "Time", ylab = "Pop",cex.lab=2, cex.axis=2)
  lines(out$time, out$N2, col = "goldenrod2", lwd = 2)
  {
    if (example == "a") Delta = c(4.605170185988091,36.84136148790473) # a
    else Delta = c(4.605170185988091,9.210340371976182) # b, c
  }
  p = c(1,1)
  
  {
    if (example == "a") beta = p / c(1,8) # a
    else beta = p / c(1,2) # b, c
  }
  parameters = c(m1 = m[1], m2 = m[2],
                 mu1 = mu[1], mu2 = mu[2],
                 r1 = r[1], r2 = r[2],
                 a1 = a[1], a2 = a[2],
                 b1 = b[1], b2 = b[2],
                 Delta1 = Delta[1], p1 = p[1],
                 Delta2 = Delta[2], p2 = p[2],
                 beta1 = beta[1], beta2 = beta[2])
  
  N = Ninit
  state = c(N1 = N[1], N2 = N[2],
            x1 = rep(0,parameters[["p1"]]),
            x2 = rep(0,parameters[["p2"]]),
            f1 = 0, f2 = 0,
            y1 = rep(0,parameters[["p1"]]),
            y2 = rep(0,parameters[["p2"]]))
  
  times = seq(0.01, T, 0.01)
  
  out = lsoda(y = state, times = times, func = f, parms = parameters)
  out = as.data.frame(out)
  #out = pmin(pmax(out, -1e5), 1e5)
  lines(out$time, out$N1, col = "dodgerblue2", lwd = 2)
  lines(out$time, out$N2, col = "red2", lwd = 2)
  

}

