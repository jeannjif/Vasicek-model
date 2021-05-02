taux_ZC <- function(s = 5, Tz =10,N = 10, T = 5, n = 60){#price of zeros a zeros coupon bond following Vasicek model.
  #Tz = maturity of zeros coupon
  ##-----------Evaluate free_rate
  delta_t <- T/n
  eps <- matrix(rnorm(n*N), nrow = n, ncol = N)
  r <- matrix(NA, nrow = n+1, ncol = N)
  r[1,] <- 0.01
  a <- 0.2
  theta <- 0.04
  sig <- 0.015
  
  for(i in seq(2,n+1)){
    r[i,] <- r[i-1,] + a*(theta-r[i-1,])*delta_t + sig*sqrt(delta_t)*eps[i-1,]
    
  }
  ###-----evaluate zeros coupon-----
  t <-  (n/T) + 1 #bcoz r_0 is already given
  P_ts <- matrix(NA, nrow = T, ncol = N)
  B_ts <- (1-exp(-a*10))/a # move outside, also A_ts
  A_ts <- exp((theta - (sig^2)/(2*a^2))*(B_ts - 10) - ((sig^2)/(4*a))*B_ts^2)
  for(y in seq(1,s)){
    P_ts[y,] <- A_ts*exp(-B_ts*r[t,])
    t <- t+ (n/T)                     
  }
  new_list <- list(rate = r, taux_zc = -log(P_ts)/Tz) #access by doing: price_ZC$price or rate
  return(new_list)
  
}

#####################################################################""

actualisation <- function(s,rate){ #evaluate /beta with simpson. Please note that s is month (i.e the same step as rate)
  leng = nrow(rate)
  delta_x <- s/(leng-1)
  f_a <- rate[1,]
  f_b <- rate[leng,]
  f_odd <- 0
  f_pair <- 0
  for(i in seq(2,leng-1, 2)){ 
    f_odd <- f_odd + rate[i,] 
  }
  f_odd <- 4*f_odd
  for(i in seq(3,leng-1, 2)){
    f_pair <- f_pair + rate[i,]
  }
  f_pair <- 2*f_pair
  integr <- (delta_x/3) * (f_a + f_odd + f_pair + f_b)
  return(exp(-integr))
}

############################################################"""
floating_coupon <- function(coupon,rat,T = 5,n=60){#evaluate floating coupon bond by
                                                   #Monte carlo
  K_s <- matrix(NA, nrow = nrow(coupon), ncol = ncol(rat))
  
  #------------evaluate cashflow---------
  max_year <- nrow(coupon)
  for(s in seq(1, max_year-1)){
    mi <- pmin(0.06, coupon[s,])
    K_s[s,] <- 100*pmax(mi, 0.03)
  }
  last <- nrow(coupon)
  K_s[last,] <- 100*(pmax(pmin(coupon[last,],0.06),0.03) + 1)
  
  #------Now, do the monte_carlo simulation
  t <-  (n/T) + 1 #bcoz r_0 is already given
  over_beta <- matrix(NA,nrow = nrow(coupon) , ncol = ncol(rat))
  payoff <- matrix(NA,nrow = 1, ncol = ncol(rat))
  for(yea in seq(1,T)){
    over_beta[yea,] <- actualisation(yea, rat[1:t,])
    t <- t+ (n/T)
  }
  for(sim in seq(1,ncol(rat))){
    payoff[1,sim] <- sum(K_s[,sim]*over_beta[,sim])
  }
  price <- mean(payoff)*1
  S <- sd(payoff)
  #----building confidence interval: TODO: ask if it's good
  l <- price - qnorm(0.975)*S/sqrt(ncol = ncol(rat))
  h <- price + qnorm(0.975)*S/sqrt(ncol = ncol(rat))
  inter <- c(h,l)
  new_list <- list(Ks = K_s, pri = price, IC = inter)
  return(new_list)
}

##############partie B

prix_zc_mc <- function(rate, Tz = 5){#Compare price of zero coupon evaluate by monte carlo vs vasicek
                              # Tz = maturity of zero coupon in year
  t <- ((nrow(rate)-1)/Tz) +1
  if (Tz>1){
    for(i in seq(2,Tz)){
      t<- t+ ((nrow(rate)-1)/Tz)
    }
  }
    P_0T <- actualisation(Tz,rate[1:t,])
  ###-----evaluate zeros coupon-----
  a <- 0.2
  theta <- 0.04
  sig <- 0.015
  B_ts <- (1-exp(-a*Tz))/a # move outside, also A_ts
  A_ts <- exp((theta - (sig^2)/(2*a^2))*(B_ts - Tz) - ((sig^2)/(4*a))*B_ts^2)
  P_ts <- A_ts*exp(-B_ts*rate[1,])
  valu <- list(sim = mean(P_0T), val = P_ts[1])
  return(valu)
  
}


taux <- taux_ZC()
monte <- floating_coupon(taux$taux_zc, taux$rate)
matplot(taux$rate,type="l")
compar <- prix_zc_mc(taux$rate,5)
error_zc <- abs(compar$sim-compar$val)
monte$pri
error_zc

