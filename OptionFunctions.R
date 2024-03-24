MCVD = function(S = 1, K = 1, b = 0.9, tau = 1, r = 0.015, q = 0.01, sigma = 0.13, N = 252, nSim = 1000,
                guyseed = 1930, flag = "fixed"){
  
# Input parameters (function arguments)

# S - stock price at time t
# K - strike price  
# b - barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
# tau - time to maturity in years 
# r - risk-free annual interest rate 
# q - deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
# sigma - annual volatility std dev) of the notional price (NB NOT scaleable by sqr(T) to observed price) 
# N -  number of time steps
# nSim - 1000  #number of simulations (paths) 
# guyseed - the seed fed in from main run
# flag = "fixed" for a fixed barrier of b, or "trailing" for a trailing barrier of b * HWM of Y (quotes are required)
  
set.seed(guyseed)
  
# "Fixed" parameters
  
c <- -0.0001 # to investigate limit on trading close to barrier: cannot trade when Y - b < c. 
#It turns out c has only marginal effect
#If don't want limit, set c to tiny neg number is safest (tiny discrepancies otherwise, prob from "Y = b exactly" cases).  
h <- 1 # hedging frequency - hedge every h time steps. Normally I just use 1.
epsilon <- 0.0001 # variation in X for finite difference calc of delta

#Check validity of inputs 
stopifnot(b <= min(S,K), r!=q, K!=b)

#analytic prices
# z's as in the paper
z1 <- (log(S/K) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z2 <- (log(b^2/(K*S)) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z3 <- (log(S/b) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z4 <- (log(b/S) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
theta <- 2*(r-q)/sigma^2

BS_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau)) #BS value call
BS_put <- -S*exp(-q*tau)*pnorm(-z1) + K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) #BS value put

Analytic_barrier_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau))+
  1/theta*(S*exp(-q*tau)*(b/S)^(1+theta)*pnorm(z2) - K*exp(-r*(tau))*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))

Analytic_barrier_put <- K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) - S*exp(-q*tau)*pnorm(-z1)-
  b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) + S*exp(-q*tau)*pnorm(-z3)+
  1/theta*(b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) - S*exp(-q*tau)*(b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)) - K*exp(-r*tau)*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))

Analytic_barrier_put_delta <- exp(-q*tau)*(pnorm(z1) - pnorm(z3) + (b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)))

# Compute the Monte Carlo price for option with trailing barrier - using temp (tempoorary) variables 

dt_temp_0 <- tau/N #length of each time sub interval
Z_temp_0 <- array(rnorm(nSim*N, mean=0, sd=1), dim=c(nSim, N))
#generating fresh noise here - appropriately different from what you're trying to hedge, 
#but same for both valuations in the finite difference calculation (p1274 P Wilmott on Quant Finance)
dW_temp_0 <- Z_temp_0 * sqrt(dt_temp_0) #Brownian motion increments (nSim_temp simulations) x (N_temp increment
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_trailing <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_T <- vector(length=nSim) 

Option_payoff_bought_put_temp_0 <-array(dim=c(nSim))
Monte_Carlo_put_price_temp_0 <- array(dim=c(nSim, 2))


for (n in 1:2){

  X[,1] <- S  + (n-1) * epsilon 
  Y[,1] <- X[,1]
  b_as_matrix[,1] <- b 
  b_trailing[,1] <- b
  B_running_max[,1] <- b_trailing[,1]/X[,1]
  Y_running_max[,1] <-  Y[,1]

    for(j in 2:N){
  
      if (flag == "fixed"){
        #This bit for fixed barrier
      X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt_temp_0 + sigma*dW_temp_0[,j])
      B[,j] <- b_as_matrix[,1]/X[,j] #b never changes - doesn't really need to be a matrix (but does when we come to b-trailing)
      B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
      Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
      
        }else{
  
      #This bit for trailing barrier
      X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt_temp_0 + sigma*dW_temp_0[,j])
      B[,j] <- b_trailing[,j-1]/X[,j] #Have to use b_trailing[j-1] because don't know new b_trailing before calculating latest Y
      B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
      Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
      Y_running_max[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], Y[,j], Y_running_max[,j-1])
      b_trailing[,j] <- b * Y_running_max[,j] # ifelse(Y[,j] > Y_running_max[,j-1], b * Y[,j], b_trailing[,j-1])
      }
      
    } # end of (j:N) steps of j, 
  
  Option_payoff_bought_put_temp_0 <- pmax(K - Y[,N],0) 
  
  Monte_Carlo_put_price_temp_0[n] <- exp(-r*tau)*mean(Option_payoff_bought_put_temp_0)
  
  
} #end of(1 to 2) steps of n, for FD_delta calc



Thomas_delta_bought_put_by_FD_0 <- (Monte_Carlo_put_price_temp_0[2] - Monte_Carlo_put_price_temp_0[1])/epsilon

#That gives delta wrt to X. At time 0, X = Y, so it's also delta wrt Y.


return(c(Monte_Carlo_put_price_temp_0[1], Thomas_delta_bought_put_by_FD_0))

} # end function
    




