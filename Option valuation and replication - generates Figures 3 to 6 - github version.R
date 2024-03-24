#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purposes.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at  R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

source("OptionFunctions.R")


p <- 4  # Setting number of instances of k simulation. Was 9 in old paper. Can prob only do 5 here.


St_dev_hedging_errors_put_by_FD <- matrix(nrow=0, ncol=2) 
St_dev_hedging_errors_put_analytic <- matrix(nrow=0, ncol=2)

for(m in 4:p){  
  
  
# Setting parameters 

barrier_type <- "fixed"    # "fixed" or "trailing".  NEEDS QUOTES - defines a string def to MCVD() function.

S <- 1.0 # Asset price at time 0
K <- 1.0  # strike price 
b <- 0.8 # barrier 
tau <- 5  #time to maturity (in years) 
r <- 0.015 #risk-free annual continuous rate of interest

mu <- 0.03 #real-world asset drift, used to simulate the payoff we have to replicate. But need to change r to mu
               # in the main ( in 2:N) loop below to use this. Normally I leave it as r - that makes available  
               # a larger number of drift-neutralised simulations to compute the MC option price.

q <- 0.01 #deferment rate (yield) ta)
sigma <- 0.13 #annual volatility of the asset price (standard deviation)

c <-  -0.001 # to investigate limit on trading close to barrier: cannot trade when Y - b_trailing < c. 
#It turns out c has only marginal effect
#If don't want limit, set c to tiny neg number is safe (tiny discrepancies otherwise, prob from "Y = b_trailing exactly" cases).  
h <- 1 # hedging frequency - hedge every h time steps. Normally I just use 1.
epsilon <- 0.0001 # finite difference. Use 4th root of machine precision x initial stock price (normally 1) = 0.00012220703,  
                                    #  using .Machine command in R to find machine precision (Jackel 2002 book, p138).

guyseed <- 1930 # set like this so can send it to external function
set.seed(guyseed) #set the seed (if not set, it's taken from the computer clock)

N <-   60 # * 2^(m-1) # nuber of time steps. ^2(m-1) uses more time steps in each instance, to plot convergence
FD_delta_N_multiplier <- 2 # multiple of N for number of time steps in (longest) FD delta calculations

nSim <- 100 #number of simulations (paths) #)
FD_delta_nSim_multiplier <- 1  #multiple of nSim for number nSim_temp of simulations in FD delta calcs 


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

Analytic_forward <- S * exp(-q * tau) - K * exp(-r *tau)


# Compute the Monte Carlo prices 

dt <- tau/N #length of each time sub interval
Z <-  matrix(rnorm(nSim*N, mean=0, sd=1),nrow = nSim, ncol = N) #standard normal sample of N elements
dW <- Z*sqrt(dt) #Brownian motion increments (nSim simulations) x (N increments)
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_trailing <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_T <- vector(length=nSim) 
Intrinsic_value_of_call <-matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Intrinsic_value_of_put <-matrix(numeric(nSim*N), nrow = nSim, ncol = N)



X[,1] <- S 
Y[,1] <- X[,1]
b_as_matrix[,] <- b
b_trailing[,1] <- b
B_running_max[,1] <- b_trailing[,1]/X[,1]
Y_running_max[,1] <- S

Option_payoff_bought_put <- numeric(nSim)


Thomas_final_cash <-numeric(nSim)
Thomas_final_stock <-numeric(nSim)
Thomas_hedging_error <- numeric(nSim)

Thomas_call_final_cash <-numeric(nSim)
Thomas_call_final_stock <-numeric(nSim)
Thomas_call_hedging_error <- numeric(nSim)


#Black-Scholes replication of barrier put    

BS_delta_bought_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_delta_bought_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_before_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

BS_delta_bought_put[,1] <- exp(-q*tau)*(pnorm(z1)-1) #negative number
BS_stock_position_after_trade[,1] <- BS_delta_bought_put[1] 
BS_cash_account[,1] <-  -S * BS_delta_bought_put[1] + BS_put 
BS_delta_bought_put_most_recent_permitted_trade[,1] <-BS_delta_bought_put[,1] 
#opening the position always permitted
BS_replicating_portfolio[,1] <- BS_stock_position_after_trade[,1] * Y[,1] + BS_cash_account[,1]
BS_trade_size[,1] <- 0
  
#Thomas put replication
Thomas_delta_bought_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Thomas_synthetic_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Difference_put_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Option_payoff_bought_put_temp <-array(dim=c(nSim, N * FD_delta_N_multiplier, nSim * FD_delta_nSim_multiplier,2)) # add extra 2 dimensions for debugging
Thomas_MC_put_price_temp <- array(dim=c(nSim, N, 2)) # extra 2-dimension because need 2 estimates to get FD delta
Thomas_delta_bought_put_by_FD <- array(dim=c(nSim, N)) 

Thomas_delta_bought_put_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_put_most_recent_permitted_trade_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_stock_position_after_trade_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_trade_size_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_cash_account_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_replicating_portfolio_analytic <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)


Thomas_z1 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z2 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z3 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z4 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
# Now call function for the initial MC valuation and FD delta

Initial_values <- MCVD(S, K, b, tau, r, q, sigma, N, nSim * 100, guyseed, flag = barrier_type) # Use nSim *100 & N*10 for accurate initial price & delta

# May prefer to overrtide with Initial_values from a long run of MCVD() in OptionFunctions.R - specimen values below.

# for 10 yrs, b = 0.6, sigma = 0.20, K = 1.0
#Initial_values <- if(barrier_type == "fixed"){ c(0.106891364, - 0.100434506) } else {
#  c(0.04853437, -0.24549251)} 

# for 10 yrs, b = 0.6, sigma = 0.20, K = 1.25
# if(barrier_type == "fixed"){ Initial_values <- c(0.241158256, -0.198679881) } else {
#                     Initial_values <-   c(0.1473610, -0.4620159)}  

# for tau = 5, b =0.6, sigma = 0.20 
#if(barrier_type == "fixed") { Initial_values <- c(0.115606759, -0.21091623) } else {
#                     Initial_values <-   c( 0.08074558, -0.37161128)}  

# for tau = 5, b =0.6, sigma = 0.13 
#  if(barrier_type == "fixed") { Initial_values <- c(0.09138134, -0.3397917) } else {
#                     Initial_values <- (0.08484505, -0.39559939)}  
    



Thomas_MC_put_price_temp[,1,1] <-  Initial_values[1] 
Thomas_delta_bought_put_by_FD[,1] <- Initial_values[2] 

Thomas_stock_position_after_trade[,1] <- Thomas_delta_bought_put_by_FD[,1] # same as above
Thomas_cash_account[,1] <- - S * Thomas_delta_bought_put_by_FD[,1] + Thomas_MC_put_price_temp[1] 
#first term is positive, becasue for a bought put we go short the asset. 
Thomas_replicating_portfolio[,1] <- Thomas_stock_position_after_trade[,1] * Y[,1] + 
                                   Thomas_cash_account[,1]
Thomas_delta_bought_put_most_recent_permitted_trade[,1] <-Thomas_delta_bought_put_by_FD[,1] 
#opening the position always permitted
Thomas_trade_size[,1] <- 0

Thomas_delta_bought_put_analytic[,1] <- exp(-q*tau)*(pnorm(z1) - pnorm(z3) + (b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)))
Thomas_stock_position_after_trade_analytic[,1] <- Thomas_delta_bought_put_analytic[,1] # same as above
Thomas_cash_account_analytic[,1] <- - S * Thomas_delta_bought_put_analytic[,1] + Analytic_barrier_put 
Thomas_replicating_portfolio_analytic[,1] <- Thomas_stock_position_after_trade_analytic[,1] * Y[,1] + 
  Thomas_cash_account_analytic[,1]
Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,1] <-Thomas_delta_bought_put_analytic[,1] 
Thomas_trade_size_analytic[,1] <- 0


Thomas_synthetic_call_replicating_portfolio[,1] <-  Thomas_stock_position_after_trade[,1] *Y[,1]+
  Thomas_cash_account[,1] +
  Y[,1]*exp(-q*tau) - K*exp(-r*tau)


#Thomas call replication

Thomas_delta_bought_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_call_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Difference_call_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

#Forward F' and synthetic forward F, replication

Forward_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Forward_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Synthetic_forward_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Synthetic_forward_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Difference_forward_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_bought_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_call_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_synthetic_put_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Interim_rolled_up_difference_in_analytic_call_prices <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <- 
  matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <- 
  matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_bought_call[,1] <- exp(-q*tau)*(pnorm(z1) - (b/S)^(1+theta)*pnorm(z2))
#using z1, z2 from top of program...OK at time 1.
Thomas_call_stock_position_after_trade[,1] <- Thomas_delta_bought_call[,1] # same as above
Thomas_call_cash_account[,1] <- -S * Thomas_delta_bought_call[,1] + Analytic_barrier_call 
#Negative number, as for BS, for low barrier. Can be positive for high barrier (when
#the delta of the Thomas call gets low).
Thomas_call_replicating_portfolio[,1] <- Thomas_call_stock_position_after_trade[,1] * Y[,1] + 
  Thomas_call_cash_account[,1]
Thomas_delta_bought_call_most_recent_permitted_trade[,1] <- Thomas_delta_bought_call[,1] 
#opening the position always permitted
Thomas_call_trade_size[,1] <- 0

Thomas_synthetic_put_replicating_portfolio[,1] <-  Thomas_call_replicating_portfolio[,1] -
                                                     (Y[,1]*exp(-q*tau) - K*exp(-r*tau))

Difference_put_replicating_portfolios[,1] <- Thomas_synthetic_put_replicating_portfolio[,1] - 
                                                Thomas_replicating_portfolio[,1]  

Interim_rolled_up_difference_in_analytic_call_prices[,1] <- Analytic_barrier_call - 
                                                           (Analytic_forward + Analytic_barrier_put)

Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices[,1] <-
  Thomas_synthetic_call_replicating_portfolio[,1] /  Interim_rolled_up_difference_in_analytic_call_prices[,1]


#Forward F and synthetic forward F*, replication
Forward_replicating_portfolio[,1] <- Y[,1]* exp(-q * tau) - K * exp(-r*tau)
Synthetic_forward_replicating_portfolio[,1] <- Thomas_call_replicating_portfolio[,1] -
                 Thomas_replicating_portfolio[,1]
Difference_forward_replicating_portfolios[,1] <- Synthetic_forward_replicating_portfolio[,1] -
                                                     Forward_replicating_portfolio[,1]
                                   

#=============#


#Then do the loop - with hedging every h time steps


for(j in 2:N){
    
  if (barrier_type == "trailing"){
    X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt + sigma*dW[,j])
    B[,j] <- b_trailing[,j-1]/X[,j] #Have to use b_trailing[j-1] because don't know new b_trailing before calculating latest Y
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    Y_running_max[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], Y[,j], Y_running_max[,j-1])
    b_trailing[,j] <- b * Y_running_max[,j]  #ifelse(Y[,j] > Y_running_max[,j-1], b * Y[,j], b_trailing[,j-1])
   
    }else{
   
  #Bring this next bit into action instead, for fixed barrier

    X[,j] <- X[,j-1]*exp((r - q -0.5 * sigma^2)*dt + sigma*dW[,j])
    B[,j] <- b_as_matrix[,j]/X[,j]
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    
    # Then set all elements of b_trailing[,] to b
    # Needed because b_trailing[,j] is the barrier variable used in hedging calcs. 
    # Easier to force this to b than change the variable.
    
    b_trailing[,j] <- b
    
  }
  
    #==============# 
    # Do the FD delta estimation here (can't easily define functions to apply to vectors, it seems.
    # maybe the "apply, 1" function would help)
    
    tau_temp <- tau * (1-j/N) # time steps ahead on outer loop basis - we only do the FD_delta calcualtion exercise for (j in 2:N) times.
    N_temp <- FD_delta_N_multiplier * (N + 2 - j) # For multiplier = 1, starts at N; when hedging final step N is calculated, N_temp = 2
    # correct, that's counting the (N-1)th step you start at and N for the end step you go to.
    nSim_temp <- FD_delta_nSim_multiplier * nSim 
    dt_temp <-  tau_temp/N_temp # starts as before; when tau_temp ticks down by equiv of 1 step, N_temp also ticks down. so dt_temp constant.
    Z_temp <- array(rnorm(nSim*N_temp*nSim_temp, mean=0, sd=1), dim=c(nSim, N_temp, nSim_temp)) 
    #generating fresh noise here - appropriately different from what you're trying to hedge, 
    #but same for both valuations in the finite difference calculation (p1274 P Wilmott on Quant Finance)
    dW_temp <- Z_temp * sqrt(dt_temp) #Brownian motion increments (nSim_temp simulations) x (N_temp increments)
    
    # EXTRA DIMENSION, n = 1 & 2, IN ALL THE FOLLOWING FOR DEBUGGING
    X_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))
    b_temp <- array(dim=c(nSim,N_temp, nSim_temp,2)) # that's the *initial" b in the MC first-differences simulation starting at step j 
    b_trailing_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))
    b_as_matrix_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))
    B_temp <- array(dim=c(nSim, N_temp,nSim_temp,2))
    B_running_max_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))
    Y_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))
    Y_running_max_temp <- array(dim=c(nSim, N_temp, nSim_temp,2))  
   
    
    #Now MC valns 1 and 2 for finite difference calcs ("temp" variables all relate to this). 
    #I tried vectorisng this loop by using the 4th dimension (n) of the array for 
    # n of 1 and 2, but surprisingly that increased runtime of whole program by 30%!)
    
    #browser()
    
    
     for (n in 1:2){
    
       #browser(expr = {j==117 && n==1})
       
       X_temp[,1,,n] <- rep(X[1:nSim,j] + (n-1) * epsilon, length.out = nSim_temp)
       #Repeating the nSim rows of X[,j] to fill out the (possibly larger) nSim_temp desired number of FD simulations.
        #epsilon is the finite difference, set at top of page. 
        Y_temp[,1,,n] <- X_temp[,1,,n] * pmax(1,B_running_max[,j]) 
        b_as_matrix_temp [,,,n] <- b # for "fixed" barrier case
        b_trailing_temp[,1,,n] <- rep(b_trailing[1:nSim,j] + ifelse(barrier_type == "trailing", 
                                        ifelse(b_trailing[1:nSim,j] > b, b * (n-1) * epsilon * Y_running_max[,j]/X[,j], 0),0), 
                                        length.out = nSim_temp) 
                                        
        
        # That's OK as value 1 in inner loop for "fixed", because set all b_trailing[,] == b for "fixed", in the outer loop above.
      
       # Reason for the ifelse addition of (n-1) * epsilon bit:
      
        # We're importing the corresponding b_trailing[,j] from the outer loop, to use as the starting value in the inner loop. 
        # The position of b_trailing[,j] is b * Y_running_max[,j], and hence depends on the whole history of X[,j[. So the 
        # finite difference path in which X[,j] is epsilon higher almo implies a higher b_trailing[,j]. 
        # Similar arguments could possibly be made for B_running_max_temp and Y_running_max_temp below (but it's less important 
        # there, so don't bother). But we need to get the barrier right, because otherwise, deep bug: increase in B_temp 
        # (and hence B_running-max_temp), as X_temp progresses further below b_trailing_temp, exactly (to machine precision)
        # offsets the finite difference in X_temp[n=1&2]. This can lead to Y_temp[1&2] identical. Hence 
        # Option_payoff_bought_put-temp[1&2} the same; their means Thomas_MC_put_price_temp[1&2] the same; and so delta 
        # wrongly estimated as 0, even as Y_temp goes some distance above the barrier. To see the problem, remove the epsilon 
        # here and then plot delta for some cases with big hedging errors - you'll see the delta gets stuck at 0 after the 
        # first barrier touch.
         
      
        B_running_max_temp[,1,,n] <- rep(B_running_max[1:nSim,j], length.out = nSim_temp)
        # Take the above from the "outside environment", so get Y at correct starting place
        Y_running_max_temp[,1,,n] <- rep(Y_running_max[1:nSim,j], length.out = nSim_temp)  # that's the lookback max of Y at step j 
       
         for (k in 2:N_temp){
          
           if (barrier_type == "trailing"){          
             X_temp[,k,,n] <- X_temp[,k-1,,n]*exp((r - q -0.5 * sigma^2)*dt_temp + sigma*dW_temp[,k,])
             B_temp[,k,,n] <- b_trailing_temp[,k-1,,n]/X_temp[,k,,n] #Have to use b_trailing[,k-1,] because don't know new b_trailing before calculating latest Y
             B_running_max_temp[,k,,n] <- ifelse(B_temp[,k,,n] > B_running_max_temp[,k-1,,n], B_temp[,k,,n], B_running_max_temp[,k-1,,n])
             Y_temp[,k,,n] <- X_temp[,k,,n]*pmax(1,B_running_max_temp[,k,,n])
             Y_running_max_temp[,k,,n] <- ifelse(Y_temp[,k,,n] > Y_running_max_temp[,k-1,,n], Y_temp[,k,,n], Y_running_max_temp[,k-1,,n])
             b_trailing_temp[,k,,n] <- b * Y_running_max_temp[,k,,n] # The original b is still the relevant multiplier here

             }else{

           #Bring this next bit into action instead, for fixed barrier

           X_temp[,k,,n] <- X_temp[,k-1,,n]*exp((r - q -0.5 * sigma^2)*dt_temp + sigma*dW_temp[,k,])
           B_temp[,k,,n] <- b_as_matrix_temp[,k,,n]/X_temp[,k,,n] 
           B_running_max_temp[,k,,n] <- ifelse(B_temp[,k,,n] > B_running_max_temp[,k-1,,n], B_temp[,k,,n], B_running_max_temp[,k-1,,n])
           Y_temp[,k,,n] <- X_temp[,k,,n]*pmax(1,B_running_max_temp[,k,,n])

             } 
           
        } # end of (2 to N_temp) steps of k, obtaining put value to feed into FD_delta first diff calc
       
       #browser(expr = {j==73 && n==2})
       
       Option_payoff_bought_put_temp[,j,,n] <- pmax(K - Y_temp[,N_temp,,n],0) # pmax is "parallel maximum" of 2 vectors.
    
       Thomas_MC_put_price_temp[,j,n] <- exp(-r*tau_temp)*apply(Option_payoff_bought_put_temp[,j,,n], 1, mean)
       
       
       } #end of(1 to 2) steps of n, for FD_delta calc
    
                                                                                                                                  
     Thomas_delta_bought_put_by_FD[,j] <- pmin(pmax(((Thomas_MC_put_price_temp[,j,2] - Thomas_MC_put_price_temp[,j,1])/epsilon) * X[,j]/Y[,j] , -1) ,0)  
       
     
     # * X[,j /Y[,j]] is rating down because we've just calculated delta wrt X (it was X[1:nSim,j] that we added epsilon to), 
     # but want wrt Y for hedging
     # pmin,pmax restricts delta to (-1,0) to discard bad FD estimates
     # ifelse(condition, Do 1, Do 2) operates on every element of a vector if(condition){Do 1}{Do 2 (opptional)} doesn't. 
     
     
     
     #==============#
          
    #Then after calculating the FD delta, do the hedging loop - First for BS hedging of the barrier put
   
      # If(j %% h > 0) branch (j modulo h) says: if we're not at a hedging day, 
      #  just maintain the existing cash and stock positions.
      
    
      if(j %% h > 0) {
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt 
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1]
        BS_delta_bought_put_most_recent_permitted_trade[,j] <- BS_delta_bought_put_most_recent_permitted_trade[,j-1]
        
      }else{
        
        BS_delta_bought_put[,j] <- exp(-q*tau*(1-j/N))*(pnorm((log(Y[,j]/K) + 
                                     (r- q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N))))-1)
        
        BS_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, BS_delta_bought_put[,j] - 
                                      BS_delta_bought_put_most_recent_permitted_trade[,j-1])
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1] + BS_trade_size[,j] 
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + 
                                 BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt - BS_trade_size[,j] * Y[,j] 
        BS_delta_bought_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
                                                                BS_delta_bought_put_most_recent_permitted_trade[,j-1],
                                                                 BS_delta_bought_put[,j]) 
        BS_replicating_portfolio[,j] <- BS_stock_position_after_trade[,j] *Y[,j] + BS_cash_account[,j]
        
      }    


    #Then for Thomas put hedging
    
    Thomas_z1[,j] <-(log(Y[,j]/K) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z2[,j] <-(log(b^2/(K*Y[,j])) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z3[,j] <-(log(Y[,j]/b) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z4[,j] <-(log(b/Y[,j]) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
 
    
    if(j %% h > 0) {
      Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + 
                                   Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1]
      Thomas_delta_bought_put_most_recent_permitted_trade[,j] <- Thomas_delta_bought_put_most_recent_permitted_trade[,j-1]
    
      Thomas_cash_account_analytic[,j] <- Thomas_cash_account_analytic[,j-1] * exp(r*dt) + 
        Thomas_stock_position_after_trade_analytic[,j-1] *Y[,j] * q *dt
      Thomas_stock_position_after_trade_analytic[,j] <- Thomas_stock_position_after_trade_analytic[,j-1]
      Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,j] <- Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,j-1]
    
        
       }else{
    
         
           if (j!=N) { #AMENDED HERE TO USE FD_DELTA
            Thomas_delta_bought_put[,j]<- Thomas_delta_bought_put_by_FD[,j] 
            Thomas_delta_bought_put_analytic[,j] <- exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) - pnorm(Thomas_z3[,j]) +
             (b/Y[,j])^(1+theta) * (pnorm(Thomas_z4[,j]) - pnorm(Thomas_z2[,j])))
           } else {
           Thomas_delta_bought_put[,j] <- ifelse(K > Y[,j], -1, 0)
           Thomas_delta_bought_put_by_FD[,j] <- Thomas_delta_bought_put[,j] # set that one too, to make graphs of FD & analytic deltas consistent
           Thomas_delta_bought_put_analytic[,j] <- ifelse(K > Y[,j], -1, 0)
            }
    
       
       Thomas_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, Thomas_delta_bought_put[,j] - 
                                         Thomas_delta_bought_put_most_recent_permitted_trade[,j-1])
       Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1] + Thomas_trade_size[,j]
       Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + 
                                     Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                                     Thomas_trade_size[,j] * Y[,j]
    
       Thomas_delta_bought_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
             Thomas_delta_bought_put_most_recent_permitted_trade[,j-1], Thomas_delta_bought_put[,j])
    
       
       Thomas_trade_size_analytic[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, Thomas_delta_bought_put_analytic[,j] - 
                                         Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,j-1])
       Thomas_stock_position_after_trade_analytic[,j] <- Thomas_stock_position_after_trade_analytic[,j-1] + Thomas_trade_size_analytic[,j]
       Thomas_cash_account_analytic[,j] <- Thomas_cash_account_analytic[,j-1] * exp(r*dt) + 
         Thomas_stock_position_after_trade_analytic[,j-1] *Y[,j] * q *dt -
         Thomas_trade_size_analytic[,j] * Y[,j]
       
       
       Thomas_delta_bought_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
                               Thomas_delta_bought_put_most_recent_permitted_trade[,j-1], Thomas_delta_bought_put[,j])

       Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
                               Thomas_delta_bought_put_most_recent_permitted_trade_analytic[,j-1], Thomas_delta_bought_put_analytic[,j])
       
              
       #Now set the delta manually if at last step
       if (j!=N) {
         Thomas_replicating_portfolio[,j] <- Thomas_stock_position_after_trade[,j] *Y[,j] +
           Thomas_cash_account[,j]
         Thomas_replicating_portfolio_analytic[,j] <- Thomas_stock_position_after_trade_analytic[,j] *Y[,j] +
           Thomas_cash_account_analytic[,j]
       }else{
         Thomas_replicating_portfolio[,j] <-  Thomas_delta_bought_put[,j]*Y[,j] +
           Thomas_cash_account[,j]
         Thomas_replicating_portfolio_analytic[,j] <-  Thomas_delta_bought_put_analytic[,j]*Y[,j] +
           Thomas_cash_account_analytic[,j]
         
       }
       
    # Do ysnthetic call, now you have the Thomas put value   
    if (j!=N) {
      Thomas_synthetic_call_replicating_portfolio[,j] <- Thomas_stock_position_after_trade[,j]*Y[,j] +
        Thomas_cash_account[,j] + Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
    }else{
      Thomas_synthetic_call_replicating_portfolio[,j] <-  Thomas_delta_bought_put[,j]*Y[,j] +
        Thomas_cash_account[,j] + Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
      
    }
   
    
    Interim_rolled_up_difference_in_analytic_call_prices[,j] <- 
         (Analytic_barrier_call - (Analytic_forward + Analytic_barrier_put))*exp(r * tau * j/N)
    
    Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices[,j] <-
        Thomas_synthetic_call_replicating_portfolio[,j] /  Interim_rolled_up_difference_in_analytic_call_prices[,j]
    
    Intrinsic_value_of_call[,j] <- pmax(Y[,j]- K,0)
    Intrinsic_value_of_put[,j] <- pmax(K - Y[,j],0)
    
    Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call[,j] <-
      Thomas_synthetic_call_replicating_portfolio[,j] - Intrinsic_value_of_call[,j]
      
    }
    
     #Then for Thomas call hedging 
    
    # Already calculated the required column of z_1[,j] and column of z_2[,j] above

    if(j %% h > 0) {
      Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1]
      Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- Thomas_delta_bought_call_most_recent_permitted_trade[,j-1]    
     
       }else{
      
      
          if (j!=N) {
          Thomas_delta_bought_call[,j] <- exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) -
                                                             (b/Y[,j])^(1+theta) * pnorm(Thomas_z2[,j]))
           } else {
          Thomas_delta_bought_call[,j] <- ifelse(K < Y[,j], 1, 0)
          }
    # Setting the final delta manually avoids the deep bug of z2, z3, z4 evaluating as 0/0 if Y_T[i,j]] = b 
    #   and time remaining = 0. 
    
       
         Thomas_call_trade_size[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 0, Thomas_delta_bought_call[,j] - Thomas_delta_bought_call_most_recent_permitted_trade[,j-1])
         Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1] + Thomas_call_trade_size[,j]
         Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                                          Thomas_call_trade_size[,j] * Y[,j]
         
         Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b_trailing[,j] < c, 
                 Thomas_delta_bought_call_most_recent_permitted_trade[,j-1], Thomas_delta_bought_call[,j])        
         
         #Now set the delta manually if at last step
         if (j!=N) {
         Thomas_call_replicating_portfolio[,j] <- Thomas_call_stock_position_after_trade[,j] *Y[,j] +
                                                      Thomas_call_cash_account[,j]
         }else{
          Thomas_call_replicating_portfolio[,j] <-  Thomas_delta_bought_call[,j]*Y[,j] +
                                                       Thomas_call_cash_account[,j]
          }
         
         if (j!=N) {
           Thomas_synthetic_put_replicating_portfolio[,j] <- Thomas_call_replicating_portfolio[,j] -
                                                               (Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N)))
         }else{
           Thomas_synthetic_put_replicating_portfolio[,j] <-  Thomas_delta_bought_call[,j]*Y[,j] +
             Thomas_call_cash_account[,j] - (Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N)))
         } 
         #Synthetic put has to go here, where call calc has just been done
         
         Difference_put_replicating_portfolios[,j] <- Thomas_synthetic_put_replicating_portfolio[,j] - 
                                                         Thomas_replicating_portfolio[,j] 
         
         Difference_call_replicating_portfolios[,j] <- Thomas_call_replicating_portfolio[,j] - 
                                                         Thomas_synthetic_call_replicating_portfolio[,j]  
         
         
         Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call[,j] <-
           Thomas_call_replicating_portfolio[,j] - Intrinsic_value_of_call[,j]
         
         
         
       }     
    
    
    #Forward F and synthetic forward F*, replication
    
    Forward_replicating_portfolio[,j] <- Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
    Synthetic_forward_replicating_portfolio[,j] <-  Thomas_call_replicating_portfolio[,j] -
      Thomas_replicating_portfolio[,j]
    
    Forward_replicating_portfolio_log_change[,j-1] <- log(Y[,j]) - log(Y[,j-1])
    # only the change in the stock position counts as volatility, not the flow of yield or interest
    Synthetic_forward_replicating_portfolio_log_change[,j-1] <- 
      log(max(abs(Thomas_call_stock_position_after_trade[,j] - Thomas_stock_position_after_trade[,j]) *Y[,j], 0.0000001)) -  
      (log(max(abs(Thomas_call_stock_position_after_trade[,j-1] - Thomas_stock_position_after_trade[,j-1]) *Y[,j-1], 0.0000001)))
    #ditto  
    
    Difference_forward_replicating_portfolios[,j] <- Synthetic_forward_replicating_portfolio[,j] -
      Forward_replicating_portfolio[,j]
    
       

} # end of N steps


#=================#
   


Y_T <- Y[,N]

Y_T_log_change <- log(Y[,N]) - log(Y[,1])
  
 
Option_payoff_bought_put <-  pmax(K - Y_T,0)# pmax is "parallel maximum" of 2 vectors.
Option_payoff_bought_call <- pmax(Y_T - K,0)
  
#First for BS replication

BS_final_cash <- BS_cash_account[,j]
BS_final_stock <-  BS_stock_position_after_trade[,j]


BS_hedging_error <- ((BS_final_cash + BS_final_stock * Y_T) - Option_payoff_bought_put)
# ("proceeds of replication scheme" - "option payoff")

#Then for Thomas bought put replication
  

Thomas_final_cash <-  Thomas_cash_account[,j]
Thomas_final_stock <-  Thomas_stock_position_after_trade[,j]
Thomas_hedging_error <- (Thomas_final_cash + Thomas_final_stock * Y_T) - Option_payoff_bought_put
  
Thomas_final_cash_analytic <-  Thomas_cash_account_analytic[,j]
Thomas_final_stock_analytic <-  Thomas_stock_position_after_trade_analytic[,j]
Thomas_hedging_error_analytic <- (Thomas_final_cash_analytic + Thomas_final_stock_analytic * Y_T) - Option_payoff_bought_put

#Then for Thomas bought call replication

Thomas_call_final_cash <- Thomas_call_cash_account[,j]
Thomas_call_final_stock <-  Thomas_call_stock_position_after_trade[,j]


Thomas_call_hedging_error <- ((Thomas_call_final_cash + Thomas_call_final_stock * Y_T) - Option_payoff_bought_call)



payoff_expiry_call <-pmax(Y_T-K,0) 
expected_payoff_call <- mean(payoff_expiry_call)


Monte_Carlo_call_price <- exp(-r*(tau))*expected_payoff_call

payoff_expiry_put <-pmax(K-Y_T,0) 

#TRICK BELOW - using all the inner loop sims at time 2, n= 1 of 2 simulations, to improve the mean payoff
#expected_payoff_put <- mean(payoff_expiry_put) 

expected_payoff_put <- mean(Option_payoff_bought_put_temp[,2,,1]) 

Monte_Carlo_put_price <- exp(-r * tau) * expected_payoff_put  * exp(r * tau/N)


#First meansn for BS 
Mean_BS_final_cash <-mean(BS_final_cash)
Max_BS_final_cash <- max(BS_final_cash)
Min_BS_final_cash <- min(BS_final_cash)
Mean_BS_hedging_error <- mean(abs(BS_hedging_error))
St_dev_BS_hedging_error <- sd(BS_hedging_error)
St_dev_BS_hedging_error_graph  <- rbind(St_dev_BS_hedging_error,
                               c(log(N),log(St_dev_BS_hedging_error)))

#Then means for Thomas put
Thomas_mean_final_cash <-mean(Thomas_final_cash)
Thomas_max_final_cash <- max(Thomas_final_cash)
Thomas_min_final_cash <- min(Thomas_final_cash)
Thomas_mean_hedging_error <- mean(abs(Thomas_hedging_error))
Thomas_st_dev_hedging_error <- sd(Thomas_hedging_error)

Thomas_mean_hedging_error_analytic <- mean(abs(Thomas_hedging_error_analytic))
Thomas_st_dev_hedging_error_analytic <- sd(Thomas_hedging_error_analytic)

#Then means for Thomas call
Thomas_call_mean_final_cash <-mean(Thomas_call_final_cash)
Thomas_call_max_final_cash <- max(Thomas_call_final_cash)
Thomas_call_min_final_cash <- min(Thomas_call_final_cash)
Thomas_call_mean_hedging_error <- mean(abs(Thomas_call_hedging_error))
Thomas_call_st_dev_hedging_error <- sd(Thomas_call_hedging_error)


St_dev_hedging_errors_put_by_FD  <- rbind(St_dev_hedging_errors_put_by_FD,
                                    c(log(N),log(Thomas_st_dev_hedging_error)))

#Same for analytic - but only if fixed barrier, no analytic results relevant for trailing barrier
if(barrier_type == "fixed"){
St_dev_hedging_errors_put_analytic <- rbind(St_dev_hedging_errors_put_analytic,
                                            c(log(N),log(Thomas_st_dev_hedging_error_analytic)))
}

} # End of m'th trial

#Code below arranges histograms in group of 4. Edit the 2,2 2 to change, or edit out whole thing.

par(mfrow=c(2,2), oma = c(1,1,1,1) + 0.2,
    mar = c(3,3,1,1) + 1, mgp=c(1,1,0))

mybins <-seq(0.0, 0.1,0.001) #This can be used to make both histograms look same, with narrow bins.
mybins2 <-seq(-0.3, 0.3, 0.01) 
mybins3 <-seq(-0.05, 0.05, 0.001)

hist(Thomas_hedging_error, xlab="", ylim=c(0, nSim/3)) #, breaks = mybins2)
hist(Thomas_hedging_error_analytic, xlab="", ylim=c(0, nSim/3))  #, breaks = mybins2)

Rolled_up_difference_in_call_valuations <- 
  (Analytic_barrier_call - (S * exp(-q*tau) - K*exp(-r*tau) + Analytic_barrier_put)) * exp(r * tau)

Mean_diff_of_analytic_less_FD_delta <- sum(Thomas_delta_bought_put_analytic[,] - Thomas_delta_bought_put_by_FD[,])/(nSim * N)

cat("Analytic B-Scholes call      :",BS_call, " Analytic B-Scholes put    :",BS_put)

cat("\nAnalytic fixed barrier call  :",Analytic_barrier_call, " Analytic fixed barrier put:",Analytic_barrier_put)

cat("\nMC fixed barrier call        :",Monte_Carlo_call_price, 
    ifelse(barrier_type == "trailing", " MC trailing barrier put   :", " MC fixed barrier put      :"),Monte_Carlo_put_price)

cat("\nMn Thms abs hg err: Call (fxd b, anaytc):",round(Thomas_call_mean_hedging_error,digits=7), "   Put ditto     :"
    ,round(Thomas_mean_hedging_error_analytic,digits=7))

cat(ifelse(barrier_type == "trailing", "\nPut ditto (trailing b, Monte Carlo & FD):", 
           "\nPut ditto (fixed b, Monte Carlo & FD)   :") ,round(Thomas_mean_hedging_error,digits=7), "   Mean actual   :",
    round(mean(Thomas_hedging_error),digits=7))

if(barrier_type == "fixed"){cat("\nMean diff of analytic less FD delta (fixed barrier)  :" ,
                              round(Mean_diff_of_analytic_less_FD_delta,digits=7))}



library(beepr)
beep()

#=====================#

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#For the GNU General Public License, see <https://www.gnu.org/licenses/>.
