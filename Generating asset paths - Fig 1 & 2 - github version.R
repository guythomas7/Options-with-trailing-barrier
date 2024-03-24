#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purposes.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at  R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

# Setting parameters 
S <- 1.0 #stock price at time t
b <- 0.5 #barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
tau <- 10 #time to maturity T - t (in years) 
r <- 0.015 #risk-free annual interest rate (convenient to set to zero)
q <- 0.01 #deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
sigma <- 0.30 #annual volatility of the stock price (standard deviation)

set.seed(1930) #set the seed (if not set, it's taken from the computer clock)
N <- 6300 #N is number of time steps, e.g. 252 x 25 = 6300 steps for every working day over 25 years. 
nSim <- 100 #number of simulations (paths) 


#Check validity of inputs
stopifnot(b <= min(S,K), r!=q, K!=b)


# Compute the Monte Carlo prices 

dt <- tau/N #length of each time sub interval
Z <-  matrix(rnorm(nSim*N, mean=0, sd=1),nrow = nSim, ncol = N) #standard normal sample of N elements
dW <- Z*sqrt(dt) #Brownian motion increments (nSim simulations) x (N increments)
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_trailing <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_T <- vector(length=nSim) 

b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_old <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)


X[,1] <- S
Y[,1] <- S
b_trailing[,1] <- b
B_running_max[,1] <- b_trailing[,1]/X[,1]
Y_running_max[,1] <- S

Y_old[,1] <- X[,1]
b_as_matrix[,] <- b
B_running_max_old[,1] <- b_as_matrix[,1]/X[,1]


for(j in 2:N){
    
    X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt + sigma*dW[,j])
    B[,j] <- b_trailing[,j-1]/X[,j] #Have to use b_trailing[j-1] because don't know new b_trailing before calculating latest Y
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    Y_running_max[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], Y[,j], Y_running_max[,j-1])
    b_trailing[,j] <- ifelse(Y[,j] > Y_running_max[,j-1], b * Y[,j], b_trailing[,j-1])
    
    B_old[,j] <- b_as_matrix[,j]/X[,j]
    B_running_max_old[,j] <- ifelse(B_old[,j] > B_running_max_old[,j-1], B_old[,j], B_running_max_old[,j-1])
    Y_old[,j] <- X[,j]*pmax(1,B_running_max_old[,j])
    
 
}  # end of N steps
 

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
