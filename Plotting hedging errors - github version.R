#Sets up charts in a quad - and they look nicer when copied from this

par(mfrow=c(2,2), oma = c(1,1,1,1) + 0.2,
    mar = c(3,3,1,1) + 1, mgp=c(1,1,0))


plot(St_dev_hedging_errors_put_by_FD, xlab = "Log (Number of time steps)", ylab="Log(SD of replication errors)", ylim=c(-7.1,-2.8), 
     xlim=c(2,8), cex.axis = 1.2, cex.lab = 1.33, mgp = c(2.5,1,0), col="red" )

#mgp sets distance of axis title, axis numbers ("labels"), axis ticks. Default is c(3,1.0). IN previous paper, you probably uses yaxt=n, yaxt = n and
# then applies axes by separate axis() statements aterwards. But the wauy here wll do. 


Best_fit_by_FD <- lm(St_dev_hedging_errors_put_by_FD[,2]~ St_dev_hedging_errors_put_by_FD[,1])
abline(Best_fit_by_FD, col="red")

if(barrier_type == "fixed"){
  points(St_dev_hedging_errors_put_analytic[,1], St_dev_hedging_errors_put_analytic[,2])
  
  Best_fit_analytic <- lm(St_dev_hedging_errors_put_analytic[,2]~ St_dev_hedging_errors_put_analytic[,1])
  
  abline(Best_fit_analytic)

  legend("bottomleft", legend=c("using analytic delta", "using finite difference delta"), pch = 1, lty = 1, col=c("red", "black")) 

  #title(main = "Fixed barrier", cex.main = 1.0)
  
  }else{ # only the data for FD delta, if trailing barrier - no relevaant analytical reusults

  legend("bottomleft", legend=c("using finite difference delta"), pch = 1, lty = 1, col=c("red")) 
    
  #title(main = "Trailing barrier", cex.main = 1.0)
  }

text( 4.9, -3.0,  ifelse(barrier_type == "fixed", "Fixed barrier", "Trailing barrier"), cex = 1.4)


