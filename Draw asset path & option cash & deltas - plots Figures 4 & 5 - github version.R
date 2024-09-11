
RunNo <- 34

# Enter Run number above that you want to graph


par(mfrow=c(1,1), oma = c(1,1,1,1) + 0.2,  # CHANGE mfrow TO (2,2) for quad-plot 
   mar = c(3,3,1,1) + 1, mgp=c(1,1,0))
    
    
matplot(cbind(X[RunNo,1:N], Y[RunNo,1:N], Thomas_cash_account[RunNo,1:N]    ,   Thomas_delta_bought_put_by_FD[RunNo, 1:N] , 
                           b,0,1,                       
                          if(barrier_type == "trailing"){b_trailing[RunNo,1:N]}, 
                          #if(barrier_type == "trailing"){pmax(1,B_running_max[RunNo,1:N])}, 
                          if(barrier_type == "fixed"){Thomas_delta_bought_put_analytic[RunNo, 1:N]}), 
        col=c( "red", "blue", "green4" , "orange" ,"black", "black", "black", "blue", "chocolate4","blue"), 
        type="l", lty=c(1,1,1,1,3,3,3, 2,2,1), lwd=0.1, xaxt = "n",
        main = ifelse(barrier_type=="fixed", "GBM reflected at fixed barrier b", "GBM reflected at trailing barrier (b x HWM)"),
        xlab="Time steps", ylab="Price\n", cex.main = 1.2, cex.axis = 1.2, 
        cex.lab = 1.2, xlim=c((0/52*N),(52/52*N))) #, ylim=c(0.55, 0.65))   

text(4/52*N ,b-0.04,"fixed barrier", cex=0.6, col="black")

axis(1, at=c(1,N), cex.axis = 1.2, labels = c("0", expression(italic("T"))))
text(4/52*N , Thomas_cash_account[2/52*N]-0.15,"   cash", cex=1.0, col="green4")
text(3/52*N , Thomas_delta_bought_put_by_FD[4/52*N]-0.15,"  FD delta", cex=1.0, col="orange")

if(barrier_type == "fixed"){
  text(4/52*N , Thomas_delta_bought_put_by_FD[4/52*N]-0.23,"   analytic delta", cex=1.0, col="purple1")
}

if(barrier_type == "trailing"){
text(4/52*N ,b_trailing[RunNo,100] + 0.05,"   trailing barrier", cex=0.6, col="blue")
#text(4/52*N ,1.05,"    Uprating factor", cex=0.6, col="chocolate4")
  
}
  
legend("topleft", ifelse(barrier_type == "fixed","FIXED","TRAILING"))

