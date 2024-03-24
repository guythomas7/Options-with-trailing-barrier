
RunNo <- 33


# After running simulation in other file, enter run number above to graph




par(mfrow=c(1,1), oma = c(1,1,1,1) + 0.2,  # CHANGE mfrow TO (2,2) for quad-plot 
   mar = c(3,3,1,1) + 1, mgp=c(1,1,0))
    
  
matplot(cbind(X[RunNo,1:N], Y_old[RunNo,1:N], Y[RunNo,1:N], b,0,1,b_trailing[RunNo,1:N], 
               pmax(1,B_running_max[RunNo,1:N])), # tempstrike),
        col=c( "red", "green4", "blue", "black", "black", "black",  "blue", "blue"), type="l", lty=c(1,1, 1, 2,2,2,2,2,1,1), lwd=0.1, xaxt = "n", yaxt = "n",
        main = "GBM reflected at trailing barrier (b x HWM)", xlab="Time steps", ylab="Price\n", cex.main = 1.6, cex.axis = 1.6, cex.lab = 1.6, xlim=c((0/52*N),(52/52*N)))   #ylim=c(0.9, 1.6))

text(4/52*N ,b-0.08,"fixed barrier", cex=1.4, col="black")
text(4/52*N ,b+0.08,"trailing barrier", cex=1.4, col="blue")
text(19/52*N ,1.5,"    Uprating factor", cex=1.4, col="blue")
text(35/52*N ,0.66,"    GBM", cex=1.4, col="red")
text(30/52*N ,2.5,"    GBM reflected at trailing barrier", cex=1.4, col="blue")

axis(1, at=c(1,N), cex.axis = 1.4, labels = c("0", expression(italic("T"))))
axis(2, at=c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), cex.axis = 1.4, labels = c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), las=2)


