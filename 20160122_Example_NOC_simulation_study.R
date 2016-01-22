########################################
######  Tamar Sofer
#####  2016-01-22
########################################
###########  Example simulation study for NOC paper
##########  the true effect size is 3.


#####  setting variables for simulations
mean.u0 <- 0
mean.u1 <- 2
sd.u <- 1.5


mean.w0 <- 1
mean.w1 <- 3
sd.w <- sd.u

sd.y <- 3
sd.n <- 1.5

alpha <- c(1,2)
beta <- c(2,3)



n.sim <- 1000
dis <- "norm" ## for (dis in dist){
n <- 200        ## 	for (n in n.samp/2){
A <- c(rep(0,n), rep(1,n))
est.eff <- rep(0, n.sim)
		
for (i in 1:n.sim){
	
		set.seed(i)
		C <- cbind(1, rnorm(n*2))
		colnames(C) <- c("intercept", "C1")

		if (dis == "norm"){
	  		U0 <- rnorm(n,mean.u0 ,sd.u)
		   U1 <- rnorm(n,mean.u1,sd.u)
		   W0 <- rnorm(n,mean.w0,sd.w)
		   W1 <- rnorm(n,mean.w1,sd.w)
		  } else {
		   		U0 <- runif(n,1,9 )
		   		U1 <- runif(n,3,13)
		   		W0 <- runif(n,2,10)
		   		W1 <- runif(n,4,14)
		   }



		Y0 <- (U0 + C[1:n,] %*% alpha)*sd.y
		Y1 <- (U1 + C[(n+1):(2*n),] %*% alpha + 1)*sd.y  ### true effect is 1*sd.y = 3. 
	
		Y <- c(Y0,Y1)

		N <- (c(W0,W1) + C %*% beta )*sd.n
		N0 <- N[1:n]
		N1 <- N[(n+1):(2*n)]

		temp <- LocScaleNOC(Y, N, A, C, verbose = F)

		est.eff[i] <- temp$alpha
	
		if (i %% 100 == 0) print(i)
}

## the bias is estimated by: 
mean(est.eff) - 3