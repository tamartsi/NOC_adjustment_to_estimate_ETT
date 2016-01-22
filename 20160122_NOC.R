########################################
######  Tamar Sofer
#####  2016-01-22
########################################
###### code for estimating the ETT of A on Y, while using the negative control outcome N to correct for unobserved confounding.
#### Y is the outcome
#### N negative control
#### A exposure
##### C a matrix of covariates
##### model.variance = T if assuming heteroesedasticity and a linear model of the magnitude of variance should be fit
##### assume.qq = T if one assumes that the residuals of Y and N in the unexposed have the same distribution.
##### verbose = T if want to print messages

LocScaleNOC <- function(Y, N, A, C, model.variance = F, assume.qq = F, verbose = F){
	
	if (!is.numeric(C)) {
		stop("The design matrix C must be numeric")
	}
	C <- as.matrix(C)
	### assumes first column in C is intercept
	if (!all(C[,1] == 1)) C <- cbind(1, C)
	
	if (!all(is.element(A, c(0,1)))) stop("A must have only 0 or 1 values")
	
	stopifnot(length(Y) == length(N), length(Y) == length(A), nrow(C) == length(Y))
	
	na.inds <- c(which(is.na(Y)), which(is.na(N)), which(is.na(A)), which( apply(C, 1, function(x) sum(is.na(x))) > 0 ) )
	if (length(na.inds) > 0){
		message(paste("Removing ", length(na.inds), " observations with missing values"))
		Y <- Y[-na.inds]
		N <- N[-na.inds]
		A <- A[-na.inds]
		C <- C[-na.inds,]
	}

	n <- length(Y)	
	
	if (verbose) message(paste("Estimating the ETT using ", n, "observations") )
	
	inds.A0 <- which(A == 0)
	inds.A1 <- which(A == 1)
	
	N0 <- N[inds.A0]
	N1 <- N[inds.A1]
	
	Y0 <- Y[inds.A0]
	Y1 <- Y[inds.A1]
	
	N0.model.mat <- C[inds.A0,]
	N1.model.mat <- C[inds.A1,]

	Y0.model.mat <- C[inds.A0,]
	Y1.model.mat <- C[inds.A1,]

	
	mod.N0 <- lm(N0 ~ -1 +  N0.model.mat)
	mod.Y0 <- lm(Y0 ~ -1 + Y0.model.mat)

	beta.N0 <- coef(mod.N0)
	beta.Y0 <- coef(mod.Y0)
	
	if (!model.variance) {
		s.n <- sd(resid(mod.N0))
		s.y <- sd(resid(mod.Y0))
		
		N1.scaled.by.N0 <- (N1 - N1.model.mat %*% beta.N0)/s.n
		
		if (assume.qq) Y1.scaled <- N1.scaled.by.N0*s.y + Y1.model.mat %*% beta.Y0 
		
	} else{
		
		mod.resid.sq.N0 <- glm(resid(mod.N0)^2 ~ -1 + N0.model.mat, family = quasi(variance = "mu", link = "log"))
		mod.resid.sq.Y0 <- glm(resid(mod.Y0)^2 ~ -1 + Y0.model.mat, family = quasi(variance = "mu", link = "log"))
		
		beta.sn <- coef(mod.resid.sq.N0)
		beta.sy <- coef(mod.resid.sq.Y0)
		
		N1.scaled.by.N0 <- (N1 - N1.model.mat %*% beta.N0)/sqrt(exp(N1.model.mat %*% beta.sn))
		
		if (assume.qq) Y1.scaled <- N1.scaled.by.N0*sqrt(exp(Y1.model.mat %*% beta.sy)) + Y1.model.mat %*% beta.Y0 
		
	}
	

	if (!assume.qq){
		
		### we need to prepare maps between values and cumulative probabilities, in the no-treatment groups. 
		### we use the residuals of Y0 and N0 after modelling. 
		
		## first: find the scaled-down residuals
		if (!model.variance){
			scaled.resids.N0 <- 	resid(mod.N0)/s.n
			scaled.resids.Y0 <- resid(mod.Y0)/s.y
		} else{
			
			scaled.resids.N0 <- resid(mod.N0)/sqrt(exp(N0.model.mat %*% beta.sn))
			scaled.resids.Y0 <- resid(mod.Y0)/sqrt(exp(Y0.model.mat %*% beta.sy))
		}
		
		## second: prepare maps between values and probabilities
		map.Y0.resids <- matrix(0, nrow = 2, ncol = length(Y0))
		map.Y0.resids[1,] <- scaled.resids.Y0[order(scaled.resids.Y0)]
		map.Y0.resids[2,] <-  sapply(map.Y0.resids[1,], function(u) {mean(map.Y0.resids[1,] <= u)})
		
		map.N0.resids <- matrix(0, nrow = 2, ncol = length(N0))
		map.N0.resids[1,] <- scaled.resids.N0[order(scaled.resids.N0)]
		map.N0.resids[2,] <-  sapply(map.N0.resids[1,], function(u) {mean(map.N0.resids[1,] <= u)})
		
		## prepare function that gives probabilities of values according to a map
		match.prob.to.val <- function(val, map){
			if (is.element(val, map[1,])) return(map[2, which(map[1,] == val)]) 
			if (map[1,1] > val) return(0)
			if (map[1,ncol(map)] < val) return(1)
			ind <- which.max(val - map[1,] < 0)
			return(map[2,ind -1] + (map[2, ind] - map[2, ind -1])*(val - map[1, ind -1])/(map[1, ind ] - map[1, ind - 1])  )
		}
		
		## prepare function that gives values for probabilities according to a map
		match.val.to.prob <- function(prob, map){		
			if (is.element(prob, map[2,])) return(map[1, which(map[2,] == prob)]) 
			if (prob < map[2,1]) return(map[1,1])
			if (prob > map[2,ncol(map)]) return(map[1,ncol(map)])
			ind <- which.min(prob - map[2,] > 0)
			return(map[1,ind -1] + (map[1, ind ] - map[1, ind -1])*(prob - map[2, ind - 1])/(map[2, ind] - map[2, ind - 1])  )
		}
 
 
 		### now match the scaled-by-N0 N1 values with scaled Y values
		prob.N1 <- sapply(N1.scaled.by.N0, function(u) {match.prob.to.val(u, map.N0.resids)})
		matched.Y1 <- sapply(prob.N1, function(u) {match.val.to.prob(u, map.Y0.resids)})
				
		if (!model.variance) Y1.scaled <- matched.Y1*s.y + Y1.model.mat %*% beta.Y0 else{ ## we modellled the variance
			
			Y1.scaled <- matched.Y1*sqrt(exp(Y1.model.mat %*% beta.sy)) + Y1.model.mat %*% beta.Y0
		}
		
		
	}
	


	if (model.variance) return(list(
	beta.y0 = beta.Y0,
	beta.n0 = beta.N0,
	beta.sy = beta.sy,
	beta.sn = beta.sn,
	alpha = mean(Y1 - Y1.scaled) 
	)) else return(list(
		beta.y0 = beta.Y0,
		beta.n0 = beta.N0,
		beta.sy = s.y,
		beta.sn = s.n,
		alpha = mean(Y1 - Y1.scaled) 
		))
	
}



