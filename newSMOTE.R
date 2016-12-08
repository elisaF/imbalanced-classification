mySMOTE <- function(data, 
					perc.over = 200, 
					k = 5, 
					perc.under = 200) 
	# INPUTS:
	# form a model formula
	# data the original training set with the unbalanced
	# perc.over / 100 is the number of new cases generated 
	#                 for each rare case. If perc.over < 100 a single case
	#				  is generated uniquely for a randomly selection perc.over
	#				  of the rare cases
	# k the number of nearest neighbors to consider as the pool from where the	  
	#   new examples are generated
	# perc.under / 100 is the numer of majority cases that are randomly selected 
	#                  selected for each smoted case
	# NOTE: this function is build specifically for the problem mentioned 
	#       in this paper
{
	# Get cases of minority class
	minExs <- which(data[,6] == 1)
	newExs <- mySMOTE.exs(data, perc.over, perc.over, k)
	
	# newExs <- new 
	
	# Undersample the majority class
	selMaj <- sample( ( 1:nrow(data) )[-minExs],
	as.integer((perc.under / 100) * nrow(newExs)), 
	replace = T)
	
	newdata <- na.omit(rbind(data[selMaj, ], data[minExs, ], newExs))
	
	
	return(newdata)
}

mySMOTE.exs <- function(data, N, k, tau2 = 0.001) 
	# INPUTS 
	# data the *entire* original dataset
	# N is the oversampling percentage
	# k the number of nearest neighbors to consider as the pool from where the	  
	#   new examples are generated
	# NOTE: this function is build specifically for the problem mentioned 
	#       in this paper
{
	# nomatr <- c()
	
	# Convert is_drug to integer
	
	# All data
	A <- data
	A[, 3] <- as.integer(A[, 3])
	
	nomatr <- 3
	
	# Choose only minority classes
	B <- A[A[, 6] == 1, ]
	
	if (N < 100) { # only take N% of data
		nB <- nrow(B)
		idx <- sample(1:nB, as.integer( (N/100) * nB) )
		B <- B[idx, ]
		N <- 100
	}

	p  <- ncol(B) 
	nB <- nrow(B)
	
	noise <- rep(0, p)
	
	ranges <- apply(A[, -6], 2, max) - apply(A[, -6], 2, min)
	nexs <- as.integer(N / 100)
	
	goal <- nrow(B) * nexs
	
	# Initialize matrix for new cases
	new <- NA
	
	kNNs.y <- matrix(nrow = nB, ncol = k)
	
	# initialize i
	# i <- 1
	# change for-loop to while-loop and increment i at end
	# if (i > nrow(B)) {i <- 1}
	
	check <- function(x) {
		if (typeof(x) != "list") {
			return(TRUE)
		} else {
			return(nrow(x) < goal)
		}
	}
	
	i <- 1
	while (check(new)) {
		if (i > nrow(B)) {
			i <- 1
		}
		
		index <- i
		i <- i + 1
		
		# print(i)
		
		# Find the k NNs of case B[, i]
		xd <- scale(A[, -6], as.numeric(B[i, -6]), ranges)
		xd[, 3] <- (xd[, 3] == 0)
		dd <- drop(xd^2 %*% rep(1, ncol(xd)))
		kNNs <- order(dd)[2:(k+1)]
		
		# print(i)
		
		# if (i == 2000) {break}
		
		# kNNs.y[i, ] <- A[kNNs, 6]
		kNNs.y <- A[kNNs, 6]
		m.x <- sum(kNNs.y == 1)
		m.y <- sum(kNNs.y == 0)
		
		if (m.x > m.y) { # Safe point
			cand <- kNNs           		 # Include all points
		} else if (m.x > 0) { # 
			d.x <- sum(xd[kNNs.y == 1])
			d.y <- sum(xd[kNNs.y == 0])
			if (d.x < d.y) {             # Boundary point
				cand <- kNNs[kNNs.y == 1]
			} else {		             # Noisy point
				cand <- c()
			}
		} else {
			cand <- c()					 # Noisy point
		}
		
		if (length(cand) == 0) { # Skip if there aren't any candidates
			next
		}
		
		noise <- rep(rnorm(5,0,0,5), p - 1)
		
		neig <- cand[sample(1:length(cand), 1)] # Choose one of the cands
		difs <- A[neig, -6] - B[i, -6] # Compute difference 
		neig.y <- A[neig, 6] # Class of chosen point
		
		if (neig.y == 1) {
			new.row <- B[i, -6] + runif(1, 0, 1) * difs + noise 
		} else {
			new.row <- B[i, -6] + runif(1, 0, 0.5) * difs + noise
		}
		
		# Handle is_drug
		new.row[, 3] <- c(A[neig, 3], B[i, 3])[1+round(runif(1),0)] 
		new.row <- cbind(new.row, B[i, 6])
		
		names(new.row) <- names(data)
		
		if (is.na(new)[1]) {
			new <- data.frame(new.row)
		} else{
			new <- rbind(new, new.row)
		}

		# print(nrow(new))
	}
	
	new[, 3] <- factor(new[, 3], levels=1:2, labels = levels(data[, 3]))
	return(new)
}

metrics <- function(pred, true) { 
	
	TP <- sum((true == 1) * (pred == 1))
	TN <- sum((true == 0) * (pred == 0))
	FP <- sum((true == 0) * (pred == 1))
	FN <- sum((true == 1) * (pred == 0))
	
	prec <- TP / (TP + FP)
	reca <- TP / (TP + FN)
	F1 <- 2 * (prec * reca) / (prec + reca)
	
	mylist <- list("precision" = prec, "recall" = reca, "F1" = F1)
	
	return(mylist)
}


ROC <- function(posterior, y) {
	quotient.post <- posterior[, 2] / (posterior[, 1] + 0.0001)
	
	cutoffs.NB <- exp(seq(-10, 4, length.out = 100))
	
	sens.vec <- rep(NA, length(cutoffs.NB)) # true positive
	spec.vec <- rep(NA, length(cutoffs.NB)) # true negative

	# For-loop for all cutoff values
	for (i in 1:length(cutoffs.NB)) {
		pred.i <- (quotient.post > cutoffs.NB[i])
	
		sens.vec[i] <- sum((pred.i == 1) * (y == 1)) / sum(y == 1)
		spec.vec[i] <- sum((pred.i == 0) * (y == 0)) / sum(y == 0)
	}
	
	sens.vec <- c(1, sens.vec)
	spec.vec <- c(0, spec.vec) 
 	 
	
	mylist <- list("sens" = sens.vec, "spec" = spec.vec)
	
	return(mylist)
}

AUC <- function(sens, spec) {
	require(caTools)
	auc <- trapz(1 - sens, spec)
	
	return(auc)
}
