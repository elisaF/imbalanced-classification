# Spencer Woody
#install.packages('DMwR')
#
library(caTools)
library(DMwR)
library(klaR)
library(ggplot2)

# Create train and test set

# Read in data, remove first column of names, create y as a factor
coref <- read.csv("coref_data.csv", header = T)
coref <- coref[, -1]
coref$y <- (coref$y == 1) + 0
coref$y <- as.factor(coref$y)
attach(coref)

# Summary of whole dataset, and for positive and negative instances
summary(coref)
summary(coref[y == 0, ])
summary(coref[y == 1, ])

head(coref)

N <- nrow(coref)
p <- ncol(coref) - 1

# Proportion of positives
sum(y == 1) / N
# [1] 0.09094212

# Proportion of negatives
sum(y == 0) / N
# [1] 0.9090579

# Create feature matrix for logistic regression
X <- as.matrix(cbind(rep(1, N), as.matrix(coref[, 1:p])))


#############################################################################
################################ Naive Bayes ################################
#############################################################################

nb <- NaiveBayes(y ~ ., data = coref, usekernel = TRUE) # Run the Naive Bayes
pred <- predict(nb, coref) # Predict values 

head(pred$class) # predicted class
head(pred$posterior) # posterior probabilities

sum(pred$class == y) / length(y) # ACCURACY
sum((pred$class == 1) * (y == 1)) / sum(y == 1) # TRUE POSITIVE (sens)
sum((pred$class == 0) * (y == 0)) / sum(y == 0) # TRUE NEGATIVE (spec)

# Create vector of cutoffs for ROC curve
cutoffs.NB <- exp(seq(-10, 4, length.out = 100))

# Create vector of ratio of posteriors (P(y = 1) / P(y = 0))
quotient.post <- pred$posterior[, 2] / (pred$posterior[, 1] + 0.0001)

# Create vectors for sensitivity and specificity
sens.NB <- rep(NA, length(cutoffs.NB)) # true positive
spec.NB <- rep(NA, length(cutoffs.NB)) # true negative

# For-loop for all cutoff values
for (i in 1:length(cutoffs.NB)) {
	pred.i <- (quotient.post > cutoffs.NB[i])
	
	sens.NB[i] <- sum((pred.i == 1) * (y == 1)) / sum(y == 1)
	spec.NB[i] <- sum((pred.i == 0) * (y == 0)) / sum(y == 0)
}

# Reverse the order (add point (1, 1))
spec.NB <- rev(c(0, spec.NB)) 
sens.NB <- rev(c(1, sens.NB))

# Create ROC curve (add point (1, 1) to plot)
plot(1 - spec.NB, sens.NB, 
	xlim = c(0, 1), 
	ylim = c(0, 1), 
	main = "ROC curve for Naive Bayes",
	type = "l",
	col = "red")
lines(c(0,1), c(0,1))
points(1 - spec.NB, sens.NB, pch = 19, col = "red")

# AUC of ROC
trapz(1 - spec.NB, sens.NB)

x <- c(0, 1)
y <- c(0, 1)

# Plot ROC curve
h <- qplot(1 - spec.NB, sens.NB, col = "firebrick3", geom = "line") + 
xlab("False Positive Rate") +
ylab("True Positive Rate") + 
ggtitle("ROC Curve for Naive Bayes") +
geom_line(aes(x = x, y = y), col = "black", size = 0.35) + 
 geom_point(aes(x = 1 - spec.NB, y = sens.NB), col = "firebrick3") +
theme(legend.position="none")

h

#############################################################################
################################## Logistic #################################
#############################################################################

# Run logistic regression
myglm <- glm(y ~ bow + tfidf + is_drug + coref_max_counts + coref_num_chains,
	data = coref, 
	family = "binomial")

# Create array of coefficients
mycoefs <- matrix(as.numeric(myglm$coef), nrow = p + 1)

# Fitted probabilities
fitted <- 1 / (1 + exp(  (- X %*% mycoefs) ))


sum((y == 1) * (fitted > 0.5)) / sum(y == 1) # True positive rate
# [1] 0.0873719

sum((fitted > 0.5) == (y == 1)) / N # Accuracy
# [1] 0.9100689


# Create ROC curve
cutoffs <- c(seq(log(0.01), log(0.1), length.out = 40), seq(log(0.11), log(0.99), length.out = 10))
cutoffs <- exp(cutoffs)

sens <- c()
spec <- c()

for (i in 1:length(cutoffs)) {
	cutoff.i <- cutoffs[i]
	
	pred.T <- (fitted > cutoff.i)
	pred.F <- (fitted <= cutoff.i)
	
	sens[i] <- sum((y == 1) * pred.T) / sum(y == 1) # true positive
	spec[i] <- sum((y == 0) * pred.F) / sum(y == 0) # true negative
}

plot(1 - spec, sens, 
	xlim = c(0, 1), 
	ylim = c(0, 1), 
	main = "ROC curve for logistic regression",
	type = "l",
	col = "red")
lines(c(0,1), c(0,1))
points(1 - spec, sens, pch = 19, col = "red")

trapz(rev(1 - spec), rev(sens))

########################################################################
################################ SMOTE #################################
########################################################################

SMOTEcoref <- SMOTE(y ~ ., data = coref, perc.over = 500, perc.under = 200)

summary(SMOTEcoref)

myglmSMOTE <- glm(y ~ bow + tfidf + is_drug + coref_max_counts + coref_num_chains, 
	data = SMOTEcoref,
	family = "binomial")

# Get coefficients
mycoefsSMOTE <- matrix(as.numeric(myglmSMOTE$coef), nrow = p + 1)

# Create estimated probabilities
fittedSMOTE <- 1 / (1 + exp( -X %*% mycoefsSMOTE ))

# True positive rate (specificity)
TP.SMOTE <- sum((fittedSMOTE > 0.5) * (y == 1)) / sum(y == 1)
print(TP.SMOTE)
# [1] 0.6459962

ACC.SMOTE <- sum((fitted2 > 0.5) == (y == 1)) / N
print(ACC.SMOTE)
# [1] 0.7823676

sensSMOTE <- c()
specSMOTE <- c()

for (i in 1:length(cutoffs)) {
	cutoff.i <- cutoffs[i]
	
	pred.T <- (fittedSMOTE > cutoff.i)
	pred.F <- (fittedSMOTE <= cutoff.i)
	
	sensSMOTE[i] <- sum((y == 1) * pred.T) / sum(y == 1) # true positive
	specSMOTE[i] <- sum((y == 0) * pred.F) / sum(y == 0) # true negative
}

plot(1 - specSMOTE, sensSMOTE, 
	xlim = c(0, 1), 
	ylim = c(0, 1), 
	main = "ROC curve for logistic regression with SMOTE",
	type = "l",
	col = "red")
lines(c(0,1), c(0,1))
points(1 - specSMOTE, sensSMOTE, pch = 19, col = "red")

trapz(rev(1 - specSMOTE), rev(sensSMOTE))



########################################################################
################################ small #################################
########################################################################

coref.pos <- coref[coref$y == 1, ]
coref.neg <- coref[coref$y == 0, ]
coref.neg <- coref.neg[1:nrow(coref.pos), ]

small <- rbind(coref.neg, coref.pos)

myglmsmall <- glm(y ~ bow + tfidf + is_drug + coref_max_counts + coref_num_chains, coref = small,
	 family = "binomial")
	
	
mymat2 <- as.matrix(cbind(rep(1, nrow(small)), small[, -p]))
mycoefs2 <- matrix(as.numeric(myglm2$coef), nrow = p + 1)

fitted <- 1 / (1 + exp(-mymat2 %*% mycoefs2 ))

# True positive rate
TP <- sum((small$y == 1) * (fitted > 0.5)) / sum(small$y == 1)
print(TP)
# > print(TP)
# [1] 0.5839847

# Accuracy
ACC <- sum((fitted > 0.5) == small$y) / (nrow(small))
print(ACC)
# [1] 0.9109851

cutoffs <- c(seq(log(0.01), log(0.1), length.out = 40), seq(log(0.11), log(0.99), length.out = 10))
cutoffs <- exp(cutoffs)

sens <- c()
spec <- c()

for (i in 1:length(cutoffs)) {
 	cutoff.i <- cutoffs[i]

 	pred.T <- (fitted > cutoff.i)
 	pred.F <- (fitted <= cutoff.i)

 	sens[i] <- sum((small$y == 1) * pred.T) / sum(small$y == 1)
 	spec[i] <- sum((small$y == 0) * pred.F) / sum(small$y == 0)
	 }

plot(1 - spec, sens, 
	xlim = c(0, 1), 
	ylim = c(0, 1), 
	main = "ROC curve for logistic regression",
	type = "l",
	col = "red")
lines(c(0,1), c(0,1))
points(1 - spec, sens, pch = 19, col = "red")

########################################################################
######################## PERF MEASURES #################################
########################################################################
get_perf <- function(true, predicted) {
  TP <- sum(true==1 & predicted==1)
  TN <- sum(true==0 & predicted==0)
  FP <- sum(predicted==1 & true==0)
  FN <- sum(predicted==0 & true==1)
  P <- TP / (TP + FP)
  R <- TP / (TP + FN)
  F1 <- 2 * (P*R) / (P+R)
  return(list(precision=P, recall=R, F1=F1))
}