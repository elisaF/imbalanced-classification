# Spencer Woody

library(caTools)
library()

# Read in data, remove first column of names, create y as a factor
coref <- read.csv("coref_data.csv", header = T)
coref <- coref[, -1]
coref$y <- (coref$y == 1) + 0
coref$y <- as.factor(coref$y)
attach(coref)

# Create feature matrix
X <- as.matrix(cbind(rep(1, N), as.matrix(coref[, 1:p])))

N <- nrow(coref)
p <- ncol(coref) - 1

# Summary of whole dataset, and for positive and negative instances
summary(coref)
summary(coref[y == 0, ])
summary(coref[y == 1, ])

head(coref)

# Proportion of positives
sum(y == 1) / N
# [1] 0.09094212

# Proportion of negatives
sum(y == 0) / N
# [1] 0.9090579

# Run logistic regression
myglm <- glm(y ~ bow + tfidf + is_drug + coref_max_counts + coref_num_chains,
	data = coref, 
	family = "binomial")

# Create array of coefficients
mycoefs <- matrix(as.numeric(myglm$coef), nrow = p + 1)

fitted <- 1 / (1 + exp(  (- X %*% mycoefs) ))

# True positive rate
TP <- sum((y == 1) * (fitted > 0.5)) / sum(y == 1)
print(TP)
# [1] 0.0873719

# Accuracy
ACC <- sum((fitted > 0.5) == (y == 1)) / N
print(ACC)
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

SMOTEcoref <- SMOTE(y ~ ., data = coref, perc.over = 300, perc.under = 100)

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

# coref.pos <- coref[coref$y == 1, ]
# coref.neg <- coref[coref$y == 0, ]
# coref.neg <- coref.neg[1:nrow(coref.pos), ]
#
# small <- rbind(coref.neg, coref.pos)
#
# myglmsmall <- glm(y ~ bow + tfidf + is_drug + coref_max_counts + coref_num_chains, coref = small,
# 	 family = "binomial")
#
#
# mymat2 <- as.matrix(cbind(rep(1, nrow(small)), small[, -p]))
# mycoefs2 <- matrix(as.numeric(myglm2$coef), nrow = p + 1)
#
# fitted <- 1 / (1 + exp(-mymat2 %*% mycoefs2 ))
#
# # True positive rate
# TP <- sum((small$y == 1) * (fitted > 0.5)) / sum(small$y == 1)
# print(TP)
# # > print(TP)
# # [1] 0.5839847
#
# # Accuracy
# ACC <- sum((fitted > 0.5) == small$y) / (nrow(small))
# print(ACC)
# # [1] 0.9109851
#
# cutoffs <- c(seq(log(0.01), log(0.1), length.out = 40), seq(log(0.11), log(0.99), length.out = 10))
# cutoffs <- exp(cutoffs)
#
# sens <- c()
# spec <- c()
#
# for (i in 1:length(cutoffs)) {
#  	cutoff.i <- cutoffs[i]
#
#  	pred.T <- (fitted > cutoff.i)
#  	pred.F <- (fitted <= cutoff.i)
#
#  	sens[i] <- sum((small$y == 1) * pred.T) / sum(small$y == 1)
#  	spec[i] <- sum((small$y == 0) * pred.F) / sum(small$y == 0)
# 	 }
#
# plot(1 - spec, sens,
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	main = "ROC curve for logistic regression",
# 	type = "l",
# 	col = "red")
# lines(c(0,1), c(0,1))
# points(1 - spec, sens, pch = 19, col = "red")
#
#