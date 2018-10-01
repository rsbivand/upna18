##########################################
#### R code to run the examples in
#### Spatio-Temporal Areal Unit Modeling in R with Conditional
#### Autoregressive Priors Using the CARBayesST Package
##########################################

set.seed(1234)

################################################
#### Section 4.1 Generating and simulating data
################################################
#### Create the regular spatial grid and time domain
n.space <- 20
N <- 20
x.easting <- 1:n.space
x.northing <- 1:n.space
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
N.all <- N * K

#### Create the binary spatial adjacency (rook) matrix W
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1

#### Create the binary temporal adjacency matrix D
distance <- as.matrix(dist(1:N))
D <- array(0, c(N, N))
D[distance == 1] <- 1
    
#### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
Q.W <- 0.8 * (diag(apply(W, 2, sum)) - W) + 0.2 * diag(rep(1, K))

#### Create the spatial covariance matrix and simulate a set of spatial random effects
Q.W.inv <- solve(Q.W)
library("mvtnorm")
phi <- rmvnorm(n = 1, mean = rep(0, K), sigma = (0.01 * Q.W.inv))
phi.long <- rep(phi, N)
    
#### Create the temporal covariance matrix and simulate a set of temporal random effects
Q.D <- 0.8 * (diag(apply(D, 2, sum)) - D) + 0.2 * diag(rep(1, N))
Q.D.inv <- solve(Q.D)
delta <- rmvnorm(n = 1, mean = rep(0, N), sigma = (0.01 * Q.D.inv))
delta.long <- kronecker(delta, rep(1, K))

#### Generate a single covariate and a set of interactions
x <- rnorm(n = N.all, mean = 0, sd = 1)
gamma <- rnorm(n = N.all, mean = 0, sd = 0.1)

#### Fix the remaining parameters in the simulation and generate a data set
beta1 <- 0
beta2 <- 0.1
n <- rep(50, N.all)
LP <- beta1 + beta2 * x + phi.long +  delta.long + gamma
theta.true <- exp(LP) / (1 + exp(LP))
Y <- rbinom(n = N.all, size = n, prob = theta.true)

#### Fit the ST.CARanova() model to the simulated data
library("CARBayesST")
model <- ST.CARanova(formula = Y ~ x, family = "binomial", trials = n,
                     W = W, burnin = 20000, n.sample = 120000, thin = 10)

#### Plot the posterior distributions of the regression parameters
colnames(model$samples$beta) <- c("beta1", "beta2")
plot(model$samples$beta)

#### View a model summary
model

###############################################
#### Section 4.2 Small simulation study
###############################################
#### This section repeats the simuation set up from above and 
#### Generate 100 data sets and computes bias and coverage probabilities of the
#### 95% uncertainty intervals

#### Specify fixed parameters
tau2 <- 0.01
rho <- 0.8

#### Matrices to save results
n.sim <- 100
results.beta <- array(NA, c(n.sim, 4))
results.rho <- array(NA, c(n.sim, 4))
results.tau2 <- array(NA, c(n.sim, 6))
results.re <- array(NA, c(n.sim, 6))
results.fitted <- array(NA, c(n.sim, 2))

#######################
#### Run the simulation
#######################
for (i in 1:n.sim) {
    ##################
    #### Generate data
    ##################
    #### Random effects
    phi <- rmvnorm(n = 1, mean = rep(0,K), sigma = (tau2 * Q.W.inv))
    phi <- phi - mean(phi)
    phi.long <- rep(phi, N)
    delta <- rmvnorm(n = 1, mean = rep(0,N), sigma = (tau2 * Q.D.inv))
    delta <- delta - mean(delta)
    delta.long <- kronecker(delta, rep(1,K))
    gamma <- rnorm(n = N.all, mean = 0, sd = sqrt(tau2))
    gamma <- gamma - mean(gamma)
    
    #### Linear predictor and data
    x <- rnorm(n = N.all, mean = 0, sd = 1)
    LP <- beta1 + beta2 * x + phi.long +  delta.long + gamma
    theta.true <- exp(LP) / (1 + exp(LP))
    Y <- rbinom(n = N.all, size = n, prob = theta.true)
    
    #### Fit the model
    model <- ST.CARanova(formula = Y ~ x, family = "binomial",
                         trials = n, W = W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
    
    #### Store the results
    results.beta[i, 1:2] <- as.numeric(apply(model$samples$beta, 2, median))
    results.beta[i, 3] <- as.numeric(quantile(model$samples$beta[, 1], 0.025) < beta1 & quantile(model$samples$beta[, 1], 0.975) > beta1)
    results.beta[i, 4] <- as.numeric(quantile(model$samples$beta[, 2], 0.025) < beta2 & quantile(model$samples$beta[, 2], 0.975) > beta2)
    
    results.rho[i, 1:2] <- as.numeric(apply(model$samples$rho, 2, median))     
    results.rho[i, 3] <- as.numeric(quantile(model$samples$rho[, 1], 0.025) < rho & quantile(model$samples$rho[, 1], 0.975) > rho)
    results.rho[i, 4] <- as.numeric(quantile(model$samples$rho[, 2], 0.025) < rho & quantile(model$samples$rho[, 2], 0.975) > rho)
    
    results.tau2[i, 1:3] <- as.numeric(apply(model$samples$tau2, 2, median))    
    results.tau2[i, 4] <- as.numeric(quantile(model$samples$tau2[, 1], 0.025) < tau2 & quantile(model$samples$tau2[, 1], 0.975) > tau2)
    results.tau2[i, 5] <- as.numeric(quantile(model$samples$tau2[, 2], 0.025) < tau2 & quantile(model$samples$tau2[, 2], 0.975) > tau2)
    results.tau2[i, 6] <- as.numeric(quantile(model$samples$tau2[, 3], 0.025) < tau2 & quantile(model$samples$tau2[, 3], 0.975) > tau2)
    
    results.re[i, 1] <- mean(apply(model$samples$phi, 2, median) - phi)
    quants <- t(apply(model$samples$phi, 2, quantile, c(0.025, 0.975)))
    results.re[i, 4] <- sum(as.numeric(quants[, 1] < phi & quants[, 2] > phi))
    
    results.re[i, 2] <- mean(apply(model$samples$delta, 2, median) - delta)
    quants <- t(apply(model$samples$delta, 2, quantile, c(0.025, 0.975)))
    results.re[i, 5] <- sum(as.numeric(quants[, 1] < delta & quants[, 2] > delta))
    
    results.re[i, 3] <- mean(apply(model$samples$gamma, 2, median) - gamma)
    quants <- t(apply(model$samples$gamma, 2, quantile, c(0.025, 0.975)))
    results.re[i, 6] <- sum(as.numeric(quants[, 1] < gamma & quants[, 2] > gamma))  
    
    results.fitted[i, 1] <- mean(model$fitted.values - theta.true * n)
    quants <- t(apply(model$samples$fitted, 2, quantile, c(0.025, 0.975)))
    results.fitted[i, 2] <- sum(as.numeric(quants[, 1] < theta.true * n & quants[, 2] > theta.true * n))  
    
    #### Remove the model object
    rm(model)
}

##########################
#### Summarise the results
##########################
#### beta
mean(results.beta[, 1]) - beta1
mean(results.beta[, 2]) - beta2
mean(results.beta[, 3])
mean(results.beta[, 4])

#### rho
mean(results.rho[, 1]) - rho
mean(results.rho[, 2]) - rho
mean(results.rho[, 3])
mean(results.rho[, 4])

#### tau2
mean(results.tau2[, 1]) - tau2
mean(results.tau2[, 2]) - tau2
mean(results.tau2[, 3]) - tau2
mean(results.tau2[, 4])
mean(results.tau2[, 5])
mean(results.tau2[, 6])

#### RE
mean(results.re[, 1])
mean(results.re[, 2])
mean(results.re[, 3])
mean(results.re[, 4]) / K
mean(results.re[, 5]) / N
mean(results.re[, 6]) / N.all

#### Fitted
mean(results.fitted[, 1])
mean(results.fitted[, 2]) / N.all

#####################################################
#### Section 4.3 timing based on different data sizes
#####################################################
#### First create a function to simulate data of various sizes
dat.sim <- function(n.space, N, rho, tau2, beta1, beta2, trials) {
    #### Create the regular spatial grid and time domain
    x.easting <- 1:n.space
    x.northing <- 1:n.space
    Grid <- expand.grid(x.easting, x.northing)
    K <- nrow(Grid)
    N.all <- N * K
  
    #### Create the binary spatial adjacency (rook) matrix W
    distance <- as.matrix(dist(Grid))
    W <- array(0, c(K, K))
    W[distance == 1] <- 1
    
    #### Create the binary temporal adjacency matrix D
    distance <- as.matrix(dist(1:N))
    D <- array(0, c(N, N))
    D[distance == 1] <- 1
    
    #### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
    Q.W <- rho * (diag(apply(W, 2, sum)) - W) + (1 - rho) * diag(rep(1, K))
    Q.W.inv <- solve(Q.W)
    phi <- rmvnorm(n = 1, mean = rep(0, K), sigma = (tau2 * Q.W.inv))
    phi.long <- rep(phi, N)
    
    #### Create the temporal covariance matrix and simulate a set of temporal random effects
    Q.D <- rho * (diag(apply(D, 2, sum)) - D) + (1 - rho) * diag(rep(1, N))
    Q.D.inv <- solve(Q.D)
    delta <- rmvnorm(n = 1, mean = rep(0, N), sigma = (tau2 * Q.D.inv))
    delta.long <- kronecker(delta, rep(1, K))
    
    #### Generate a single covariate and a set of interactions
    x <- rnorm(n = N.all, mean = 0, sd = 1)
    gamma <- rnorm(n = N.all, mean = 0, sd = 0.1)
    
    #### Fix the remaining parameters in the simulation and generate a data set
    n <- rep(trials, N.all)
    LP <- beta1 + beta2 * x + phi.long +  delta.long + gamma
    theta.true <- exp(LP) / (1 + exp(LP))
    Y <- rbinom(n = N.all, size = n, prob = theta.true)
    
    #### Return the results
    results <- list(Y = Y, n = n, x = x, W = W)
    return(results)
}

#### 10 by 10 grid and 10 time periods
dat <- dat.sim(n.space = 10, N = 10, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 10 by 10 grid and 20 time periods
dat <- dat.sim(n.space = 10, N = 20, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 20 by 20 grid and 10 time periods
dat <- dat.sim(n.space = 20, N = 10, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 20 by 20 grid and 20 time periods
dat <- dat.sim(n.space = 20, N = 20, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 30 by 30 grid and 20 time periods
dat <- dat.sim(n.space = 30, N = 20, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 30 by 30 grid and 30 time periods
dat <- dat.sim(n.space = 30, N = 30, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 40 by 40 grid and 30 time periods
dat <- dat.sim(n.space = 40, N = 30, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 40 by 40 grid and 40 time periods
dat <- dat.sim(n.space = 40, N = 40, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

#### 50 by 50 grid and 40 time periods
dat <- dat.sim(n.space = 50, N = 40, rho = 0.8, tau2 = 0.01, beta1 = 0, beta2 = 1, trials = 50)
system.time(
    model <- ST.CARanova(formula = Y ~ x, data = dat, family = "binomial",
                         trials = dat$n, W = dat$W, burnin = 20000,
                         n.sample = 120000, thin = 10, verbose = FALSE)
)

##############################################################
#### Section 5 of the paper - air pollution and health example
##############################################################

set.seed(1234)

#### Load the data
library("CARBayesdata")
library("sp")
data("GGHB.IG", package = "CARBayesdata")
data("pollutionhealthdata", package = "CARBayesdata")

#### View the structure of the data
head(pollutionhealthdata)

#### Pairs plot of the data
pollutionhealthdata$SMR <- with(pollutionhealthdata, observed / expected)
pollutionhealthdata$logSMR <- log(pollutionhealthdata$SMR)
par(pty = "s", cex.axis = 1.5)
pairs(pollutionhealthdata[, c(9, 5:7)], pch = 19, cex = 0.5, lower.panel = NULL, panel = panel.smooth,
      labels = c("ln(SMR)", "PM10", "JSA", "Price (*100,000)"))

#### Aggregate the data to produce a map and add to the GGHB.IG object
library("dplyr")
SMR.av <- summarise(group_by(pollutionhealthdata, IG), SMR.mean = mean(SMR))
GGHB.IG@data$SMR <- SMR.av$SMR.mean

#### Plot the map of the average SMR
l1 <- list("SpatialPolygonsRescale", layout.north.arrow(),
           offset = c(220000, 647000), scale = 4000)
l2 <- list("SpatialPolygonsRescale", layout.scale.bar(),
           offset = c(225000, 647000), scale = 10000, fill = c("transparent", "black"))
l3 <- list("sp.text", c(225000,649000), "0")
l4 <- list("sp.text", c(230000,649000), "5000 m")
breakpoints <- seq(min(SMR.av$SMR.mean) - 0.1, max(SMR.av$SMR.mean) + 0.1, length.out = 11)
spplot(GGHB.IG, "SMR", sp.layout = list(l1, l2, l3, l4),
       xlab = "Easting", ylab = "Northing", 
       scales = list(draw = TRUE), at = breakpoints, 
       col.regions = terrain.colors(n = length(breakpoints) - 1),
       par.settings = list(fontsize = list(text = 20)))

#### Create the spatial information
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = SMR.av$IG)
W.list <- nb2listw(W.nb, style = "B")
W <- nb2mat(W.nb, style = "B")

#### Run a simple Poisson model
formula <- observed ~ offset(log(expected)) + jsa + price + pm10
model1 <- glm(formula = formula, family = "quasipoisson", data = pollutionhealthdata)
resid.glm <- residuals(model1)
coef(summary(model1))
summary(model1)$dispersion

#### Conduct a Moran's I test on the residuals
moran.mc(x = resid.glm[1:271], listw = W.list, nsim = 10000)

#### Run the CARBayesST model
library("CARBayesST")
model2 <- ST.CARar(formula = formula, family = "poisson", data = pollutionhealthdata,
                   W = W, burnin = 20000, n.sample = 220000, thin = 10)

#### Print the summary of the model
model2

#### Show how the samples are stored
summary(model2$samples)

#### Plot the relative risk distributions
par(pty = "m")
colnames(model2$samples$beta) <- c("Intercept", "JSA", "Price", "PM10")
plot(exp(model2$samples$beta[, -1]))

#### Summarise the relative risk distributions
library("CARBayes")
parameter.summary <- summarise.samples(exp(model2$samples$beta[, -1]),
                                       quantiles = c(0.5, 0.025, 0.975))
round(parameter.summary$quantiles, 3)

####################################################
#### Section 6 of the paper - housing market example
####################################################

set.seed(1234)

#### Load the data
library("CARBayesdata")
library("sp")
data("GGHB.IG", package = "CARBayesdata")
data("salesdata", package = "CARBayesdata")

#### Draw boxplots of the data
salesdata$salesprop <- with(salesdata, sales / stock)
boxplot(salesprop ~ year, data = salesdata, range = 0, xlab = "Year",
        ylab = "Property sales rate", col = "darkseagreen", border = "navy")

#### Draw a map of the average spatial pattern
library("dplyr")
salesprop.av <- summarise(group_by(salesdata,IG), salesprop.mean = mean(salesprop))
GGHB.IG@data$sales <- salesprop.av$salesprop.mean
l1 <- list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000, 647000), 
           scale = 4000)
l2 <- list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000, 647000), 
          scale = 10000, fill = c("transparent","black"))
l3 <- list("sp.text", c(225000, 649000), "0")
l4 <- list("sp.text", c(230000, 649000), "5000 m")
breakpoints <- c(0, quantile(salesprop.av$salesprop.mean, seq(0.1, 0.9, 0.1)), 0.1)
spplot(GGHB.IG, "sales", sp.layout = list(l1, l2, l3, l4),
       xlab = "Easting", ylab = "Northing", 
       scales = list(draw  =  TRUE), at = breakpoints,
       col.regions = terrain.colors(n = length(breakpoints) - 1),
       par.settings = list(fontsize = list(text = 20)))

#### Create the neighbourhood matrix
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = salesprop.av$salesprop.mean)
W <- nb2mat(W.nb, style = "B")

#### Fit the model
library("CARBayesST")
formula <- sales ~ offset(log(stock))
model1 <- ST.CARsepspatial(formula = formula, family = "poisson", data = salesdata,
                           W = W, burnin = 20000, n.sample = 220000, thin = 10)

#### Calculate the temporal trends
trend.mean <- array(NA, c(11, 3))
trend.sd <- array(NA, c(11, 3))
for (i in 1:11) {
    posterior <- exp(model1$samples$phi[, ((i-1)*271 + 1):(i*271)] +
                     matrix(rep(model1$samples$beta + model1$samples$delta[, i], 271),
                            ncol = 271, byrow = FALSE))
    trend.mean[i, ] <- quantile(apply(posterior, 1, mean), c(0.5, 0.025, 0.975))
    trend.sd[i, ] <- quantile(apply(posterior, 1, sd), c(0.5, 0.025, 0.975))
}

#### Plot the temporal trends
par(mfrow = c(2, 1))
plot(jitter(salesdata$year), salesdata$salesprop, pch = 19, cex = 0.2, col = "blue",
     xlab = "Year", main = "(a)", ylab = "Average sales rate", ylim = c(0, 0.11),
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
lines(2003:2013, trend.mean[, 1], col = "red", type = "l")
lines(2003:2013, trend.mean[, 2])
lines(2003:2013, trend.mean[, 3])

plot(2003:2013, trend.sd[, 1], col = "red", type = "l", xlab = "Year", main = "(b)",
     ylab = "Spatial standard deviation", ylim = c(0, 0.06), cex.axis = 1.5,
     cex.lab = 1.5, cex.main = 1.5)
lines(2003:2013, trend.sd[, 2])
lines(2003:2013, trend.sd[, 3])

#### Map the fitted proportions for the odd numbered years
rate.est <- matrix(model1$fitted.values / salesdata$stock, nrow = nrow(W), byrow = FALSE)
rate.est <- as.data.frame(rate.est)
colnames(rate.est) <- c("r2003", "r2004", "r2005", "r2006", "r2007", "r2008", "r2009",
                        "r2010", "r2011", "r2012", "r2013")
GGHB.IG@data <- data.frame(GGHB.IG@data, rate.est)
spplot(GGHB.IG, c("r2011", "r2013", "r2007", "r2009", "r2003", "r2005"),
       names.attr = c("Rate 2011", "Rate 2013", "Rate 2007", "Rate 2009", "Rate 2003", "Rate 2005"),
       sp.layout = list(l1, l2, l3, l4), xlab = "Easting", ylab = "Northing",
       scales = list(draw  =  TRUE), at = breakpoints,
       col.regions = terrain.colors(n = length(breakpoints - 1)),
       par.settings = list(fontsize = list(text = 13)))

