library(brms)
# requires brms 1.10.0 or higher
theme_set(theme_default())

### ----- Example 1: Catching fish -----

# Load the data
zinb <- read.csv("http://stats.idre.ucla.edu/stat/data/fish.csv")
zinb$camper <- factor(zinb$camper, labels = c("no", "yes"))
head(zinb)

# Fit a basic zero-inflated model
fit_zinb1 <- brm(count ~ persons + child + camper, 
                 data = zinb, family = zero_inflated_poisson())

# Obtain summaries
summary(fit_zinb1)
plot(marginal_effects(fit_zinb1), ask = FALSE)

# Fit a more advanced zero-inflated model
fit_zinb2 <- brm(bf(count ~ persons + child + camper, zi ~ child), 
                 data = zinb, family = zero_inflated_poisson())

# Obtain summaries
summary(fit_zinb2)
plot(marginal_effects(fit_zinb2), ask = FALSE)

# Compare model fit
LOO(fit_zinb1, fit_zinb2)


### ----- Example 2: Housing Rents -----

# Load the data
data("rent99", package = "gamlss.data")
head(rent99)

# Only model the mean of the response distribution
fit_rent1 <- brm(rentsqm ~ t2(area, yearc) + (1|district),
                 data = rent99, chains = 2, cores = 2)

# summarize and vizualize results
summary(fit_rent1)
marginal_effects(fit_rent1, surface = TRUE)


# Model both mean and SD of the response distribution
bform <- bf(rentsqm ~ t2(area, yearc) + (1|ID1|district),
            sigma ~ t2(area, yearc) + (1|ID1|district))
fit_rent2 <- brm(bform, data = rent99, chains = 2, cores = 2)

# summarize and vizualize results
summary(fit_rent2)
marginal_smooths(fit_rent2)

# Compare model fit
LOO(fit_rent1, fit_rent2)


### ----- Example 3: Insurance loss payments -----

# Load the data
url <- paste0("https://raw.githubusercontent.com/mages/",
              "diesunddas/master/Data/ClarkTriangle.csv")
loss <- read.csv(url)
head(loss)

# Fit the model
nlform <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
             ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1, 
             nl = TRUE)

nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
             prior(normal(1, 2), nlpar = "omega"),
             prior(normal(45, 10), nlpar = "theta")) 

fit_loss1 <- brm(formula = nlform, data = loss, 
                 family = gaussian(), prior = nlprior,
                 control = list(adapt_delta = 0.9))

# Obtain summaries
summary(fit_loss1)
marginal_effects(fit_loss1)

# Visualize cumulative insurance loss separately for each year
conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)
me_year1 <- marginal_effects(fit_loss1, conditions = conditions, 
                             re_formula = NULL, method = "predict")
plot(me_year1, ncol = 5, points = TRUE)

# Let omega and theta als vary by year
nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
              ult + omega + theta ~ 1 + (1|ID1|AY),
              nl = TRUE)
fit_loss2 <- update(fit_loss1, formula = nlform2,
                    control = list(adapt_delta = 0.90))

# Compare model fit
LOO(fit_loss1, fit_loss2)


### ----- Example 4: Performance of school children -----

# Define the simulation function
sim_multi_mem <- function(nschools, nstudents = nschools * 100, 
                          change = 0.1, ymean = 20,
                          var_overall = 20, icc = 0.3,
                          w1 = 0.5, w2 = 0.5, seed = 1234) {
  # simulate data for a simple multi-membership model
  # Args:
  #   nschools: total number of schools
  #   nstudents: total number of students
  #   change: percentage of students changing school during the year
  #   ymean: mean performance of students across schools
  #   var_overall: overall variance between students
  #   icc: intra-class-correlation; percentage of overall variance 
  #        explained by schools
  #   w1, w2: used to weight schools
  #   seed: used by set.seed to make results reproducible
  # Returns:
  #   a data.frame with columns s1, s2, w1, w2, and y (the performance)
  stopifnot(icc >= 0, icc <= 1)
  stopifnot(change >= 0, change <= 1)
  set.seed(seed)
  var_schools <- var_overall * icc
  var_resid <- var_overall - var_schools
  eff_schools <- rnorm(nschools, ymean, sqrt(var_schools))
  nstudents_change <- round(nstudents * change)
  students_change <- sample(1:nstudents, nstudents_change)
  students_stay <- setdiff(1:nstudents, students_change)
  schools_change <- t(replicate(nstudents_change, sample(1:nschools, 2)))
  schools_stay <- sample(1:nschools, length(students_stay), replace = TRUE)
  schools_stay <- cbind(schools_stay, schools_stay)
  data <- as.data.frame(rbind(schools_change, schools_stay))
  names(data) <- c("s1", "s2")
  data$w1 <- w1
  data$w2 <- w2
  data$y <- data$w1 * eff_schools[data$s1] + 
    data$w2 * eff_schools[data$s2] + 
    rnorm(nstudents, 0, sqrt(var_resid))
  data
}

# Simulate some data
data_mm <- sim_multi_mem(nschools = 10, nstudents = 1000, change = 0.1)
head(data_mm)
data_mm[101:106, ]

# Fit a multi-membership model
fit_mm1 <- brm(y ~ 1 + (1|mm(s1, s2)), data = data_mm)

# Summarize the model and visualize model fit
summary(fit_mm1)
pp_check(fit_mm1)

# Specify non-equal weights for the two schools
data_mm[1:100, "w1"] <- runif(100, 0, 1)
data_mm[1:100, "w2"] <- 1 - data_mm[1:100, "w1"]
head(data_mm)

# Fit a second multi-membership model explicitly incorporating weights
fit_mm2 <- brm(y ~ 1 + (1|mm(s1, s2, weights = cbind(w1, w2))), 
               data = data_mm)
summary(fit_mm2)
