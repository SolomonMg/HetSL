# End-to-end example using the paper dataset.

library(HetEffects)

data("het_experiment")

covars <- c("pid3", "ideo3", "gender", "race", "educ", "inc", "byear")
dat <- het_experiment[, c("approval", "party_treat", covars)]

fit <- het_ensemble(
  approval ~ pid3 + ideo3 + gender + race + educ + inc + byear | party_treat,
  data = dat,
  learners = "lasso",
  nfolds = 5,
  family = "gaussian",
  verbose = TRUE
)

print(fit)
summary(fit)

effects <- treatment_effects(fit, moderators = "pid3")
print(effects$ate)
print(effects$mcate$pid3)
