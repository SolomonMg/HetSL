test_that("end-to-end workflow runs on het_experiment dataset", {
  skip_if_not(exists("het_experiment", where = asNamespace("HetEffects"), inherits = FALSE))

  data("het_experiment", package = "HetEffects", envir = environment())

  covars <- c("pid3", "ideo3", "gender", "race", "educ", "inc", "byear")
  dat <- het_experiment[, c("approval", "party_treat", covars)]
  dat <- dat[complete.cases(dat), ]
  dat <- dat[seq_len(min(200, nrow(dat))), ]

  fit <- het_ensemble(
    approval ~ pid3 + ideo3 + gender + race + educ + inc + byear | party_treat,
    data = dat,
    learners = "lasso",
    nfolds = 3,
    family = "gaussian",
    verbose = FALSE
  )

  expect_s3_class(fit, "het_ensemble")
  expect_true(all(is.finite(predict(fit))))
  expect_true(is.finite(ate(fit)$estimate))
})
