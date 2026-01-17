test_that("het_ensemble fits and predicts with formula interface", {
  set.seed(123)
  n <- 80
  dat <- data.frame(
    y = rnorm(n),
    treat = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  fit <- het_ensemble(
    y ~ x1 + x2 | treat,
    data = dat,
    learners = "lasso",
    nfolds = 3,
    family = "gaussian",
    verbose = FALSE
  )

  expect_s3_class(fit, "het_ensemble")
  expect_equal(sum(fit$weights), 1, tolerance = 1e-6)

  pred <- predict(fit)
  expect_length(pred, n)

  cf <- predict(fit, type = "counterfactual")
  expect_equal(nrow(cf), n)
  expect_true(all(c("y0", "y1", "cate") %in% names(cf)))

  ate_res <- ate(fit)
  expect_true(is.finite(ate_res$estimate))
})

test_that("het_ensemble fits with explicit arguments", {
  set.seed(456)
  n <- 60
  dat <- data.frame(
    y = rnorm(n),
    treat = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  fit <- het_ensemble(
    data = dat,
    outcome = "y",
    treatment = "treat",
    covariates = c("x1", "x2"),
    learners = "lasso",
    nfolds = 3,
    family = "gaussian",
    verbose = FALSE
  )

  expect_s3_class(fit, "het_ensemble")
  cate_vals <- cate(fit)
  expect_length(cate_vals, n)
})

test_that("estimate_propensity returns valid scores", {
  set.seed(789)
  n <- 100
  dat <- data.frame(
    treat = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  ps <- estimate_propensity(treat ~ x1 + x2, data = dat, method = "logit")
  expect_true(is.list(ps))
  expect_length(ps$propensity, n)
  expect_true(all(ps$propensity > 0 & ps$propensity < 1))
})
