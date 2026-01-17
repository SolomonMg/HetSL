# Replication sketch for paper results.
# This script is a best-effort approximation using the package API.

library(HetEffects)

data("het_experiment_clean")
data("blame_preds")

# Inspect original ensemble weights from the replication materials.
print(blame_preds$Weights)

# Use a similar (but not identical) learner set supported by the package.
# Note: the original replication includes learners not yet implemented here
# (e.g., FindIt). Expect differences in weights and predictions.

# Follow RepCode.R treatment construction (credit claiming experiment).
svdat <- het_experiment_clean

type_levels <- sort(unique(as.character(svdat$cond.type)))
type_mat <- matrix(0, nrow = nrow(svdat), ncol = length(type_levels))
colnames(type_mat) <- type_levels
for (i in seq_len(nrow(svdat))) {
  type_mat[i, which(colnames(type_mat) == svdat$cond.type[i])] <- 1
}
type_mat_final <- type_mat[, -1, drop = FALSE]

num_levels <- c("control", "$20 million", "$50 thousand")
num_mat <- matrix(0, nrow = nrow(svdat), ncol = length(num_levels))
colnames(num_mat) <- num_levels
for (i in seq_len(nrow(svdat))) {
  num_mat[i, which(colnames(num_mat) == svdat$cond.money[i])] <- 1
}
num_mat_final <- num_mat[, -1, drop = FALSE]

stage_levels <- c("control", "requested", "secured", "will request")
stage_mat <- matrix(0, nrow = nrow(svdat), ncol = length(stage_levels))
colnames(stage_mat) <- stage_levels
for (i in seq_len(nrow(svdat))) {
  stage_mat[i, which(colnames(stage_mat) == svdat$cond.stage[i])] <- 1
}
stage_mat_final <- stage_mat[, -1, drop = FALSE]

party_levels <- c("control", "a Republican", "a Democrat")
party_mat <- matrix(0, nrow = nrow(svdat), ncol = length(party_levels))
colnames(party_mat) <- party_levels
for (i in seq_len(nrow(svdat))) {
  party_mat[i, which(colnames(party_mat) == svdat$cond.party[i])] <- 1
}
party_mat_final <- party_mat[, -1, drop = FALSE]

along_levels <- c("control", "alone", "w/ Rep", "w/ Dem")
along_mat <- matrix(0, nrow = nrow(svdat), ncol = length(along_levels))
colnames(along_mat) <- along_levels
for (i in seq_len(nrow(svdat))) {
  along_mat[i, which(colnames(along_mat) == svdat$cond.alongWith[i])] <- 1
}
along_mat_final <- along_mat[, -1, drop = FALSE]

dem <- ifelse(svdat$pid3l == "Dem", 1, 0)
dem[is.na(dem)] <- 0
rep <- ifelse(svdat$pid3l == "Rep", 1, 0)
rep[is.na(rep)] <- 0
lib <- ifelse(svdat$ideo3 %in% c(4, 5), 1, 0)
lib[is.na(lib)] <- 0
cons <- ifelse(svdat$ideo3 < 3, 1, 0)
cons[is.na(cons)] <- 0

covs <- cbind(Dem = dem, Rep = rep, Lib = lib, Cons = cons)

type_mat_final <- data.frame(type_mat_final)
colnames(type_mat_final) <- c("PlanParent", "Parks", "Gun_Range", "Fire", "Police", "Roads")
num_mat_final <- data.frame(num_mat_final)
colnames(num_mat_final) <- c("mil_20", "thou_50")
stage_mat_final <- data.frame(stage_mat_final)
colnames(stage_mat_final) <- c("request", "secure", "will")
party_mat_final <- data.frame(party_mat_final)
colnames(party_mat_final) <- c("rep_rep", "dem_rep")
along_mat_final <- data.frame(along_mat_final)
colnames(along_mat_final) <- c("alone", "w_rep", "w_dem")

treats <- cbind(
  type_mat_final,
  num_mat_final[, 1],
  stage_mat_final[, 1:2],
  party_mat_final[, 1],
  along_mat_final[, 1:2],
  type_mat_final[, 1:5] * num_mat_final[, 1],
  type_mat_final[, 1:5] * stage_mat_final[, 1],
  type_mat_final[, 1:5] * stage_mat_final[, 2],
  type_mat_final[, 1:5] * party_mat_final[, 1],
  type_mat_final[, 1:5] * along_mat_final[, 1],
  type_mat_final[, 1:5] * along_mat_final[, 2],
  num_mat_final[, 1] * stage_mat_final[, 1],
  num_mat_final[, 1] * stage_mat_final[, 2],
  num_mat_final[, 1] * party_mat_final[, 1],
  num_mat_final[, 1] * along_mat_final[, 1],
  num_mat_final[, 1] * along_mat_final[, 2],
  stage_mat_final[, 1:2] * party_mat_final[, 1],
  stage_mat_final[, 1:2] * along_mat_final[, 1],
  stage_mat_final[, 1:2] * along_mat_final[, 2],
  party_mat_final[, 1] * along_mat_final[, 1],
  party_mat_final[, 1] * along_mat_final[, 2]
)

type_short <- c("PlanParent", "Parks", "Gun_Range", "Fire", "Police", "Roads")
num_short <- c("mil_20", "thou_50")
stage_short <- c("request", "secure", "will")
party_short <- c("rep_rep", "dem_rep")
along_short <- c("alone", "w_rep", "w_dem")

colnames(treats)[1:12] <- c(type_short, num_short[1], stage_short[1:2], party_short[1], along_short[1:2])
colnames(treats)[13:17] <- paste(type_short[1:5], num_short[1], sep = "_x_")
colnames(treats)[18:22] <- paste(type_short[1:5], stage_short[1], sep = "_x_")
colnames(treats)[23:27] <- paste(type_short[1:5], stage_short[2], sep = "_x_")
colnames(treats)[28:32] <- paste(type_short[1:5], party_short[1], sep = "_x_")
colnames(treats)[33:37] <- paste(type_short[1:5], along_short[1], sep = "_x_")
colnames(treats)[38:42] <- paste(type_short[1:5], along_short[2], sep = "_x_")
colnames(treats)[43] <- paste(num_short[1], stage_short[1], sep = "_x_")
colnames(treats)[44] <- paste(num_short[1], stage_short[2], sep = "_x_")
colnames(treats)[45] <- paste(num_short[1], party_short[1], sep = "_x_")
colnames(treats)[46] <- paste(num_short[1], along_short[1], sep = "_x_")
colnames(treats)[47] <- paste(num_short[1], along_short[2], sep = "_x_")
colnames(treats)[48:49] <- paste(stage_short[1:2], party_short[1], sep = "_x_")
colnames(treats)[50:51] <- paste(stage_short[1:2], along_short[1], sep = "_x_")
colnames(treats)[52:53] <- paste(stage_short[1:2], along_short[2], sep = "_x_")
colnames(treats)[54] <- paste(party_short[1], along_short[1], sep = "_x_")
colnames(treats)[55] <- paste(party_short[2], along_short[2], sep = "_x_")

dat <- data.frame(approve_bi = svdat$approve_bi, covs, treats)

# Replication controls
fast_mode <- FALSE
nfolds <- if (fast_mode) 1 else 10

learners <- c(
  "lasso", "elastic_net_0.75", "elastic_net_0.5", "elastic_net_0.25",
  "findit", "bayesglm", "glmboost", "bart", "randomforest",
  "glm", "krls", "svm_smo"
)

# Optional profiling to identify and skip slow learners.
profile_sample <- 250
max_fit_seconds <- 120

profile_learners <- function(learners, dat, covs, treats) {
  set.seed(123)
  idx <- seq_len(min(profile_sample, nrow(dat)))
  X <- as.matrix(dat[idx, covs, drop = FALSE])
  Tmat <- as.matrix(dat[idx, treats, drop = FALSE])
  y <- dat$approve_bi[idx]

  timings <- data.frame(
    learner = learners,
    fit_seconds = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(learners)) {
    lname <- learners[i]
    spec <- HetEffects:::.get_learner(lname)
    t <- system.time({
      model <- tryCatch(
        spec$fit_fn(X, y, Tmat, family = "binomial"),
        error = function(e) NULL
      )
    })
    timings$fit_seconds[i] <- t[["elapsed"]]
  }

  timings
}

timings <- profile_learners(learners, dat, colnames(covs), colnames(treats))
print(timings)

if (!is.null(max_fit_seconds)) {
  keep <- timings$fit_seconds <= max_fit_seconds
  learners <- timings$learner[keep]
  message("Keeping learners: ", paste(learners, collapse = ", "))
}

fit <- het_ensemble(
  data = dat,
  outcome = "approve_bi",
  treatment = colnames(treats),
  covariates = colnames(covs),
  learners = learners,
  nfolds = nfolds,
  family = "binomial",
  parallel = TRUE,
  verbose = TRUE
)

print(fit$weights)

# MCATE is not defined for multi-treatment models in this API.
print(fit$weights)
