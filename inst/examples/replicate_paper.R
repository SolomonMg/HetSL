# Replication sketch for paper results.
# For exact reproduction, use `replication/replicate_paper_results.R` and the
# pre-computed predictions bundled with the replication materials.

suppressPackageStartupMessages(library(HetEffects))

data("het_experiment_clean", package = "HetEffects")
data("blame_preds", package = "HetEffects")

# Inspect original ensemble weights from the replication materials.
print(blame_preds$Weights)

clean_names <- function(values) {
  values <- tolower(values)
  values <- gsub("[^a-z0-9]+", "_", values)
  values <- gsub("(^_+|_+$)", "", values)
  values
}

make_levels <- function(values, control = "control") {
  levels <- sort(unique(as.character(values)))
  if (control %in% levels) {
    c(control, setdiff(levels, control))
  } else {
    levels
  }
}

make_indicator_matrix <- function(values, levels, drop_first = TRUE) {
  values <- ifelse(is.na(values), levels[1], values)
  values <- factor(values, levels = levels)
  design <- stats::model.matrix(~ values - 1)
  colnames(design) <- levels
  if (drop_first && ncol(design) > 1) {
    design <- design[, -1, drop = FALSE]
  }
  as.data.frame(design)
}

make_interactions <- function(left, right) {
  left_matrix <- as.matrix(left)
  right_matrix <- as.matrix(right)
  interactions <- lapply(seq_len(ncol(left_matrix)), function(index) {
    left_matrix[, index] * right_matrix
  })
  interaction_matrix <- do.call(cbind, interactions)
  colnames(interaction_matrix) <- as.vector(
    outer(
      colnames(left_matrix),
      colnames(right_matrix),
      function(left_name, right_name) paste(left_name, right_name, sep = "_x_")
    )
  )
  as.data.frame(interaction_matrix)
}

build_covariates <- function(data) {
  dem <- as.integer(data$pid3l == "Dem")
  dem[is.na(dem)] <- 0L
  rep <- as.integer(data$pid3l == "Rep")
  rep[is.na(rep)] <- 0L
  lib <- as.integer(data$ideo3 %in% c(4, 5))
  lib[is.na(lib)] <- 0L
  cons <- as.integer(data$ideo3 < 3)
  cons[is.na(cons)] <- 0L

  data.frame(dem = dem, rep = rep, lib = lib, cons = cons)
}

build_treatments <- function(data) {
  type_levels <- make_levels(data$cond.type)
  num_levels <- c("control", "$20 million", "$50 thousand")
  stage_levels <- c("control", "requested", "secured", "will request")
  party_levels <- c("control", "a Republican", "a Democrat")
  along_levels <- c("control", "alone", "w/ Rep", "w/ Dem")

  type_matrix <- make_indicator_matrix(data$cond.type, type_levels)
  num_matrix <- make_indicator_matrix(data$cond.money, num_levels)
  stage_matrix <- make_indicator_matrix(data$cond.stage, stage_levels)
  party_matrix <- make_indicator_matrix(data$cond.party, party_levels)
  along_matrix <- make_indicator_matrix(data$cond.alongWith, along_levels)

  colnames(type_matrix) <- clean_names(colnames(type_matrix))
  colnames(num_matrix) <- clean_names(colnames(num_matrix))
  colnames(stage_matrix) <- clean_names(colnames(stage_matrix))
  colnames(party_matrix) <- clean_names(colnames(party_matrix))
  colnames(along_matrix) <- clean_names(colnames(along_matrix))

  num_main <- num_matrix[, "20_million", drop = FALSE]
  stage_main <- stage_matrix[, c("requested", "secured"), drop = FALSE]
  party_main <- party_matrix[, "a_republican", drop = FALSE]
  along_main <- along_matrix[, c("alone", "w_rep"), drop = FALSE]

  # Mirror RepCode.R by omitting the final type from interactions.
  interaction_types <- type_matrix[, seq_len(min(5L, ncol(type_matrix))), drop = FALSE]

  interaction_terms <- cbind(
    make_interactions(interaction_types, num_main),
    make_interactions(interaction_types, stage_main[, "requested", drop = FALSE]),
    make_interactions(interaction_types, stage_main[, "secured", drop = FALSE]),
    make_interactions(interaction_types, party_main),
    make_interactions(interaction_types, along_main[, "alone", drop = FALSE]),
    make_interactions(interaction_types, along_main[, "w_rep", drop = FALSE]),
    make_interactions(num_main, stage_main),
    make_interactions(num_main, party_main),
    make_interactions(num_main, along_main),
    make_interactions(stage_main, party_main),
    make_interactions(stage_main, along_main),
    make_interactions(party_main, along_main)
  )

  treatments <- cbind(
    type_matrix,
    num_main,
    stage_main,
    party_main,
    along_main,
    interaction_terms
  )

  treatments
}

svdat <- het_experiment_clean
covariates <- build_covariates(svdat)
treatments <- build_treatments(svdat)

dat <- data.frame(approve_bi = svdat$approve_bi, covariates, treatments)

# Replication controls.
fast_mode <- FALSE
nfolds <- if (fast_mode) 1L else 10L
requested_learners <- c(
  "lasso",
  "elastic_net_0.5",
  "elastic_net_0.25",
  "findit",
  "bayesglm",
  "bart",
  "randomforest",
  "krls",
  "svm_smo"
)

installed_learners <- available_learners()$name
learners <- intersect(requested_learners, installed_learners)
missing_learners <- setdiff(requested_learners, learners)
if (length(missing_learners) > 0) {
  message("Skipping unavailable learners: ", paste(missing_learners, collapse = ", "))
}

set.seed(123)
fit <- het_ensemble(
  data = dat,
  outcome = "approve_bi",
  treatment = colnames(treatments),
  covariates = colnames(covariates),
  learners = learners,
  nfolds = nfolds,
  family = "binomial",
  parallel = TRUE,
  verbose = TRUE
)

print(fit$weights)

# MCATE is not defined for multi-treatment models in this API.
