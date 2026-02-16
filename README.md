# HetEffects

Ensemble methods for estimating heterogeneous treatment effects,
implementing the approach described in [Grimmer, Messing & Westwood
(2017)](https://doi.org/10.1017/pan.2017.2).

The package combines multiple machine learning methods using
cross-validation to optimize ensemble weights, then estimates treatment
effects that vary across subgroups defined by observed covariates. It
supports both randomized experiments and observational studies (via
propensity score weighting), and provides analytic standard errors via
influence functions as well as bootstrap inference.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("solomonmessing/HetEffects")
```

## Quick start

```r
library(HetEffects)

# Load included experiment data
data("het_experiment")

# Fit an ensemble using the formula interface
# outcome ~ covariates | treatment
fit <- het_ensemble(
  approval ~ pid3 + ideo3 + gender + race + educ + inc + byear | party_treat,
  data = het_experiment,
  learners = c("lasso", "elastic_net_0.5"),
  nfolds = 5,
  family = "gaussian"
)

# Average treatment effect with analytic SEs
ate(fit)

# Treatment effects moderated by partisanship
effects <- treatment_effects(fit, moderators = "pid3")
effects$mcate$pid3

# Conditional average treatment effects (per-observation)
cate_vals <- cate(fit)
```

## Available learners

The ensemble can draw from any combination of the following learners.
Learners whose dependencies are not installed are automatically skipped.

| Learner | Package | Category |
|---------|---------|----------|
| LASSO | glmnet | Core |
| Elastic net (0.75, 0.5, 0.25) | glmnet | Core |
| Bayesian GLM | arm | Classic |
| BART | dbarts / BayesTree | Classic |
| Random forest | randomForest | Classic |
| SVM | e1071 | Classic |
| KRLS | KRLS | Classic |
| FindIt | FindIt | Classic |
| GBM | gbm | Classic |
| GLMBoost | mboost | Classic |
| Causal forest | grf | Modern |
| XGBoost | xgboost | Modern |

Check which learners are available in your current installation:

```r
available_learners()
```

Register custom learners with `register_learner()`.

## Learner tuning

Pass multiple parameter configurations per learner to create distinct
ensemble candidates:

```r
fit <- het_ensemble(
  y ~ x1 + x2 | treat,
  data = mydata,
  learners = c("lasso", "svm"),
  learner_params = list(
    svm = list(
      list(cost = 0.1, kernel = "radial"),
      list(cost = 1,   kernel = "radial"),
      list(cost = 10,  kernel = "polynomial")
    )
  )
)
```

Each configuration is treated as its own learner variant and receives a
separate ensemble weight.

## Inference

```r
# Analytic standard errors (influence-function based)
analytic_se(fit, effect_type = "ate")
analytic_se(fit, effect_type = "mcate", moderator = "pid3")

# Bootstrap inference
boot <- bootstrap_effects(fit, R = 500)
confint(fit)
```

## Observational studies

For non-randomized data, supply propensity scores to adjust for
confounding:

```r
fit <- het_ensemble(
  y ~ x1 + x2 | treat,
  data = obs_data,
  propensity = ~ x1 + x2 + x3
)
```

Or provide a pre-estimated propensity score column:

```r
fit <- het_ensemble(
  y ~ x1 + x2 | treat,
  data = obs_data,
  propensity = "pscore"
)
```

## Interactive exploration

Launch a Shiny app to explore fitted results:

```r
view_effects(fit)
```

## Included datasets

- **`het_experiment`** / **`het_experiment_clean`** -- Survey experiment
  on government spending approval (Grimmer, Messing & Westwood 2017).
- **`blame_preds`** -- Predicted effects from a blame attribution
  experiment.

## Citation

If you use this package, please cite:

> Grimmer, J., Messing, S., & Westwood, S. J. (2017). Estimating
> heterogeneous treatment effects and the effects of heterogeneous
> treatments with ensemble methods. *Political Analysis*, 25(4),
> 413--434. <https://doi.org/10.1017/pan.2017.2>

## License

MIT
