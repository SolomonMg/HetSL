# Prepare package datasets from the paper replication files.

wd <- getwd()
root <- normalizePath(file.path(wd, ".."))
if (basename(wd) == "data-raw") {
  root <- normalizePath(file.path(wd, ".."))
} else if (basename(wd) == "HetEffects") {
  root <- normalizePath(file.path(wd, ".."))
} else {
  root <- normalizePath(wd)
}

source_dir <- root

blame_path <- file.path(source_dir, "BlamePreds.RData")
het_path <- file.path(source_dir, "Het_Experiment.RData")

if (!file.exists(blame_path)) {
  blame_path <- file.path(source_dir, "HET_rep", "BlamePreds.RData")
}
if (!file.exists(het_path)) {
  het_path <- file.path(source_dir, "HET_rep", "Het_Experiment.RData")
}

load(blame_path)       # loads preds
load(het_path)         # loads svdat

# Rename for package datasets
blame_preds <- preds
het_experiment <- svdat

# Add a binary treatment indicator for party cue vs control.
het_experiment$party_treat <- ifelse(het_experiment$cond.party == "control", 0, 1)

# Create a cleaned dataset following RepCode.R preprocessing.
svdat <- het_experiment
names(svdat)[23:51] <- c(
  "preq1", "preq2", "preq3",
  "nextc", "contr", "nextt", "treat",
  "approval", "therm", "fiscRespbl",
  "bringMoneyEff", "passLegEff",
  "secReqMC", "likGetM", "daysGetM", "break",
  "gender", "race", "byear", "ntvEnglish",
  "ideo3", "voted", "pid3", "pidCloser", "educ",
  "inc", "finalinst", "howLong", "comments"
)

approv <- agrep("I pay attention", max.distance = 0.3, svdat$comments)
approv2 <- agrep("I PAY ATTENTION", max.distance = 0.3, svdat$comments)
approv <- c(approv, approv2)
svdat <- svdat[approv, ]

svdat$cond.type[svdat$contr == 1] <- "control"
svdat$cond.type <- stats::relevel(factor(svdat$cond.type), ref = "control")
svdat$cond.money[svdat$contr == 1] <- "control"
svdat$cond.money <- stats::relevel(factor(svdat$cond.money), ref = "control")
svdat$cond.stage[svdat$contr == 1] <- "control"
svdat$cond.stage <- stats::relevel(factor(svdat$cond.stage), ref = "control")
svdat$cond.party[svdat$contr == 1] <- "control"
svdat$cond.party <- stats::relevel(factor(svdat$cond.party), ref = "control")
svdat$cond.alongWith[svdat$contr == 1] <- "control"
svdat$cond.alongWith <- stats::relevel(factor(svdat$cond.alongWith), ref = "control")
levels(svdat$cond.alongWith) <- c("control", "alone", "w/ Dem", "w/ Rep")

svdat$pid3l <- factor(c("Dem", "Rep", "Ind/Oth", "Ind/Oth")[svdat$pid3])
svdat$pid3l <- stats::relevel(svdat$pid3l, ref = "Ind/Oth")

svdat$approve_bi <- ifelse(svdat$approval < 3, 1, 0)
svdat$party_treat <- ifelse(svdat$cond.party == "control", 0, 1)

het_experiment_clean <- svdat

data_dir <- file.path(root, "HetEffects", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

save(blame_preds, file = file.path(data_dir, "blame_preds.rda"))
save(het_experiment, file = file.path(data_dir, "het_experiment.rda"))
save(het_experiment_clean, file = file.path(data_dir, "het_experiment_clean.rda"))
