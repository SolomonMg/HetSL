# TODO: Fit a super learned model
# 
# Author: solomon
###############################################################################
# setwd('~/Dropbox/creditClaimingProjects/het/HetPackage')
set.seed(23432)
  
load('tests/Het_Experiment.RData')

names(svdat)[23:51] = c("preq1", "preq2", "preq3", 
                      "nextc", "contr", "nextt", "treat",
                      "approval", "therm", "fiscRespbl", 
                      "bringMoneyEff", "passLegEff",
                      "secReqMC", "likGetM", "daysGetM", "break", 
                      "gender", "race", "byear", "ntvEnglish",
                      "ideo3", "voted", "pid3", "pidCloser", "educ",
                      "inc", "finalinst", "howLong", "comments")
approv = agrep("I pay attention", max.distance=.3,
             svdat$comments)
approv2 = agrep("I PAY ATTENTION", max.distance=.3,
              svdat$comments)
approv = c(approv,approv2)
#svdat$comments[-approv]
svdat = svdat[approv,]

svdat$cond.type[which(svdat$contr==1)] = "control" 
svdat$cond.type = relevel(factor(svdat$cond.type), ref="control") 
svdat$cond.money[which(svdat$contr==1)] = "control" 
svdat$cond.money = relevel(factor(svdat$cond.money), ref="control") 
svdat$cond.stage[which(svdat$contr==1)] = "control" 
svdat$cond.stage = relevel(factor(svdat$cond.stage), ref="control") 
svdat$cond.party[which(svdat$contr==1)] = "control" 
svdat$cond.party = relevel(factor(svdat$cond.party), ref="control") 
svdat$cond.alongWith[which(svdat$contr==1)] = "control" 
svdat$cond.alongWith = relevel(factor(svdat$cond.alongWith), ref="control") 
levels(svdat$cond.alongWith) = c("control", "alone", "w/ Dem", "w/ Rep") 

# Fix up pid3
svdat$pid3l = factor(c("Dem", "Rep", "Ind/Oth", "Ind/Oth")[svdat$pid3])
svdat$pid3l = relevel(svdat$pid3l, ref="Ind/Oth")

with<- rep(0, nrow(svdat))
with[grep('w/', as.character(svdat$cond.along))]<- 1

cons<- ifelse(svdat$ideo3<3, 1, 0)
lib<- ifelse(svdat$ideo3==4|svdat$ideo3==5, 1, 0)

##setting up the conditions
types<- sort(unique(as.character(svdat$cond.type)))
type.num<- match(svdat$cond.type, types)
number<- c('control', '$20 million', '$50 thousand')
amount.num<- match(svdat$cond.money, number)
request<- c('control', 'requested', 'secured', 'will request')
stage.num<- match(svdat$cond.stage, request)
party<- c('control', 'a Republican', 'a Democrat')
party.num<- match(svdat$cond.party, party)
along<- c('control', 'alone', 'w/ Rep', 'w/ Dem')
along.num<- match(svdat$cond.alongWith, along)

approve_bi<- ifelse(svdat$approval<3, 1, 0)
svdat$pid3 <- as.factor(svdat$pid3l)

# params for testing: 
formula = approve_bi ~ 0 + cond.type * pid3l + cond.party * pid3l
formula = approve_bi ~ 0 + cond.type + cond.party + pid3l
treatments = c("cond.type", "cond.party")
data = svdat
X_intercept = FALSE
bootstrap = FALSE
max_covariates = 1e10
SL.library = c("SL.glm", "SL.mean")
method = 'method.NNLS'
family = binomial()
id = NULL
verbose = TRUE
weights = NULL
control = list()
cvControl = list()


source('R/HetEffects.R')

# Models to implement: 
models <- c('lasso', 'e_net_0.75', 'e_net_0.5', 'e_net_0.25', 'FindIt', 
                  'BayesGLM', 'GLMBoost', 'BART', 'RandomForest', 'glm', 'KRLS', 'SVM-SMO')

# build Library 
SL.glmnet_a_1_lasso <- function(..., alpha = 1){ SL.glmnet(..., alpha = 1)}
SL.glmnet_a_0.75 <- function(..., alpha = 0.75){ SL.glmnet(..., alpha = 0.75)}
SL.glmnet_a_0.5 <- function(..., alpha = 0.5){ SL.glmnet(..., alpha = 0.5)}
SL.glmnet_a_0.25 <- function(..., alpha = 0.25){ SL.glmnet(..., alpha = 0.25)}

test.SL.library = c("SL.mean", "SL.glm", "SL.glmnet_a_1_lasso", "SL.glmnet_a_0.5"), 

# Set params for initial debugging
plotdat <- HetEffects(
  formula = approve_bi ~ 0 + cond.type * pid3l + cond.party * pid3l, 
  treatments = c("cond.type", "cond.party"), 
  data = svdat, 
  bootstrap = FALSE, 
  SL.library = test.SL.library, 
  method = 'method.NNLS', 
  family = binomial(), 
  id = NULL, 
  verbose = TRUE, 
  weights = NULL, 
  control = list(), 
  cvControl = list()
)

plotdat$weights


### CONVERT DUMMIES TO FACTOR ###
dummies = plotdat$effectX
header <- strsplit(colnames(dummies), ':')
header <- unlist(strsplit(colnames(dummies), '___'))[2 * (1:ncol(dummies))]
species <- factor(dummies %*% 1:ncol(dummies), labels = header)









n <- 500
p <- 2
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
#X <- data.frame(X)
T1 <- c("C", "T")[rbinom(n, 1, .5) +1]
T2 <- c("C", "T")[rbinom(n, 1, .5) +1]
Y <- as.numeric(as.factor(T1)) * X[, 1] + 
    as.numeric(as.factor(T2))  
data = data.frame(Y, T1, T2, X)
formula = Y ~ T1*T2*X
treatments = c('T1', 'T2')
weights = NULL
method = 'method.NNLS'
id = NULL 
verbose = FALSE
na.action = NULL
control = list()
cvControl = list()
family=gaussian()
R = 2

# build Library and run Super Learner
SL.glmnet_0.5 <- function(..., alpha = 0.5){ SL.glmnet(..., alpha = 0.5)}


SL.library <- c("SL.glm", "SL.mean")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.randomForest")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5", "SL.randomForest")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5", "SL.randomForest", "SL.gam")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights





##alright, beginning to 
SL.glmnet_0.5<- function(..., alpha = 0.5){ SL.glmnet(..., alpha = 0.5)}
ee<- SuperLearner(Y, X, SL.library = c("SL.glm", "SL.randomForest","SL.krls",  "SL.glmnet_0.5"))



n <- 500
p <- 50
t_p<- 3
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
X <- data.frame(X)
Treat<- matrix(NA, nrow=n, ncol=t_p)
for(z in 1:n){
	for(k in 1:3){
		Treat[z,k]<- rbinom(1, size = 1, p = 0.5)
		}
	}
	
Y <- pnorm(X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + 0.5*Treat[,1] + 0.25*Treat[,2] + -2*Treat[,3] + 5*X[,1]*Treat[,3]  + rnorm(n))

Y_bin<- c()
for(z in 1:length(Y)){
	Y_bin[z]<- rbinom(1, size =1, prob = Y[z])
	}

Y_bin[which(Y_bin==0)]<- -1	

mkstand<- function(x){
	x<- x - mean(x, na.rm=T)
	if(sd(x, na.rm=T)>0){
		x<- x/sd(x, na.rm=T)
		}
	return(x)
	}
	


SDsToRescaleX<- apply(X, 2, sd)
SDsToRescaleX[SDsToRescaleX==0]<- 1
Xstd<- apply(X, 2, mkstand)
colnames(Treat)<- paste('Treat', 1:3)
dd<- FindIt(Y_bin, Xstd, Treat, type='multiple', scale.c = SDsToRescaleX, search.lambdas=TRUE, fit.glmnet=TRUE)

ert<- cbind(1, Xstd, Treat)%*%dd$coefs/2 + 0.5

  treatments <- attr(Xmat, "treatments")
  Treat_columns = c()
  for(t in treatments){
    Treat_columns <- c(Treat_columns, grep(paste('^', treatments[i],sep=''),  colnames(Xmat))) 
  }
  Treat_columns = Treat_columns[!duplicated(Treat_columns)]
  attr(Xmat, "Treat_columns") <- Treat_columns


