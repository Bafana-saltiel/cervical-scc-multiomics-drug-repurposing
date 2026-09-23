# ============================================================
# 03_five_marker_classifier.R
# Fit the final five-marker elastic-net classifier and classify
# pathology-labelled training tumours.
# ============================================================

suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble);library(glmnet)})

SEED <- 20260916
set.seed(SEED)
INPUT <- "results/classifier/SCC_K_SCC_NK_training_cohort.csv"
OUTDIR <- "results/classifier/five_marker"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
MARKERS <- c("ZAP-70","EVI1","EPPK1","ANNEXIN1","CD171")
ALPHA_GRID <- c(0,0.25,0.5,0.75,1)

if(!file.exists(INPUT)) stop("Run script 01 first.")
dat <- readr::read_csv(INPUT,show_col_types=FALSE)
x <- scale(as.matrix(dat[,MARKERS]))
y <- dat$y

# Internal tuning is performed within the labelled SCC-K/SCC-NK cohort.
fold_id <- integer(length(y))
for(cls in sort(unique(y))){ii<-which(y==cls);fold_id[ii]<-sample(rep(1:min(5,length(ii)),length.out=length(ii)))}

fits <- lapply(ALPHA_GRID,function(a) glmnet::cv.glmnet(x,y,family="binomial",alpha=a,foldid=fold_id,type.measure="deviance",standardize=FALSE,nlambda=100,maxit=100000))
metrics <- bind_rows(lapply(seq_along(fits),function(i){f<-fits[[i]];j<-which.min(abs(f$lambda-f$lambda.min));tibble(alpha=ALPHA_GRID[i],lambda=f$lambda.min,cv_deviance=f$cvm[j],cv_se=f$cvsd[j])})) %>% arrange(cv_deviance,cv_se,desc(alpha))
best <- metrics$alpha[1]
fit <- fits[[which(ALPHA_GRID==best)[1]]]

coef <- as.matrix(coef(fit,s="lambda.min"))
coef_tbl <- tibble(feature=rownames(coef),coefficient=as.numeric(coef[,1]),retained=as.numeric(coef[,1])!=0)

prob <- as.numeric(predict(fit,newx=x,s="lambda.min",type="response"))
training_predictions <- dat %>% mutate(keratinizing_probability=prob,RPPA_keratinization_score=qlogis(pmin(pmax(prob,1e-6),1-1e-6)),predicted_group=ifelse(prob>=0.5,"SCC-K","SCC-NK"))

saveRDS(list(model=fit,markers=MARKERS,alpha=best,lambda=fit$lambda.min,scale_center=attr(x,"scaled:center"),scale_scale=attr(x,"scaled:scale")),file.path(OUTDIR,"five_marker_elastic_net_model.rds"))
write_csv(metrics,file.path(OUTDIR,"alpha_tuning.csv")); write_csv(coef_tbl,file.path(OUTDIR,"final_coefficients.csv")); write_csv(training_predictions,file.path(OUTDIR,"training_predictions.csv"))

message("Final five-marker elastic-net model saved. Markers: ",paste(MARKERS,collapse=", "))
