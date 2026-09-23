# ============================================================
# 04_nested_cross_validation.R
# Repeated nested cross-validation of the five-marker classifier
# ============================================================

suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble);library(glmnet);library(pROC)})

SEED <- 20260916; set.seed(SEED)
INPUT <- "results/classifier/SCC_K_SCC_NK_training_cohort.csv"
OUTDIR <- "results/classifier/nested_cv"; dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
MARKERS <- c("ZAP-70","EVI1","EPPK1","ANNEXIN1","CD171")
ALPHA_GRID <- c(0,0.25,0.5,0.75,1); N_REPEATS <- 20L; N_OUTER <- 5L

dat <- readr::read_csv(INPUT,show_col_types=FALSE); x0<-as.matrix(dat[,MARKERS]); y<-dat$y
stratified <- function(y,k){z<-integer(length(y));for(cls in sort(unique(y))){ii<-sample(which(y==cls));z[ii]<-rep(seq_len(k),length.out=length(ii))};z}
auc_safe <- function(y,p) as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE)))

preds<-list(); idx<-0L
for(r in seq_len(N_REPEATS)){
 set.seed(SEED+r); folds<-stratified(y,N_OUTER)
 for(f in seq_len(N_OUTER)){
  tr<-which(folds!=f); va<-which(folds==f)
  center<-colMeans(x0[tr,,drop=FALSE]); sc<-apply(x0[tr,,drop=FALSE],2,sd); sc[!is.finite(sc)|sc==0]<-1
  xt<-sweep(sweep(x0[tr,,drop=FALSE],2,center,"-"),2,sc,"/"); xv<-sweep(sweep(x0[va,,drop=FALSE],2,center,"-"),2,sc,"/")
  inner<-stratified(y[tr],min(5,min(table(y[tr])))); fits<-lapply(ALPHA_GRID,function(a)glmnet::cv.glmnet(xt,y[tr],family="binomial",alpha=a,foldid=inner,type.measure="deviance",standardize=FALSE,nlambda=100,maxit=100000))
  cv<-sapply(fits,function(z)min(z$cvm)); ai<-which.min(cv); fit<-fits[[ai]]
  p<-as.numeric(predict(fit,newx=xv,s="lambda.min",type="response")); idx<-idx+1L
  preds[[idx]]<-tibble(repeat_id=r,outer_fold=f,Patient=dat$Patient[va],observed_binary=y[va],observed_group=as.character(dat$classifier_group[va]),keratinizing_probability=p,selected_alpha=ALPHA_GRID[ai])
 }
}

pred<-bind_rows(preds)
patient_oof<-pred%>%group_by(Patient,observed_binary,observed_group)%>%summarise(keratinizing_probability=mean(keratinizing_probability),.groups="drop")
roc<-pROC::roc(patient_oof$observed_binary,patient_oof$keratinizing_probability,levels=c(0,1),direction="<",quiet=TRUE); auc<-as.numeric(pROC::auc(roc)); ci<-as.numeric(pROC::ci.auc(roc,method="delong"))
perf<-tibble(AUC=auc,AUC_CI_low=ci[1],AUC_CI_high=ci[3],repeats=N_REPEATS,outer_folds=N_OUTER,n=nrow(patient_oof))
write_csv(pred,file.path(OUTDIR,"nested_cv_predictions.csv"));write_csv(patient_oof,file.path(OUTDIR,"patient_oof_predictions.csv"));write_csv(perf,file.path(OUTDIR,"nested_cv_performance.csv"))
print(perf)
