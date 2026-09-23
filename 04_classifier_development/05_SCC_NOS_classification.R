# ============================================================
# 05_SCC_NOS_classification.R
# Apply the final five-marker classifier to SCC-NOS and construct
# the five-state resolved SCC cohort.
# ============================================================

suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble);library(glmnet)})

MODEL <- "results/classifier/five_marker/five_marker_elastic_net_model.rds"
TRAIN <- "results/classifier/SCC_K_SCC_NK_training_cohort.csv"
SCC_NOS <- "data/processed/SCC_NOS_RPPA_classifier_data.csv"
OUTDIR <- "results/classifier/SCC_NOS"; dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
MARKERS <- c("ZAP-70","EVI1","EPPK1","ANNEXIN1","CD171")
NK_CUTOFF <- 0.35; K_CUTOFF <- 0.65

if(!file.exists(MODEL)) stop("Run 03_five_marker_classifier.R first.")
if(!file.exists(SCC_NOS)) stop("SCC-NOS prediction table not found: ",SCC_NOS)
obj<-readRDS(MODEL); model<-obj$model
nos<-readr::read_csv(SCC_NOS,show_col_types=FALSE)
if(!all(c("Patient",MARKERS)%in%names(nos))) stop("SCC-NOS table must contain Patient and all five markers.")
for(m in MARKERS) nos[[m]]<-as.numeric(nos[[m]])
if(anyNA(nos[,MARKERS])) stop("Missing values occur in SCC-NOS markers.")

x<-as.matrix(nos[,MARKERS]); x<-sweep(sweep(x,2,obj$scale_center,"-"),2,obj$scale_scale,"/")
p<-as.numeric(predict(model,newx=x,s="lambda.min",type="response"))

pred<-nos%>%mutate(keratinizing_probability=p,RPPA_keratinization_score=qlogis(pmin(pmax(p,1e-6),1-1e-6)),resolved_group=case_when(p<=NK_CUTOFF~"Nonkeratinizing-like",p>=K_CUTOFF~"Keratinizing-like",TRUE~"Indeterminate"))
summary<-pred%>%count(resolved_group)%>%mutate(percent=100*n/sum(n))

train<-readr::read_csv(TRAIN,show_col_types=FALSE)%>%mutate(resolved_group=ifelse(y==1,"Keratinizing","Nonkeratinizing"))
complete<-bind_rows(train%>%select(Patient,all_of(MARKERS),resolved_group),pred%>%select(Patient,all_of(MARKERS),resolved_group))

write_csv(pred,file.path(OUTDIR,"SCC_NOS_classification.csv"));write_csv(summary,file.path(OUTDIR,"SCC_NOS_state_summary.csv"));write_csv(complete,file.path(OUTDIR,"complete_resolved_SCC_cohort.csv"))
print(summary)
message("SCC-NOS classification complete. Cutoffs: <=",NK_CUTOFF," Nonkeratinizing-like; >=",K_CUTOFF," Keratinizing-like; intermediate = Indeterminate.")
