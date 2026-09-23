# 01_RPPA_vs_RNA_classifier.R
# Orthogonal comparison of the fixed five-protein RPPA panel with RNA panels.
suppressPackageStartupMessages({library(readr);library(dplyr);library(pROC)})
OUT <- "results/validation"; dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
RPPA_FILE <- Sys.getenv("RPPA_CLASSIFIER_FILE","data/processed/labelled_SCC_RPPA_classifier_data.csv")
RNA_FILE <- Sys.getenv("RNA_VALIDATION_FILE","data/processed/SCC_RNA_validation_matrix.csv")
if(!file.exists(RPPA_FILE)) stop("Missing: ",RPPA_FILE); if(!file.exists(RNA_FILE)) stop("Missing: ",RNA_FILE)
rppa<-read_csv(RPPA_FILE,show_col_types=FALSE); rna<-read_csv(RNA_FILE,show_col_types=FALSE)
RPPA5<-c("ZAP-70","EVI1","EPPK1","ANNEXIN1","CD171"); RNA4<-c("KRT10","KRT13","KRT16","IVL")
if(!all(c("Patient","classifier_group",RPPA5)%in%names(rppa))) stop("RPPA file must contain Patient, classifier_group and the five RPPA markers.")
if(!all(c("Patient",RPPA5,RNA4)%in%names(rna))) stop("RNA file must contain Patient, five marker genes and KRT10/KRT13/KRT16/IVL.")
rppa<-rppa%>%filter(classifier_group%in%c("SCC-K","SCC-NK","Keratinizing","Nonkeratinizing"))%>%mutate(label=as.integer(classifier_group%in%c("SCC-K","Keratinizing")))
d<-inner_join(rppa%>%select(Patient,label,all_of(RPPA5)),rna%>%select(Patient,all_of(c(RPPA5,RNA4))),by="Patient")
if(nrow(d)<10) stop("Too few matched patients: ",nrow(d))
fit<-function(dat,features,name){cc<-complete.cases(dat[,features]); x<-dat[cc,,drop=FALSE]; m<-glm(label~.,data=x[,c("label",features)],family=binomial()); s<-as.numeric(predict(m,newdata=x,type="response")); ro<-roc(x$label,s,quiet=TRUE,direction="<"); ci<-ci.auc(ro); tibble(Patient=x$Patient,label=x$label,model=name,score=s,AUC=as.numeric(auc(ro)),CI_low=ci[1],CI_high=ci[3])}
res<-bind_rows(fit(d,RPPA5,"RPPA_5_protein"),fit(d,RPPA5,"RNA_5_gene"),fit(d,RNA4,"RNA_conventional_4_gene"))
perf<-res%>%group_by(model)%>%summarise(n=n(),AUC=first(AUC),CI_low=first(CI_low),CI_high=first(CI_high),.groups="drop")
write_csv(res,file.path(OUT,"RPPA_vs_RNA_classifier_predictions.csv"));write_csv(perf,file.path(OUT,"RPPA_vs_RNA_classifier_performance.csv"));print(perf)
