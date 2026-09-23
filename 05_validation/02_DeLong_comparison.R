# 02_DeLong_comparison.R
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(pROC)})
IN<-"results/validation/RPPA_vs_RNA_classifier_predictions.csv";OUT<-"results/validation";if(!file.exists(IN))stop("Run 01_RPPA_vs_RNA_classifier.R first.")
p<-read_csv(IN,show_col_types=FALSE);w<-p%>%select(Patient,label,model,score)%>%pivot_wider(names_from=model,values_from=score)
cmp<-function(a,b){z<-w%>%filter(complete.cases(.data[[a]],.data[[b]],label));r1<-roc(z$label,z[[a]],quiet=TRUE,direction="<");r2<-roc(z$label,z[[b]],quiet=TRUE,direction="<");t<-roc.test(r1,r2,paired=TRUE,method="delong");c1<-ci.auc(r1);c2<-ci.auc(r2);tibble(model_1=a,model_2=b,n=nrow(z),AUC_1=auc(r1),CI1_low=c1[1],CI1_high=c1[3],AUC_2=auc(r2),CI2_low=c2[1],CI2_high=c2[3],Delta_AUC=auc(r1)-auc(r2),Z=as.numeric(t$statistic),P_value=t$p.value)}
out<-bind_rows(cmp("RPPA_5_protein","RNA_conventional_4_gene"),cmp("RPPA_5_protein","RNA_5_gene"));write_csv(out,file.path(OUT,"RPPA_vs_RNA_DeLong_test.csv"));print(out)
