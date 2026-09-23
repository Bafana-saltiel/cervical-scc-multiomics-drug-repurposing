# ============================================================
# 03_RNA_consequence_analysis.R
# Proteomic consequences of gene-associated RNA abundance
# ============================================================

suppressPackageStartupMessages({library(limma);library(data.table);library(dplyr);library(ggplot2)})

INPUT_RDS <- "data/processed/CESC_multiomics_harmonised.rds"
ANNOTATION_FILE <- "data/processed/CESC_sample_annotations.csv"
OUTDIR <- "results/feature_discovery/RNA"
DRIVER_GENES <- c("CREBBP","DMD","KMT2D","LRP2","RYR2","SYNE1","USH2A")
STATE_FILTER <- "SCC_NOS"
LOW_QUANTILE <- 0.25
HIGH_QUANTILE <- 0.75
FDR_CUTOFF <- 0.05
EFFECT_CUTOFF <- 0.20

dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(INPUT_RDS); ann <- fread(ANNOTATION_FILE,data.table=FALSE)
stopifnot(all(c("Sample","State")%in%names(ann)))
get_component <- function(x,candidates){if(is.list(x)){hit<-candidates[candidates%in%names(x)];if(length(hit))return(x[[hit[1]]])};NULL}
rppa<-get_component(obj,c("RPPA","rppa","RPPA_matrix","rppa_matrix")); rna<-get_component(obj,c("RNA","rna","RNA_matrix","rna_matrix","logCPM"))
if(is.null(rppa)||is.null(rna)) stop("Harmonised RDS must contain RPPA and RNA components.")
rppa<-as.matrix(rppa);rna<-as.matrix(rna);if(!is.numeric(rppa))mode(rppa)<-"numeric";if(!is.numeric(rna))mode(rna)<-"numeric"
samples<-intersect(colnames(rppa),colnames(rna));if(!is.null(STATE_FILTER))samples<-intersect(samples,ann$Sample[ann$State==STATE_FILTER]);rppa<-rppa[,samples,drop=FALSE];rna<-rna[,samples,drop=FALSE]

results<-list()
for(gene in DRIVER_GENES){
 if(!gene%in%rownames(rna)) next
 x<-as.numeric(rna[gene,samples]);q<-quantile(x,probs=c(LOW_QUANTILE,HIGH_QUANTILE),na.rm=TRUE,names=FALSE)
 keep<-which(!is.na(x)&x<=q[1]|!is.na(x)&x>=q[2]);if(length(keep)<20)next
 grp<-ifelse(x[keep]>=q[2],1,0);if(sum(grp==0)<10||sum(grp==1)<10)next
 design<-model.matrix(~factor(grp,levels=c(0,1)));colnames(design)<-c("Intercept","RNAHigh")
 fit<-eBayes(lmFit(rppa[,keep,drop=FALSE],design));tt<-topTable(fit,coef="RNAHigh",number=Inf,sort.by="none")
 tt$protein<-rownames(tt);tt$driver_gene<-gene;tt$n_low<-sum(grp==0);tt$n_high<-sum(grp==1);tt$RPPA_difference<-tt$logFC;tt$gene_specific_FDR<-tt$adj.P.Val
 results[[gene]]<-tt[,c("driver_gene","protein","n_low","n_high","RPPA_difference","P.Value","gene_specific_FDR")]
}
res<-bind_rows(results);if(!nrow(res))stop("No RNA driver met the analysis requirements.")
res$global_FDR<-p.adjust(res$P.Value,method="BH");res$FDR_effect_supported<-res$gene_specific_FDR<FDR_CUTOFF&abs(res$RPPA_difference)>=EFFECT_CUTOFF;res$global_FDR_effect_supported<-res$global_FDR<FDR_CUTOFF&abs(res$RPPA_difference)>=EFFECT_CUTOFF
fwrite(res,file.path(OUTDIR,"RPPA_RNA_consequences_all.csv"));fwrite(filter(res,FDR_effect_supported),file.path(OUTDIR,"RPPA_RNA_consequences_gene_FDR_effect.csv"));fwrite(filter(res,global_FDR_effect_supported),file.path(OUTDIR,"RPPA_RNA_consequences_global_FDR_effect.csv"))
summary_tbl<-res%>%group_by(driver_gene)%>%summarise(n_low=first(n_low),n_high=first(n_high),nominal=sum(P.Value<0.05),gene_FDR=sum(gene_specific_FDR<FDR_CUTOFF),gene_FDR_effect=sum(FDR_effect_supported),global_FDR=sum(global_FDR<FDR_CUTOFF),global_FDR_effect=sum(global_FDR_effect_supported),.groups="drop")
fwrite(summary_tbl,file.path(OUTDIR,"RPPA_RNA_consequence_summary.csv"))
pdat<-res%>%group_by(driver_gene)%>%slice_min(global_FDR,n=10,with_ties=FALSE)%>%ungroup();p<-ggplot(pdat,aes(x=-log10(global_FDR),y=reorder(protein,-global_FDR)))+geom_point(aes(size=abs(RPPA_difference),shape=global_FDR_effect_supported))+facet_wrap(~driver_gene,scales="free_y")+theme_bw()+labs(x="-log10(global FDR)",y="RPPA feature",size="|RPPA difference|",shape="Global FDR + effect");ggsave(file.path(OUTDIR,"RPPA_RNA_consequences_top10.png"),p,width=12,height=9,dpi=600)
message("Completed RNA-consequence analysis: ",nrow(res)," tests.")
