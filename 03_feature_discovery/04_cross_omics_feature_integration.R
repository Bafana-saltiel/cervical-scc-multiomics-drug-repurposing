# ============================================================
# 04_cross_omics_feature_integration.R
# Integrate significant RPPA consequences across alteration types
# ============================================================

suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(ggplot2)})
BASE <- "results/feature_discovery"
OUTDIR <- file.path(BASE,"Integrated")
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

read_supported <- function(path, alteration){
 if(!file.exists(path)) return(NULL)
 x<-fread(path,data.table=FALSE);if(!nrow(x))return(NULL)
 x$alteration_type<-alteration; x
}
files <- list(
 Mutation=file.path(BASE,"RPPA/RPPA_mutation_consequences_gene_FDR_effect.csv"),
 Methylation=file.path(BASE,"Methylation/RPPA_methylation_consequences_gene_FDR_effect.csv"),
 RNA=file.path(BASE,"RNA/RPPA_RNA_consequences_gene_FDR_effect.csv"))
res<-bind_rows(Map(read_supported,files,names(files)))
if(!nrow(res))stop("No supported consequence tables found. Run scripts 01-03 first.")

fwrite(res,file.path(OUTDIR,"cross_omics_supported_RPPA_features.csv"))

# Protein-by-driver recurrence across alteration types.
recurrence<-res%>%count(driver_gene,protein,alteration_type,name="n")%>%distinct(driver_gene,protein,.keep_all=TRUE)%>%group_by(driver_gene,protein)%>%summarise(n_alteration_types=n_distinct(alteration_type),alteration_types=paste(sort(unique(alteration_type)),collapse=";"),.groups="drop")%>%arrange(desc(n_alteration_types),driver_gene,protein)
fwrite(recurrence,file.path(OUTDIR,"cross_omics_protein_recurrence.csv"))

# 21-classifier-style feature matrix when all seven genes and three alteration types are represented.
wide<-res%>%select(driver_gene,alteration_type,protein,RPPA_difference)%>%mutate(classifier=paste(driver_gene,alteration_type,sep="__"))%>%select(protein,classifier,RPPA_difference)%>%pivot_wider(names_from=classifier,values_from=RPPA_difference,values_fill=0)
fwrite(wide,file.path(OUTDIR,"cross_omics_RPPA_feature_matrix.csv"))

# Top recurrent proteins across the supported signatures.
plot_dat<-recurrence%>%filter(n_alteration_types>=1)%>%group_by(n_alteration_types)%>%slice_head(n=15)%>%ungroup()
p<-ggplot(plot_dat,aes(x=n_alteration_types,y=reorder(paste(driver_gene,protein,sep=": "),n_alteration_types)))+geom_point()+theme_bw()+labs(x="Number of alteration types",y="Driver: protein")
ggsave(file.path(OUTDIR,"cross_omics_feature_recurrence.png"),p,width=9,height=8,dpi=600)

message("Integrated ",nrow(res)," supported RPPA consequence records across alteration types.")
