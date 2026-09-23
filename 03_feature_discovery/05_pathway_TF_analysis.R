# ============================================================
# 05_pathway_TF_analysis.R
# Reactome pathway and transcription-factor enrichment of RPPA signatures
# ============================================================

suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(ggplot2)})

INPUT <- "results/feature_discovery/Integrated/cross_omics_RPPA_feature_matrix.csv"
OUTDIR <- "results/feature_discovery/Pathway_TF"
SIGNATURE_FDR <- 0.05
MIN_SIGNATURE_SIZE <- 3
ENRICHR_WEBSITE <- "https://maayanlab.cloud/Enrichr/"

dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(INPUT)) stop("Run 04_cross_omics_feature_integration.R first.")
mat<-fread(INPUT,data.table=FALSE)
if(!"protein"%in%names(mat))stop("Feature matrix must contain protein column.")

# Map RPPA labels to gene symbols using an optional local mapping. This is
# deliberately explicit because RPPA labels can contain phosphosite/protein aliases.
MAP_FILE <- "data/processed/RPPA_feature_to_gene_symbol.csv"
if(!file.exists(MAP_FILE)){
 message("No RPPA_feature_to_gene_symbol.csv found. Writing a template and stopping before enrichment.")
 fwrite(data.frame(RPPA_feature=character(),gene_symbol=character()),MAP_FILE)
 stop("Populate ",MAP_FILE," with RPPA feature -> gene_symbol mappings, then rerun.")
}
map<-fread(MAP_FILE,data.table=FALSE);stopifnot(all(c("RPPA_feature","gene_symbol")%in%names(map)))
mat<-mat%>%left_join(map,by=c("protein"="RPPA_feature"))

# Build one gene list per alteration-defined RPPA signature.
classifier_cols<-setdiff(names(mat),c("protein","gene_symbol"))
signatures<-lapply(classifier_cols,function(cl){
 x<-mat[[cl]];genes<-mat$gene_symbol[!is.na(x)&x!=0&!is.na(mat$gene_symbol)]
 unique(genes)
});names(signatures)<-classifier_cols
signatures<-signatures[lengths(signatures)>=MIN_SIGNATURE_SIZE]

sig_table<-bind_rows(lapply(names(signatures),function(n)data.frame(signature=n,gene_symbol=signatures[[n]])))
fwrite(sig_table,file.path(OUTDIR,"RPPA_signature_gene_lists.csv"))

# Enrichment is performed through Enrichr-compatible gene lists. The script
# writes upload-ready files rather than silently depending on an online API.
for(n in names(signatures)) writeLines(signatures[[n]],file.path(OUTDIR,paste0("genes_",gsub("[^A-Za-z0-9]+","_",n),".txt")))

# Optional local enrichment input: users can place exported Enrichr results
# in results/feature_discovery/Pathway_TF/enrichr/ with columns
# signature,term,adjusted_p,combined_score,library.
ENRICH_DIR<-file.path(OUTDIR,"enrichr")
dir.create(ENRICH_DIR,recursive=TRUE,showWarnings=FALSE)
exported<-list.files(ENRICH_DIR,pattern="\\.(csv|tsv)$",full.names=TRUE)
if(length(exported)){
 enr<-bind_rows(lapply(exported,function(f){x<-fread(f,data.table=FALSE);x$source_file<-basename(f);x}))
 fwrite(enr,file.path(OUTDIR,"enrichment_results_combined.csv"))
 if(all(c("signature","term","adjusted_p")%in%names(enr))){
  pdat<-enr%>%filter(adjusted_p<SIGNATURE_FDR)%>%group_by(signature)%>%slice_min(adjusted_p,n=10,with_ties=FALSE)%>%ungroup()
  if(nrow(pdat)){
   p<-ggplot(pdat,aes(x=-log10(adjusted_p),y=reorder(term,-adjusted_p)))+geom_point(aes(size=combined_score))+facet_wrap(~signature,scales="free_y")+theme_bw()+labs(x="-log10 adjusted P",y="Enriched term",size="Combined score")
   ggsave(file.path(OUTDIR,"pathway_TF_enrichment_top10.png"),p,width=14,height=10,dpi=600)
  }
 }
} else {
 message("No exported Enrichr tables found. Gene lists were generated for Reactome/ChEA/ENCODE/TRRUST analysis.")
}

writeLines(c(
 "Enrichment resources used in the study design:",
 "- Reactome pathway enrichment for RPPA signature proteins.",
 "- Consolidated transcription-factor resources including ChEA, ENCODE and TRRUST (via Enrichr or equivalent export).",
 "Interpretation: enrichment identifies representation of pathways/regulators in the protein signatures; it does not establish pathway activation, TF binding, directionality or causality."
),file.path(OUTDIR,"README_enrichment.txt"))

message("Prepared ",length(signatures)," RPPA signatures for pathway/TF enrichment.")
