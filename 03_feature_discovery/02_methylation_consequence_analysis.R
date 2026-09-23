# ============================================================
# 02_methylation_consequence_analysis.R
# Proteomic consequences of gene-associated methylation states
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ------------------------- CONFIG ----------------------------
INPUT_RDS <- "data/processed/CESC_multiomics_harmonised.rds"
ANNOTATION_FILE <- "data/processed/CESC_sample_annotations.csv"
OUTDIR <- "results/feature_discovery/Methylation"
DRIVER_GENES <- c("CREBBP", "DMD", "KMT2D", "LRP2", "RYR2", "SYNE1", "USH2A")
STATE_FILTER <- "SCC_NOS"
METHYLATION_FEATURE_MAP <- "data/processed/CESC_gene_methylation_features.csv"
LOW_QUANTILE <- 0.25
HIGH_QUANTILE <- 0.75
FDR_CUTOFF <- 0.05
EFFECT_CUTOFF <- 0.20
# -------------------------------------------------------------

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(INPUT_RDS)
ann <- fread(ANNOTATION_FILE, data.table = FALSE)
stopifnot(all(c("Sample", "State") %in% names(ann)))

get_component <- function(x, candidates) {
  if (is.list(x)) {
    hit <- candidates[candidates %in% names(x)]
    if (length(hit)) return(x[[hit[1]]])
  }
  NULL
}
rppa <- get_component(obj, c("RPPA","rppa","RPPA_matrix","rppa_matrix"))
meth <- get_component(obj, c("Methylation","methylation","meth","methylation_matrix"))
if (is.null(rppa) || is.null(meth)) stop("Harmonised RDS must contain RPPA and Methylation components.")
rppa <- as.matrix(rppa); meth <- as.matrix(meth)
if (!is.numeric(rppa)) mode(rppa) <- "numeric"
if (!is.numeric(meth)) mode(meth) <- "numeric"

samples <- intersect(colnames(rppa), colnames(meth))
if (!is.null(STATE_FILTER)) samples <- intersect(samples, ann$Sample[ann$State == STATE_FILTER])
rppa <- rppa[, samples, drop=FALSE]; meth <- meth[, samples, drop=FALSE]

# Required mapping: one gene-associated methylation feature per driver gene.
# Accepted columns are DriverGene and MethylationFeature. This avoids silently
# choosing an arbitrary CpG when several probes map to the same gene.
if (!file.exists(METHYLATION_FEATURE_MAP)) stop("Create ", METHYLATION_FEATURE_MAP, " with columns DriverGene and MethylationFeature.")
map <- fread(METHYLATION_FEATURE_MAP, data.table=FALSE)
stopifnot(all(c("DriverGene","MethylationFeature") %in% names(map)))

results <- list()
for (gene in DRIVER_GENES) {
  feat <- map$MethylationFeature[match(gene, map$DriverGene)]
  if (is.na(feat) || !feat %in% rownames(meth)) next
  x <- as.numeric(meth[feat, samples])
  q <- quantile(x, probs=c(LOW_QUANTILE,HIGH_QUANTILE), na.rm=TRUE, names=FALSE)
  keep <- which(!is.na(x) & x <= q[1] | !is.na(x) & x >= q[2])
  if (length(keep) < 20) next
  grp <- ifelse(x[keep] >= q[2], 1, 0)
  if (sum(grp==0) < 10 || sum(grp==1) < 10) next
  design <- model.matrix(~ factor(grp, levels=c(0,1)))
  colnames(design) <- c("Intercept","MethylationHigh")
  fit <- eBayes(lmFit(rppa[, keep, drop=FALSE], design))
  tt <- topTable(fit, coef="MethylationHigh", number=Inf, sort.by="none")
  tt$protein <- rownames(tt)
  tt$driver_gene <- gene
  tt$methylation_feature <- feat
  tt$n_low <- sum(grp==0); tt$n_high <- sum(grp==1)
  tt$RPPA_difference <- tt$logFC
  tt$gene_specific_FDR <- tt$adj.P.Val
  results[[gene]] <- tt[,c("driver_gene","methylation_feature","protein","n_low","n_high","RPPA_difference","P.Value","gene_specific_FDR")]
}
res <- bind_rows(results)
if (!nrow(res)) stop("No methylation driver met the analysis requirements.")
res$global_FDR <- p.adjust(res$P.Value, method="BH")
res$FDR_effect_supported <- res$gene_specific_FDR < FDR_CUTOFF & abs(res$RPPA_difference) >= EFFECT_CUTOFF
res$global_FDR_effect_supported <- res$global_FDR < FDR_CUTOFF & abs(res$RPPA_difference) >= EFFECT_CUTOFF
fwrite(res, file.path(OUTDIR,"RPPA_methylation_consequences_all.csv"))
fwrite(filter(res,FDR_effect_supported), file.path(OUTDIR,"RPPA_methylation_consequences_gene_FDR_effect.csv"))
fwrite(filter(res,global_FDR_effect_supported), file.path(OUTDIR,"RPPA_methylation_consequences_global_FDR_effect.csv"))

summary_tbl <- res %>% group_by(driver_gene) %>% summarise(
 n_low=first(n_low), n_high=first(n_high), nominal=sum(P.Value<0.05),
 gene_FDR=sum(gene_specific_FDR<FDR_CUTOFF), gene_FDR_effect=sum(FDR_effect_supported),
 global_FDR=sum(global_FDR<FDR_CUTOFF), global_FDR_effect=sum(global_FDR_effect_supported), .groups="drop")
fwrite(summary_tbl,file.path(OUTDIR,"RPPA_methylation_consequence_summary.csv"))

pdat <- res %>% group_by(driver_gene) %>% slice_min(global_FDR,n=10,with_ties=FALSE) %>% ungroup()
p <- ggplot(pdat,aes(x=-log10(global_FDR),y=reorder(protein,-global_FDR)))+
 geom_point(aes(size=abs(RPPA_difference),shape=global_FDR_effect_supported))+
 facet_wrap(~driver_gene,scales="free_y")+theme_bw()+
 labs(x="-log10(global FDR)",y="RPPA feature",size="|RPPA difference|",shape="Global FDR + effect")
ggsave(file.path(OUTDIR,"RPPA_methylation_consequences_top10.png"),p,width=12,height=9,dpi=600)
message("Completed methylation-consequence analysis: ",nrow(res)," tests.")
