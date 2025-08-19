suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
  library(rtracklayer)
})

outdir <- "/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/tagged_alignments/aligned_bams/diffbind"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1) load
dbaObj <- dba(sampleSheet = file.path("/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/data/inputs/diffbind_samplesheet.csv"))

# 2) count reads in consensus peaks
dbaObj <- dba.count(dbaObj, summits = 50)  #2x50bp standard window

# 3) model: ~ Conc + Condition 
use_design <- "design" %in% names(formals(dba.contrast))

dbaObj <- dba.contrast(
    dbaObj,
    design = "~ Conc + Condition",
    reorderMeta = list(Condition = c("C","Hsp"))
)

# 4) analyse (DESeq2 backend)
dbaObj <- dba.analyze(dbaObj, method = DBA_DESEQ2)

# 5) report Condition effect (Hsp vs C)
res <- dba.report(dbaObj, method = DBA_DESEQ2, th = 0.05, bCalled = TRUE)
tbl <- as.data.frame(res)
write.csv(tbl, file.path(outdir, "diffbind_Hsp_vs_C.csv"), row.names = FALSE)

# 6) export significant peaks as BED (FDR <= 0.05)
sig <- res[res$FDR <= 0.05]
if (length(sig) > 0) {
  export(sig, con = file.path(outdir, "diffbind_Hsp_vs_C.FDR05.bed"), format = "BED")
}

# 7) some quick plots
pdf(file.path(outdir, "qc_plots.pdf"))
plot(dbaObj)                    # PCA / Correlation heatmap
dba.plotMA(dbaObj, contrast=1)  # MA-plot for Hsp vs C
dba.plotVolcano(dbaObj, contrast=1)
dev.off()