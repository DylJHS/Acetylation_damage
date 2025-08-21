suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
  library(rtracklayer)
})

outdir <- "/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/tagged_alignments/aligned_bams/diffbind"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1) load
dbaObj <- dba(sampleSheet = file.path("/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/data/inputs/diffbind_samplesheet_3.csv"))

# 2) count reads in consensus peaks
dbaObj <- dba.count(dbaObj, summits = 50)  #2x50bp standard window

info <- dba.show(dbaObj)
libsizes <- cbind(LibReads=info$Reads, FRIP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$SampleID
cat("\nLibrary sizes: \n")
print(libsizes)

# 3) model: ~ Conc + Condition 
dbaObj <- dba.contrast(
  dbaObj,
  design="~Antibody + Condition",
  reorderMeta=list(Condition="Hsp")
)

cat("\nContrasts design: \n")
print(dbaObj)

# 4) analyse (DESeq2 backend)
dbaObj <- dba.analyze(dbaObj, method = DBA_DESEQ2)

# 5) report Condition effect (Hsp vs C)
res <- dba.report(dbaObj, method = DBA_DESEQ2, th = 0.05, bCalled = TRUE)
tbl <- as.data.frame(res)
write.csv(tbl, file.path(outdir, "diffbind_Hsp_vs_C_cofactor.csv"), row.names = FALSE)
cat("\nSignificant peaks (FDR <= 0.05):\n")
print(tbl)

# 6) export significant peaks as BED (FDR <= 0.05)
sig <- res[res$FDR <= 0.05]
if (length(sig) > 0) {
  export(sig, con = file.path(outdir, "diffbind_Hsp_vs_C_cofactor.FDR05.bed"), format = "BED")
}

# Correlations
corr <- try(dba.plotHeatmap(dbaObj, correlations = TRUE, contrast = 1, bRetrieve = TRUE), silent = TRUE)
if (!inherits(corr, "try-error")) {
  write.csv(corr, file.path(outdir, "sample_correlations_cofactor.csv"))
}
qc_post <- capture.output(dba.show(dbaObj))
writeLines(qc_post, file.path(outdir, "dba_show_after_analyze_cofactor.txt"))

# 7) save plots: one multi-page PDF + individual PNGs
pdf_file <- file.path(outdir, "qc_plots_cofactor.pdf")

pdf(pdf_file, width = 8, height = 6, onefile = TRUE, useDingbats = FALSE)

# Each plot in try() so a failure doesn’t abort the run
try(plot(dbaObj), silent = TRUE)                      # correlation heatmap / PCA
try(dba.plotMA(dbaObj, contrast = 1), silent = TRUE)
try(dba.plotVolcano(dbaObj, contrast = 1), silent = TRUE)
try(dba.plotPCA(dbaObj, label = DBA_CONDITION), silent = TRUE)
try(dba.plotHeatmap(dbaObj, contrast = 1, correlations = FALSE), silent = TRUE)
dev.off()