suppressPackageStartupMessages({
  library(DiffBind)
  library(GenomicRanges)
  library(rtracklayer)
})

#load in the Arguments
args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
sample_sheet <- args[2]
cofactor <- args[3]

# 1) load
dbaObj <- dba(sampleSheet = sample_sheet)

# 2) count reads in consensus peaks
dbaObj <- dba.count(dbaObj, summits = 50)  #2x50bp standard window

info <- dba.show(dbaObj)
libsizes <- cbind(LibReads=info$Reads, FRIP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$SampleID
cat("\nLibrary sizes: \n")
print(libsizes)

# 3) Normalise the data
norm <- dba.normalize(dbaObj)
cat("\nNormalised library sizes: \n")
print(dba.show(dbaObj))

normlibsizes <- cbind(
  FullLibSize=norm$lib.sizes, 
  NormFactor=norm$norm.factors, 
  NormLibSize=round(norm$lib.sizes * norm$norm.factors)
)
rownames(normlibsizes) <- info$SampleID
cat("\nNormalised library sizes: \n")
print(normlibsizes)

# 4) model: ~ cofactor + Condition 
if (cofactor != "") {
  design_formula <- paste("~", cofactor, "+ Condition", sep = "")
  dbaObj <- dba.contrast(
    dbaObj,
    design =  design_formula,
    reorderMeta=list(Condition="Hsp")
  )
} else {
  dbaObj <- dba.contrast(
    dbaObj,
    reorderMeta=list(Condition="Hsp")
  )
}

cat("\nContrasts design: \n")
print(dbaObj)

# 4) analyse (DESeq2 backend)
dbaObj <- dba.analyze(dbaObj, method = DBA_DESEQ2)

# 5) report Condition effect (Hsp vs C)
res <- dba.report(dbaObj, method = DBA_DESEQ2, th = 0.05, bCalled = TRUE)
tbl <- as.data.frame(res)
write.csv(tbl, file.path(outdir, "diffbind_Hsp_vs_C.csv"), row.names = FALSE)
cat("\nSignificant peaks (FDR <= 0.05):\n")
print(tbl)

# 6) export significant peaks as BED (FDR <= 0.05)
sig <- res[res$FDR <= 0.05]
if (length(sig) > 0) {
  export(sig, con = file.path(outdir, "diffbind_Hsp_vs_C.FDR05.bed"), format = "BED")
}

# Correlations
corr <- try(dba.plotHeatmap(dbaObj, correlations = TRUE, contrast = 1, bRetrieve = TRUE), silent = TRUE)
if (!inherits(corr, "try-error")) {
  write.csv(corr, file.path(outdir, "sample_correlations.csv"))
}
qc_post <- capture.output(dba.show(dbaObj))
writeLines(qc_post, file.path(outdir, "dba_show_after_analyze.txt"))

# 7) save plots: one multi-page PDF + individual PNGs
pdf_file <- file.path(outdir, "qc_plots.pdf")

pdf(pdf_file, width = 8, height = 6, onefile = TRUE, useDingbats = FALSE)

# Each plot in try() so a failure doesn’t abort the run
try(plot(dbaObj), silent = TRUE)                      # correlation heatmap / PCA
try(dba.plotMA(dbaObj, contrast = 1), silent = TRUE)
try(dba.plotVolcano(dbaObj, contrast = 1), silent = TRUE)
try(dba.plotBox(dbaObj), silent = TRUE)
try(dba.plotVenn(dbaObj, contrast = 1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE), silent = TRUE)
try(dba.plotPCA(dbaObj, label = DBA_CONDITION), silent = TRUE)
try(dba.plotHeatmap(dbaObj, contrast = 1, correlations = FALSE), silent = TRUE)
dev.off()