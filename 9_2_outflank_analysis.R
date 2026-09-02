#!/usr/bin/env Rscript
##############################################################################
# outflank_analysis.R
#
# Runs OutFLANK on a VCF to identify FST outlier SNPs.
#
# EDIT THE "USER SETTINGS" BLOCK BELOW BEFORE RUNNING.
##############################################################################

## ---- USER SETTINGS --------------------------------------------------------
vcf_path     <- "./phased_reheadered_filtered_sorted.vcf.gz"          # path to your VCF file
popmap_path  <- "popmap.txt"         # 2-column, no header: IndividualID  PopulationID
out_prefix   <- "outflank"           # prefix for output files
LeftTrimFrac  <- 0.05                # OutFLANK defaults - adjust if needed
RightTrimFrac <- 0.05
Hmin          <- 0.1
qthreshold    <- 0.05                # q-value threshold for calling outliers
## ---------------------------------------------------------------------------

## ---- Packages --------------------------------------------------------------

library(vcfR)
library(remotes)
library(OutFLANK)
library(qvalue)

## ---- Read data --------------------------------------------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_path)

geno <- extract.gt(vcf, element = "GT")

# Recode genotypes to OutFLANK format: 0, 1, 2, 9(missing)
recode_gt <- function(x) {
  x[x %in% c("0/0", "0|0")] <- 0
  x[x %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  x[x %in% c("1/1", "1|1")] <- 2
  x[!(x %in% c(0, 1, 2))] <- 9
  as.numeric(x)
}

G <- apply(geno, 2, recode_gt)

G <- t(G)  # rows = individuals, columns = loci

locusNames <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")
colnames(G) <- NULL  # OutFLANK doesn't need colnames on the matrix itself

popmap <- read.table(popmap_path, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("Ind", "Pop"))

# sanity check: sample order in G must match popmap order
stopifnot(nrow(G) == nrow(popmap))

## ---- Run OutFLANK --------------------------------------------------------
cat("Calculating FST matrix...\n")
fst_mat <- MakeDiploidFSTMat(G, locusNames = locusNames, popNames = popmap$Pop)

cat("Running OutFLANK...\n")
out <- OutFLANK(fst_mat,
                LeftTrimFraction  = LeftTrimFrac,
                RightTrimFraction = RightTrimFrac,
                Hmin = Hmin,
                NumberOfSamples = length(unique(popmap$Pop)),
                qthreshold = qthreshold)

## ---- Save diagnostic plot --------------------------------------------------
png(paste0(out_prefix, "_diagnostic.png"), width = 1000, height = 700)
OutFLANKResultsPlotter(out, withOutliers = TRUE, NoCorr = TRUE,
                       Hmin = Hmin, binwidth = 0.005, Zoom = FALSE,
                       RightZoomFraction = 0.05, titletext = NULL)
dev.off()

## ---- Save full results and outliers ---------------------------------------
results <- out$results
results$CHROM <- sub("_[0-9]+$", "", results$LocusName)
results$POS   <- as.numeric(sub("^.*_", "", results$LocusName))

write.table(results, paste0(out_prefix, "_full_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

outliers <- results[which(results$OutlierFlag == TRUE), ]
write.table(outliers, paste0(out_prefix, "_outliers.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also write outliers as a simple BED file for downstream intersection
# (BED is 0-based half-open, VCF POS is 1-based, so start = POS-1)
bed <- data.frame(chrom = outliers$CHROM,
                  start = outliers$POS - 1,
                  end   = outliers$POS,
                  name  = outliers$LocusName)
write.table(bed, paste0(out_prefix, "_outliers.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Done. Outliers written to:", paste0(out_prefix, "_outliers.txt"),
    "and", paste0(out_prefix, "_outliers.bed"), "\n")
