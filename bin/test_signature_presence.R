#!/opt/conda/bin/Rscript

library(mSigAct)

args <- commandArgs(trailingOnly = TRUE)

spectra_file <- args[1]
catalog_file <- args[2]
target_sig   <- args[3]
cpus         <- as.integer(args[4])
output_name  <- args[5]

df_spectra <- read.table(
    spectra_file,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)

df_catalog <- read.table(
    catalog_file,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
rownames(df_spectra) <- gsub("[^A-Za-z]", "", rownames(df_spectra))

rownames(df_catalog) <- gsub("[^A-Za-z]", "", rownames(df_catalog))

common_rows <- intersect(rownames(df_spectra),rownames(df_catalog))

if (length(common_rows) != 96) {
    stop(
        "Row names do not match between df_spectra and df_catalog."
    )
}

spectra_mat <- as.matrix(df_spectra)
mode(spectra_mat) <- "numeric"

sigs_mat <- as.matrix(df_catalog)
mode(sigs_mat) <- "numeric"

if (nrow(spectra_mat) != nrow(sigs_mat)) {
    stop(paste("Mismatch: spectra has", nrow(spectra_mat), "rows, but sigs has", nrow(sigs_mat), "rows."))
}

sig.presence.test.out <- SignaturePresenceTest(
    spectra = spectra_mat,
    sigs = sigs_mat,
    target.sig.index = target_sig,
    mc.cores = cpus
)

save(
    sig.presence.test.out,
    file = file.path(
        output_name
    )
)