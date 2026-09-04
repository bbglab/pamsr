#!/opt/conda/bin/Rscript

library(mSigAct)

args <- commandArgs(trailingOnly = TRUE)

spectra_file <- args[1]
catalog_file <- args[2]
target_sig   <- args[3]

raw_bg_sigs <- args[4]
clean_bg_sigs <- gsub("[\\[\\]\"\\' ]", "", raw_bg_sigs)
background_sigs <- unlist(strsplit(clean_bg_sigs, ","))

cpus         <- as.integer(args[5])
output_name  <- args[6]

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

# Keep only samples/signatures specified in background_sigs
missing_sigs <- setdiff(background_sigs, colnames(sigs_mat))

if (length(missing_sigs) > 0) {
    stop(
        paste(
            "The following background signatures are not present in sigs_mat:",
            paste(missing_sigs, collapse = ", ")
        )
    )
}

# FIX MATRIX SUBSETTING
bg_sigs_mat <- sigs_mat[, background_sigs, drop = FALSE]
# print("Columns in sigs_mat:")
# print(colnames(sigs_mat))
# print("Columns in bg_sigs_mat:")
# print(colnames(bg_sigs_mat))

# print(bg_sigs_mat)
# print(spectra_mat)

sig.presence.test.out <- SignaturePresenceTest(
    spectra = spectra_mat,
    sigs = bg_sigs_mat, 
    target.sig.index = target_sig,
    mc.cores = cpus
)

save(
    sig.presence.test.out,
    file = output_name
)