#!/opt/conda/bin/Rscript

library(mSigAct)

args <- commandArgs(trailingOnly = TRUE)

spectra_file <- args[1]
catalog_file <- args[2]
target_sig   <- args[3]
cpus         <- as.integer(args[4])
output_name  <- args[5]

df_equi <- read.csv(
    spectra_file,
    row.names = 1
)

df_p <- read.table(
    catalog_file,
    row.names = 1,
    header = TRUE
    )

rownames(df_equi) <- gsub(
    "[^A-Za-z]",
    "",
    rownames(df_equi)
)

rownames(df_p) <- gsub(
    "[^A-Za-z]",
    "",
    rownames(df_p)
)

common_rows <- intersect(
    rownames(df_equi),
    rownames(df_p)
)

if (length(common_rows) != 96) {
    stop(
        "Row names do not match between df_equi and df_p."
    )
}

df_equi <- df_equi[
    common_rows,
    ,
    drop = FALSE
]

df_p <- df_p[
    common_rows,
    ,
    drop = FALSE
]

spectra_mat <- as.matrix(df_equi)

mode(spectra_mat) <- "numeric"

sigs_mat <- as.matrix(df_p)

mode(sigs_mat) <- "numeric"

if (!target_sig %in% colnames(sigs_mat)) {
    stop(
        paste(
            "Target signature",
            target_sig,
            "is not present in signature columns:",
            paste(
                colnames(sigs_mat),
                collapse = ", "
            )
        )
    )
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