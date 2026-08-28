#!/usr/bin/env Rscript

library(dplyr)
library(purrr)
library(tibble)

rdata_file <- commandArgs(trailingOnly = TRUE)[1]
duplicate_id <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
iteration <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
n_injected <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
summary_file <- commandArgs(trailingOnly = TRUE)[5]
exp_with_file <- commandArgs(trailingOnly = TRUE)[6]
exp_without_file <- commandArgs(trailingOnly = TRUE)[7]

load(rdata_file)

sample_names <- names(sig.presence.test.out)

df_summary <- map_dfr(sample_names, function(s_name) {

    
res <- sig.presence.test.out[[s_name]]

tibble(
    duplicate_id = duplicate_id,
    iteration = iteration,
    sample = s_name,
    loglh_with = res$loglh.with,
    loglh_without = res$loglh.without,
    statistic = res$statistic,
    chisq_p = res$chisq.p,
    n_reconstructed_with = sum(res$exp.with),
    n_reconstructed_without = sum(res$exp.without),
    n_injected = n_injected
)
    

})

df_exp_with <- map_dfr(sample_names, function(s_name) {

    
vec <- sig.presence.test.out[[s_name]]$exp.with

as_tibble(as.list(vec)) %>%
    add_column(
        duplicate_id = duplicate_id,
        iteration = iteration,
        sample = s_name,
        n_injected = n_injected,
        .before = 1
    )
    

})

df_exp_without <- map_dfr(sample_names, function(s_name) {

    
vec <- sig.presence.test.out[[s_name]]$exp.without

as_tibble(as.list(vec)) %>%
    add_column(
        duplicate_id = duplicate_id,
        iteration = iteration,
        sample = s_name,
        n_injected = n_injected,
        .before = 1
    )
    

})

print(df_summary)
print(df_exp_with)
print(df_exp_without)

write.csv(df_summary, summary_file, row.names = FALSE)
write.csv(df_exp_with, exp_with_file, row.names = FALSE)
write.csv(df_exp_without, exp_without_file, row.names = FALSE)
