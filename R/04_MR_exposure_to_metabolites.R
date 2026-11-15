# 4. Two-sample MR: loneliness/social isolation -> metabolites
library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)

# Set OPENGWAS_JWT locally; never paste the token into this script.
if (Sys.getenv("OPENGWAS_JWT")=="") stop("OpenGWAS JWT is not configured.")

lon_map <- fread("data/Loneliness_Index_MR_metabolite_mapping.csv")
si_map <- fread("data/Social_isolation_Index_MR_metabolite_mapping.csv")

# Genetic instruments
lon_exp <- extract_instruments("ukb-b-8476",p1=5e-8) %>% clump_data()
si_file <- "data/raw/Social_isolation_genome_wide_sig_summary_statistics_20230925.txt"
si_exp <- fread(si_file) %>% as.data.frame() %>% clump_data()

run_mr <- function(exp_dat, mapping, exposure_name) {
  res <- list()
  for (i in seq_len(nrow(mapping))) {
    out_dat <- extract_outcome_data(
      snps=exp_dat$SNP, outcomes=mapping$GWAS_ID[i])
    if (is.null(out_dat) || nrow(out_dat)==0) next

    dat <- harmonise_data(exp_dat,out_dat) %>% filter(mr_keep)
    if (nrow(dat)<3) next

    fit <- mr(dat,method_list=c(
      "mr_ivw","mr_weighted_median","mr_egger_regression"))
    fit$Exposure <- exposure_name
    fit$Metabolite <- mapping$Metabolite[i]
    fit$GWAS_ID <- mapping$GWAS_ID[i]
    res[[i]] <- fit
  }
  bind_rows(res)
}

lon_res <- run_mr(lon_exp,lon_map,"Loneliness")
si_res <- run_mr(si_exp,si_map,"Social isolation")
mr_results <- bind_rows(lon_res,si_res)
write.csv(mr_results,"exposure_to_metabolite_MR.csv",row.names=FALSE)

ivw <- filter(mr_results,method=="Inverse variance weighted") %>%
  group_by(Exposure) %>% mutate(FDR=p.adjust(pval,"BH")) %>% ungroup()
write.csv(ivw,"exposure_to_metabolite_MR_IVW_FDR.csv",row.names=FALSE)
