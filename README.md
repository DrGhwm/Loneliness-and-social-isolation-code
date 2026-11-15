# Loneliness and Social Isolation: Analysis Code

Four concise R scripts are provided:

1. `R/01_metabolite_selection_elastic_net.R` — metabolite screening and elastic-net scores.
2. `R/02_cox_proportional_hazards.R` — Cox regression for incident brain-health outcomes.
3. `R/03_linear_regression_brain_phenotypes.R` — linear regression for movement, cognition and mental-health phenotypes.
4. `R/04_MR_exposure_to_metabolites.R` — two-sample MR from loneliness/social isolation to candidate metabolites.

The scripts expect locally prepared analysis objects (e.g. `metabolic`, `metabolic_wave1` and `metabolism0`). Individual-level UK Biobank data are not distributed.

GWAS mapping tables are included in `data/`. Obtain the Bristol social-isolation summary-statistics file from its original repository and place it in `data/raw/`. Restricted data and API tokens must not be uploaded.

Configure `OPENGWAS_JWT` in your local R environment before running the OpenGWAS queries. Do not hard-code credentials in the scripts.
