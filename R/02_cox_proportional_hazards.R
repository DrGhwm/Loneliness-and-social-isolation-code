# 2. Cox proportional-hazards regression
library(survival)
library(broom)
library(dplyr)

# Required local object: metabolic
exposures <- c("standardized_metabolic_score_loneliness",
               "standardized_metabolic_score_isolation",
               "loneliness_status","isolation_status",
               "loneliness_index","isolation_index")
covars <- c("age","sex","edu3","Ethnic_group","Towns","units",
            "smoke","BMI","Diabetes","hyper","CVD")

outcomes <- list(
  dementia=c("dementia","dementia_durations"),
  vascular_dementia=c("dementia_vascular","dementia_vascular_durations"),
  alzheimer=c("dementia_alzheimer","dementia_alzheimer_durations"),
  stroke=c("stroke","stroke_durations"),
  parkinson=c("pd","durations"),
  depression=c("depression","depression_durations"),
  anxiety=c("anxiety","anxiety_durations")
)

res <- list()
for (expo in exposures) {
  for (name in names(outcomes)) {
    y <- outcomes[[name]][1]
    time <- outcomes[[name]][2]
    adj <- if (grepl("dementia",name)) c(covars,"apoe_carrier") else covars

    fit <- coxph(as.formula(paste0(
      "Surv(",time,",",y,") ~ ",paste(c(expo,adj),collapse="+")
    )), data=metabolic)

    z <- tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==expo)
    res[[length(res)+1]] <- data.frame(
      Exposure=expo, Outcome=name, HR=z$estimate,
      LCL=z$conf.low, UCL=z$conf.high, P=z$p.value
    )
  }
}
cox_results <- bind_rows(res) %>%
  group_by(Exposure) %>% mutate(FDR=p.adjust(P,"BH")) %>% ungroup()
write.csv(cox_results,"cox_results.csv",row.names=FALSE)
