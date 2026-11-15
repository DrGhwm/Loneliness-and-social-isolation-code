# 1. Metabolite screening and elastic-net score
library(dplyr)
library(glmnet)

# Required local objects: metabolic, metabolism0
covars <- c("age","sex","edu3","Ethnic_group","Towns","units",
            "BMI","smoke","Diabetes","hyper","CVD")
exposures <- c("loneliness_index","isolation_index",
               "loneliness_status","isolation_status")
metabolites <- setdiff(names(metabolism0), "id")

# Linear regression screening
res <- list()
for (expo in exposures) {
  for (met in metabolites) {
    fit <- lm(reformulate(c(expo,covars), response=met), data=metabolic)
    z <- coef(summary(fit))[expo,]
    res[[length(res)+1]] <- data.frame(Exposure=expo, Metabolite=met,
                                      Beta=z[1], SE=z[2], P=z[4])
  }
}
screen <- bind_rows(res) %>%
  group_by(Exposure) %>%
  mutate(FDR=p.adjust(P,"BH"),
         Bonferroni=p.adjust(P,"bonferroni")) %>%
  ungroup()
write.csv(screen,"metabolite_screening.csv",row.names=FALSE)

# Elastic-net scores based on Bonferroni-selected metabolites
lon_met <- filter(screen, Exposure=="loneliness_index",
                  Bonferroni<0.05) %>% pull(Metabolite)
si_met <- filter(screen, Exposure=="isolation_index",
                 Bonferroni<0.05) %>% pull(Metabolite)

set.seed(2023)
fit_lon <- cv.glmnet(as.matrix(metabolic[,lon_met]),
                     metabolic$loneliness_index, alpha=0.5, nfolds=10)
fit_si <- cv.glmnet(as.matrix(metabolic[,si_met]),
                    metabolic$isolation_index, alpha=0.5, nfolds=10)

metabolic$score_loneliness <- as.numeric(
  predict(fit_lon,as.matrix(metabolic[,lon_met]),s="lambda.1se"))
metabolic$score_isolation <- as.numeric(
  predict(fit_si,as.matrix(metabolic[,si_met]),s="lambda.1se"))
metabolic_wave1$score_loneliness <- as.numeric(
  predict(fit_lon,as.matrix(metabolic_wave1[,lon_met]),s="lambda.1se"))
metabolic_wave1$score_isolation <- as.numeric(
  predict(fit_si,as.matrix(metabolic_wave1[,si_met]),s="lambda.1se"))
