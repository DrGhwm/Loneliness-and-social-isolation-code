# 3. Linear regression for brain-related phenotypes
library(broom)
library(dplyr)

# Required local object: metabolic (phenotypes and covariates prepared)
exposures <- c("standardized_metabolic_score_loneliness",
               "standardized_metabolic_score_isolation",
               "loneliness_status","isolation_status",
               "loneliness_index","isolation_index")
covars <- c("age","sex","edu3","Ethnic_group","Towns","units",
            "BMI","smoke","Diabetes","hyper","CVD")

movement <- c("Hand.grip.strength..left....Instance.2",
              "Hand.grip.strength..right....Instance.2",
              "Usual.walking.pace...Instance.2",
              "Falls.in.the.last.year...Instance.2")
cognition <- c("Maximum.digits.remembered.correctly...Instance.2",
               "Fluid.intelligence.score...Instance.2",
               "Prospective.memory.result...Instance.2",
               "Number.of.puzzles.correct...Instance.2",
               "Number.of.symbol.digit.matches.made.correctly...Instance.2",
               "Number.of.word.pairs.correctly.associated...Instance.2",
               "Mean.time.to.correctly.identify.matches...Instance.2",
               "Number.of.incorrect.matches.in.round...Instance.2...Array.1",
               "Duration.to.complete.numeric.path..trail..1....Instance.2",
               "Duration.to.complete.alphanumeric.path..trail..2....Instance.2")
mental_health <- c("GMS_2","PHQ4_2","PHQ9","CIDI","GAD7",
                   "LS_2","Sub_w","upe","TE")

res <- list()
for (expo in exposures) {
  for (y in c(movement,cognition,mental_health)) {
    fit <- lm(reformulate(c(expo,covars),response=y),data=metabolic)
    z <- tidy(fit) %>% filter(term==expo)
    res[[length(res)+1]] <- data.frame(
      Exposure=expo, Outcome=y, Beta=z$estimate,
      SE=z$std.error, P=z$p.value, N=nobs(fit)
    )
  }
}
linear_results <- bind_rows(res) %>%
  group_by(Exposure) %>% mutate(FDR=p.adjust(P,"BH")) %>% ungroup()
write.csv(linear_results,"linear_brain_phenotypes.csv",row.names=FALSE)

# Brain-structure traits can be analysed with the same lm() model
# after preparing the brain-structure data locally.
