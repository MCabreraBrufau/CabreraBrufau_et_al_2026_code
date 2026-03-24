
#Example of "custom weighting scheme" using emmeans to layer the weighting: emmean is proportional to vegpresence within each season, but each season has the same weight in final emmean. 



library(lme4)
library(emmeans)
library(dplyr)

set.seed(123)

# Base grid
dat <- expand.grid(
  status = c("A","B"),
  season = c("spring","summer","fall","winter"),
  site = 1:20
)

# Assign unequal vegpresence probabilities per status × season
veg_probs <- expand.grid(
  status = c("A","B"),
  season = c("spring","summer","fall","winter")
) %>%
  mutate(p_T = runif(n(), 0.3, 0.8))

dat <- dat %>%
  left_join(veg_probs, by = c("status","season")) %>%
  mutate(
    vegpresence = ifelse(runif(n()) < p_T, "T", "F")
  ) %>%
  dplyr::select(-p_T)

# Simulate response
dat$flux <- rnorm(nrow(dat),
                  mean = ifelse(dat$status=="A", 5, 7) +
                    ifelse(dat$vegpresence=="T", 2, -1) +
                    as.numeric(factor(dat$season)),
                  sd = 1)

# Fit mixed model
model <- lmer(flux ~ status * season * vegpresence + (1|site), data = dat)

emm_season <- emmeans(
  model,
  ~ status | season,
  weights = "proportional"
)
emm_season
annual_emm <- emmeans(
  emm_season,
  ~ status,
  weights = "equal"
)

annual_emm

equal_emm<- emmeans(
  model,
  ~ status,
  weights = "equal"
)

prop_emm<-  emmeans(
  model,
  ~ status,
  weights = "proportional"
)

annual_emm
equal_emm
prop_emm

dat %>%
  group_by(status, season, vegpresence) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(status, season) %>%
  mutate(prop = n/sum(n))

pairs(annual_emm)
pairs(equal_emm)
pairs(prop_emm)



emmeans(model, ~ status | season * vegpresence, weights = "show.levels")

