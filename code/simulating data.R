library(lavaan)
library(simsem)
library(tidyverse)
popNoninvariance1 <- "
f1 =~ 0.70*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5 + 0.7*y6
f1 ~~ 1*f1
y1 ~~ 0.51*y1
y2 ~~ 0.51*y2
y3 ~~ 0.51*y3
y4 ~~ 0.51*y4
y5 ~~ 0.51*y5
y6 ~~ 0.51*y6
f1 ~ 1
y1 ~ -0.2*1
y2 ~ 0.1*1
y3 ~ 0*1
y4 ~ 0.1*1
y5 ~ 0.2*1
y6 ~ 0.3*1
"

popDataNoninvariance1 <- simulateData(popNoninvariance1,
                                      sample.nobs = 1000,
                                      seed = 2711,
                                      model.type = "cfa")

out <- 
lapply(seq(-1, 1, .10), function(x){
popNoninvariance2 <- paste0("
f1 =~ ",x,"*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5 + 0.7*y6
f1 ~~ 1*f1
y1 ~~ 0.51*y1
y2 ~~ 0.51*y2
y3 ~~ 0.51*y3
y4 ~~ 0.51*y4
y5 ~~ 0.51*y5
y6 ~~ 0.51*y6
f1 ~ 1
y1 ~ -0.2*1
y2 ~ 0.1*1
y3 ~ 0*1
y4 ~ 0.1*1
y5 ~ 0.2*1
y6 ~ 0.3*1
"
)

popNoninvariance2 <- paste0("
f1 =~ ",.3,"*y1 + 0.2*y2 + -0.7*y3 + 0.7*y4 + 0.7*y5 + 0.7*y6
f1 ~~ 1*f1
y1 ~~ 0.51*y1
y2 ~~ 0.51*y2
y3 ~~ 0.51*y3
y4 ~~ 0.51*y4
y5 ~~ 0.51*y5
y6 ~~ 0.51*y6
f1 ~ 1
y1 ~ -0.2*1
y2 ~ 0.1*1
y3 ~ 0*1
y4 ~ 0.1*1
y5 ~ 0.2*1
y6 ~ 0.3*1
"
)
popDataNoninvariance2 <- simulateData(popNoninvariance2,
                                      sample.nobs = 1000,
                                      seed = 2711,
                                      model.type = "cfa")
popDataNoninvariance <- rbind(popDataNoninvariance1, popDataNoninvariance2)
popDataNoninvariance <- data.frame(popDataNoninvariance, group = rep(c(1, 2, 3, 4 ,5 ), each = 400))



multi_group_eff(popDataNoninvariance, "group", paste0("y", 1:6))
test <- cfa('f1 =~ y1 + y2 + y3 +y4 +y5 +y6', popDataNoninvariance, group = "group", std.lv = T)
dmacs::lavaan_dmacs(test)
out_df <- 
boot_inv_eff(popDataNoninvariance, n = 20, n_sample = 500, items = paste0("y", 1:6), group = "group", seed = 11)

out_df %>%
  ggplot() +
  aes(x = item, y = av, fill = eff) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = l, ymax = h), position = "dodge")
