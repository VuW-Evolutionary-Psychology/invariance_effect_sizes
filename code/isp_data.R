library(tidyverse)
library(lavaan)

### https://osf.io/96stw/
ispdata <- read_csv("isp.csv")

ispdata$SHS4 <- 8-ispdata$SHS4

selected_countries <- 
ispdata %>%
  group_by(.$country) %>%
  count() %>%
  filter(., n >= 100) %>%
  .[[1]]

ispdata <- filter(ispdata, country %in% selected_countries)

multi_eff <- 
multi_group_eff(ispdata, "country", paste0("SHS", 1:4))

multi_average <- multi_eff$average
multi_average$direction <- if_else(multi_average$item == "SHS4", "neg", "pos")


multi_average <- rename(multi_average, "sd_av" =  "sd")
sd_sum <- Rmisc::summarySE(multi_average, measurevar = "sd_av", groupvars = c("direction", "eff"))
av_sum <- Rmisc::summarySE(multi_average, measurevar = "av", groupvars = c("direction", "eff"))

sd_sum %>%
  filter(., eff %in% c("UDI2", "WUDI", "dmacs")) %>%
  ggplot() +
  aes(y = sd_av, x = direction, fill = eff) +
  geom_bar(stat = "identity", position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, vjust = .30))


av_sum %>%
  filter(., eff %in% c("UDI2", "WUDI", "dmacs")) %>%
  ggplot() +
  aes(y = av, x = direction, fill = eff) +
  geom_bar(stat = "identity", position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, vjust = .30))
