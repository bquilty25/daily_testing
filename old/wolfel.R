# wolfel data

library(tidyverse)

list(`1` = c(3, 4, 5, 6, 7, 8, 8),  # still infectious
     `0` = c(5, 7, 8, 8, 9, 10, 11, 11, 13)) %>% # no longer infectious
  map_df(~data.frame(day = .x), .id = "y") %>%
  mutate(y = 1 - as.integer(y)) -> wolfel

# have they yet become non-infectious?
wolfel_glm <- glm(data=wolfel, formula = y ~ day, family = "binomial")

wolfel_pred <- data.frame(day = seq(0,30,by=0.01)) %>%
  mutate(y = predict(wolfel_glm, ., type = 'response')) %>%
  mutate(y = ifelse(day == 0, 0, y))

