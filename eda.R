# Initial EDA

library(tidyverse)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

covid <- read_csv('data/international_controls.csv')

# Number of countries in analysis by continent
covid %>%
  distinct(location, .keep_all = TRUE) %>%
  group_by(continent) %>%
  summarize(n = n())

# Number of observations for Taiwan
covid %>%
  filter(location == 'Taiwan') %>%
  nrow()

# Visualize 
covid %>%
  mutate(Taiwan = ifelse(location == 'Taiwan', 'Yes', 'No')) %>%
  ggplot(aes(x = time, y = jitter(int_travel_controls), group = location,
             col = Taiwan)) +
  geom_line(alpha = 0.5, size = 1) +
  facet_wrap(~ continent, ncol = 3) +
  scale_color_brewer(type = 'qual', palette = 3) +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'International Travel Control Level') +
  theme_bw()

# Taiwan response across time
covid %>%
  filter(location == 'Taiwan') %>%
  ggplot(aes(x = time, y = int_travel_controls)) +
  geom_line() +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'International Travel Control Level') +
  theme_bw()

# Visualize
covid %>%
  mutate(taiwan = ifelse(location == 'Taiwan', 'Yes', 'No')) %>%
  ggplot(aes(x = time, y = log(avg_new_cases), group = location,
             col = taiwan)) +
  geom_line() +
  facet_wrap(~ continent, ncol = 3) +
  scale_color_brewer(type = 'qual', palette = 3) +
  theme_bw()
