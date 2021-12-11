# Initial EDA

library(tidyverse)
library(cowplot)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Set ggplot theme
old_theme <- theme_set(theme_bw())
# Picture saving settings
pic_width <- 7
pic_height <- 4
pic_unit <- 'in'

# Read in data
covid <- read_csv('data/international_controls.csv')

# Number of countries in analysis by continent
covid %>%
  group_by(continent) %>%
  summarize(n_countries = length(unique(location)),
            n = n())

# Visualize control levels by continent/country
p <- covid %>%
  mutate(Taiwan = ifelse(location %in% c('Taiwan', 'United States'), 'Yes', 'No')) %>%
  ggplot(aes(x = time, y = jitter(int_travel_controls), group = location,
             col = Taiwan)) +
  geom_line(alpha = 0.5, size = 1) +
  facet_wrap(~ continent, ncol = 3) +
  scale_color_manual(values = c('grey', 'red')) +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'International Travel Control Level') +
  theme(legend.position = 'none')

ggsave('plots/eda_control_levels.png',
       plot = p,
       device = 'png',
       dpi = 300,
       width = pic_width, 
       height = pic_height,
       units = pic_unit)


# Taiwan analysis
control_tw <- covid %>%
  filter(location %in% c('United States', 'Taiwan'))
p1 <- control_tw %>%
  ggplot(aes(time, log(avg_new_cases), col = location)) +
  geom_line() +
  scale_color_manual(values = c('red', 'royalblue')) +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'Log(Average\nNew Cases)',
       col = 'Country')
p2 <- control_tw %>%
  ggplot(aes(time, perc_fully_vaccinated, col = location)) +
  geom_line() +
  scale_color_manual(values = c('red', 'royalblue')) +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'Proportion\nFully Vaccinated',
       col = 'Country')
p3 <- control_tw %>%
  ggplot(aes(time, int_travel_controls, col = location)) +
  geom_line() +
  scale_color_manual(values = c('red', 'royalblue')) +
  labs(x = 'Month (1 = Jan 2020)',
       y = 'Int. Travel\nControl Status',
       col = 'Country') +
  scale_y_continuous(breaks= c(1, 2, 3, 4, 5))

# Save plot
plot_grid(p1, p2, p3, nrow = 3) %>%
  ggsave('plots/eda_taiwan.png',
         plot = .,
         device = 'png',
         dpi = 300,
         width = pic_width, 
         height = pic_height,
         units = pic_unit)
