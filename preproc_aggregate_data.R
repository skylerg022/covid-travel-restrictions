## The purpose of this file is to create coarser data to model.
## daily reported numbers will be summarized as monthly means.

library(tidyverse)
library(lubridate)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

covid <- read_csv('data/owid-covid-data.csv') %>%
  filter(!is.na(continent)) %>% # Remove any observations that are not country-specific
  filter(date <= '2021-10-31') %>% # Remove observations after Oct 2021
  inner_join(read_csv('data/international-travel-covid.csv') %>%
              rename(location = Entity,
                     date = Day),
            by = c('location', 'date')) %>%
  rename(int_travel_controls = international_travel_controls) %>%
  select(location, continent, date, new_cases, new_deaths,
         people_fully_vaccinated, population,
         population_density, median_age, aged_65_older, gdp_per_capita,
         int_travel_controls)

mode <- function(x) {
  counts <- table(x)
  mode_idx <- which.max(counts)
  val <- as.numeric(names(counts)[mode_idx])
  
  return(val)
}

covid2 <- covid %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(location, continent, year, month,
           population, population_density,
           median_age, aged_65_older,
           gdp_per_capita) %>%
  summarize(avg_new_cases = mean(new_cases, na.rm = TRUE),
            avg_new_deaths = mean(new_deaths, na.rm = TRUE),
            mid_people_fully_vaccinated = median(people_fully_vaccinated, na.rm = TRUE),
            int_travel_controls = mode(int_travel_controls) + 1) %>%
  mutate(vaccination_info = ifelse(is.na(mid_people_fully_vaccinated), 0, 1)) %>%
  replace_na(list(mid_people_fully_vaccinated = 0,
                  avg_new_cases = 0,
                  avg_new_deaths = 0)) %>%
  # Taiwan has missing values; give ballpark estimates
  mutate(population_density = case_when(location == 'Taiwan' ~ 652, # From population table I pulled
                                        TRUE                 ~ population_density),
         # GDP PPP of 2020 adjusted for inflation to get approx. 2011 constant international dollars
         gdp_per_capita = case_when(location == 'Taiwan' ~ 48500,
                                    TRUE                 ~ gdp_per_capita)) %>%
  group_by(location, continent) %>%
  mutate(time = row_number()) %>%
  ungroup() %>%
  # Calculate percent population fully vaccinated
  mutate(perc_fully_vaccinated = mid_people_fully_vaccinated / population) %>%
  # Reorganize variables
  select(location, continent, time, year, month,
         population_density, median_age, gdp_per_capita, 
         avg_new_cases, avg_new_deaths,
         vaccination_info, perc_fully_vaccinated,
         int_travel_controls) %>%
  # Remove remaining 300+ observations without pop dens, median age, or gdp
  filter(!is.na(population_density), 
         !is.na(median_age), 
         !is.na(gdp_per_capita))


# Write to csv file
write_csv(covid2, 'data/international_controls.csv')
