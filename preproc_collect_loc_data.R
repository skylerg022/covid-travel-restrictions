library(tidyverse)
library(rvest)

# set working directory if necessary
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}



# Read in continent data --------------------------------------------------

url <-  "https://www.newworldencyclopedia.org/entry/List_of_countries_by_continent"

tables <- url %>% 
  read_html() %>% 
  html_elements('table') %>%
  html_table() %>%
  .[1:6] %>%
  map(function(x) rbind(x[,1:2], x[,3:4]))

continents <- c('Africa', 'Asia', 'Europe',
                'North America', 'South America',
                'Oceania')

continent <- tibble(continent = continents,
                    table = tables) %>%
  unnest(table) %>%
  rename(country = Country) %>%
  select(continent, country) %>%
  # Clean names
  mutate(country = case_when(str_detect(country, 'Taiwan') ~ 'Taiwan',
                             str_detect(country, 'North Korea') ~ 'North Korea',
                             str_detect(country, 'South Korea') ~ 'South Korea',
                             str_detect(country, 'Swaziland') ~ 'Eswatini',
                             str_detect(country, 'East Timor') ~ 'East Timor',
                             str_detect(country, 'Palestinian territories') ~ 'Palestine',
                             str_detect(country, 'Former Yugoslav Republic of Macedonia') ~ 'North Macedonia',
                             TRUE                          ~ country),
         country = country %>%
           str_remove_all(' ?\\(.+\\)') %>% # Remove parentheses
           str_remove_all(' ?\\[.+\\]') %>% # Remove links
           str_replace_all(' {2,}', ' ') %>% # Transform multiple spaces to 1
           str_remove(',.*') %>% # Commas adn whatever comes after them
           str_remove('\\.*') %>% # Periods
           str_remove('^ +') %>% # Leading spaces
           str_remove(' +$') %>% # Ending spaces
           stringi::stri_trans_general("Latin-ASCII") %>%
           tolower()) %>%
  filter(!is.na(country), country != '') %>%
  arrange(continent, country)


# Read in population data --------------------------------------------------

url <-  "https://en.wikipedia.org/wiki/List_of_countries_and_dependencies_by_population_density"

population <- url %>% 
  read_html() %>% 
  html_elements('table') %>%
  html_table() %>%
  # Population data is all contained in first table
  .[[1]] %>%
  .[-1,] %>%
  # Clean data
  janitor::clean_names() %>%
  select(-c(area_2, density_2)) %>%
  rename(country = country_or_territory) %>%
  # select(-c(area, density)) %>% # Keep area/density sqmi measures
  # rename(area = area_2, density = density_2,
  #        country = country_or_territory) %>%
  # Convert measures to numeric columns
  mutate(across(c(population, area, density),
                function(x) x %>%
                  str_remove_all(',') %>%
                  as.numeric())) %>%
  mutate(country = case_when(str_detect(country, 'U\\.S\\. Virgin Islands') ~ 'United States Virgin Islands',
                             str_detect(country, 'DR') ~ 'Congo',
                             str_detect(country, 'Ivory Coast') ~ 'Cote d\'Ivoire',
                             TRUE                            ~ country)) %>%
  # Clean country/territory names
  mutate(country = country %>%
           str_remove(' ?\\*.*') %>%
           str_remove(' ?\\(.*\\)') %>%
           # str_remove_all('\\.') %>% # Periods
           stringi::stri_trans_general("Latin-ASCII") %>%
           str_remove('^ +') %>% # Leading spaces
           str_remove('[^a-zA-Z]+$') %>% # Ending stuff that isn't part of the name
           str_replace_all(' {2,}', ' ') %>% # Transform multiple spaces to 1
           str_remove(',.*') %>% # Comma and everything that comes after
           tolower())

# Check performance -------------------------------------------------------

# Those remaining appear to be territories
# anti_join(continent, population, by = 'country')

# South Sudan doesn't have a match; all others have small population
# anti_join(population, continent, by = 'country')

countries <- continent %>%
  inner_join(population, by = 'country') %>%
  select(country, continent, population:sources_dates)

write_csv(countries, 'data/country_info.csv')
