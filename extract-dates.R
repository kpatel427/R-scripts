library(tidyverse)
library(lubridate)

#' @example
#' extract_days(x = "2 weeks and 1 day")
extract_days <- function(x) {
  tolower(x) %>% 
    str_remove("and ") %>% 
    duration() %>% 
    seconds_to_period() %>% 
    day()
}





