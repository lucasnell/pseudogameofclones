library(dplyr)
library(readr)

fast <- read_csv("data-raw/leslie_fast.csv", col_names = FALSE,
                 col_types = cols(.default = col_double())) %>%
    setNames(NULL) %>%
    as.matrix()
slow <- read_csv("data-raw/leslie_slow.csv", col_names = FALSE,
                 col_types = cols(.default = col_double())) %>%
    setNames(NULL) %>%
    as.matrix()

leslie <- list(fast = fast, slow = slow)

devtools::use_data(leslie)
