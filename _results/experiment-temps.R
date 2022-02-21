
library(tidyverse)
library(clonewars)

source(".Rprofile")


list.files("~/Box Sync/eco-evo_experiments/temps_humids", full.names = TRUE) %>%
    .[-1] %>% # <-- this one's redundant
    imap_dfr(~ read_csv(.x) %>% set_names(c("time", "temp", "humid")) %>%
                mutate(f = .y)) %>%
    mutate(f = factor(f)) %>%
    group_by(f) %>%
    summarize(min = min(time), max = max(time))


    ggplot(aes(time, temp, color = f)) +
    geom_line()
