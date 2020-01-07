

#'
#' These assays measured where alates went and how many died.
#' The number of alates that each rep started with weren't input, so this needs to
#' be fixed before starting.
#'

# library(clonewars)
library(tidyverse)
library(readxl)


source(".Rprofile")

fn <- "~/Dropbox/Aphid Project 2017/alate_dispersal/alate_dispersal_data_entry19July19.xlsx"

disp_df <- read_excel(fn) %>%
    filter(!grepl("side", position, ignore.case = TRUE),
           !is.na(starting_instar)) %>%
    rename(instar0 = starting_instar, pot0 = starting_pot, n_juvs = n_juveniles) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%B-%d")) %>%
    select(date, everything(), -year, -month, -day, -time, -notes) %>%
    mutate_at(vars(aphid_line, rep), factor) %>%
    mutate_at(vars(instar0, pot0, position, n_alates, n_juvs), as.integer) %>%
    mutate(dist = abs(position - pot0)) %>%
    group_by(aphid_line, rep) %>%
    mutate(date = as.integer(date - min(date))) %>%
    ungroup()


disp_df %>%
    mutate(rep = paste(rep),
           rep = case_when(
               date > 20 & rep == "1" ~ "3",
               date > 20 & rep == "2" ~ "4",
               TRUE ~ rep
               ),
           rep = factor(rep)) %>%
    group_by(aphid_line, rep) %>%
    mutate(date = as.integer(date - min(date)),
           dist = factor(dist)) %>%
    ungroup() %>%
    ggplot(aes(date, n_alates)) +
    # ggplot(aes(date, log(n_juvs))) +
    geom_line(aes(color = dist)) +
    geom_point(aes(color = dist)) +
    facet_grid(rep ~ aphid_line) +
    # theme_classic() +
    NULL
