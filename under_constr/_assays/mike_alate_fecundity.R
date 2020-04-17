


library(tidyverse)
library(clonewars)
library(viridisLite)

source(".Rprofile")

files <- sprintf("~/Box Sync/2020/aphids/alate_fecundity/%s.csv",
                 c("Al_Fecundity_ClonalLines_v1_29Mar2019",
                   "Fecundity_v2_22Apr19"))

# Apterous data:
apterous_fecund <- read_csv(files[1], col_types = cols()) %>%
    rename_all(tolower) %>%
    filter(treatment == "Control") %>%
    rename(rep = replicate, off = ap_daily) %>%
    select(line, rep, day, off) %>%
    filter(!is.na(off))

# Alate data:
alate_fecund <- read_csv(files[2], col_types = cols()) %>%
    rename_all(tolower) %>%
    rename(off = al_daily) %>%
    select(line, rep, day, off) %>%
    filter(!is.na(off))


fecund_df <- bind_rows(apterous_fecund %>% mutate(type = "apterous"),
                       alate_fecund %>% mutate(type = "alate")) %>%
    # There can't be offspring on day 0 based on the assay methods:
    filter(day > 0) %>%
    mutate(line = ifelse(line == "WI-L4_Ham+", "WI-L4 Ham+", line),
           line = ifelse(line == "WI-L4 Ham+", "WI-L4", line),
           type = factor(type, levels = c("apterous", "alate"),
                         labels = c("non-winged", "winged")))


boot_fecund <- fecund_df %>%
    split(interaction(.$line, .$type, .$day, drop = TRUE)) %>%
    map_dfr(function(.x) {
        X <- map_dbl(1:1000, function(.) {
            inds <- sample.int(nrow(.x), nrow(.x), replace = TRUE)
            mean(.x[["off"]][inds])
        })
        tibble(line = .x$line[1],
               type = .x$type[1],
               day = .x$day[1],
               low = as.numeric(quantile(X, 0.025)),
               avg = mean(.x$off),
               high = as.numeric(quantile(X, 0.975)))
    })



fecund_df %>%
    mutate(type2 = interaction(line, type)) %>%
    # filter(type == "alate") %>%
    ggplot(aes(day, off, color = line)) +
    # geom_point() +
    geom_line(data = tibble(day = 1:20,
                            off = populations$repro$low[1:20] +
                                    populations$repro$high[1:20] / 2),
              alpha = 0.5, size = 1, color = "gray70") +
    geom_line(aes(group = rep, color = type2), alpha = 0.5, size = 0.5) +
    geom_line(data = boot_fecund %>% mutate(type2 = interaction(line, type)),
              aes(y = avg, color = type2), size = 1) +
    facet_grid(line ~ type) +
    scale_color_manual(values = c(viridis(1, begin = 0.9),
                                  viridis(1, begin = 0.7),
                                  magma(1, begin = 0.2),
                                  magma(1, begin = 0.5)),
                       guide = FALSE) +
    ylab("Daily offspring") +
    xlab("Age (days)") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 12))


L <- leslie_matrix(dev_times$instar_days$lowT, populations$surv_juv$high,
                   populations$surv_adult$high, populations$repro$high)

# Average of the data from Meisner et al. 2014
avg_fecund <- populations$repro$low + populations$repro$high / 2
avg_fecund <- avg_fecund[1:nrow(L)]

fecunds <- boot_fecund %>%
    filter(type == "winged") %>%
    mutate(prop = map2_dbl(day, avg, ~ .y / avg_fecund[.x])) %>%
    filter(!is.nan(prop)) %>%
    split(.$line) %>%
    map(~ lm(prop ~ I(gtools::inv.logit(day)), .x)) %>%
    map_dfr(~ tibble(day = 1:length(avg_fecund),
                     fecund = avg_fecund *
                         predict(.x, newdata = tibble(day = 1:length(avg_fecund))))) %>%
    mutate(line = c(rep("UT3", length(avg_fecund)), rep("WI-L4", length(avg_fecund))))


fecunds %>%
    ggplot(aes(day, fecund)) +
    theme_classic() +
    geom_line(data = fecund_df %>% filter(type == "winged"),
              aes(y = off, color = line, group = rep),
              alpha = 0.5, size = 0.5) +
    geom_line(aes(color = line), size = 1) +
    ylab("Daily offspring") +
    xlab("Age (days)") +
    theme_classic() +
    facet_grid(line ~ .) +
    scale_color_manual(values = c(viridis(1, begin = 0.9),
                                  magma(1, begin = 0.2)))




#'
#' Object to use in the `leslie_matrix` function:
#'
fecunds <- fecunds %>%
    split(.$line) %>%
    map(~ .x[["fecund"]])
fecunds



