
library(tidyverse)

source(".Rprofile")


#'
#' Summary:
#'
#' * Both Clover-2017-2 and WI-L4 show response to crowding
#' * WI-L4 has much higher alate production in both treatments
#'
df1 <- read_csv("~/Desktop/alateproductionassays/AlateProduction_02Jul2019.csv") %>%
    rename_all(tolower)
df1 %>%
    group_by(line, type) %>%
    summarize(prop_sd = sd(prop),
              prop = mean(prop))



green <- c("Clover-2017-2", "Clover-2017-6", "UT3", "WI-2016-593")
#'
#' Summary:
#'
#' * UT3 (green) is super low
#' * Clover-2017-2 (green) is high
#' * R10 (red) is low
#' * WI-L4Ø (red) is pretty high (but not as much as in my plant assays)
#'
df2 <- read_csv("~/Desktop/alateproductionassays/AlateProduction_07Aug2019.csv") %>%
    rename_all(tolower) %>%
    filter(!is.na(al), !is.na(ap)) %>%
    mutate(prop = al / (al + ap),
           line = ifelse(line == "Clover-6", "Clover-2017-6", line),
           line = ifelse(line == "WI-L4-Ham-", "WI-L4", line),
           color = ifelse(line %in% green, "green", "red"))

df2_p <- df2 %>%
    ggplot(aes(line, prop, color = color, shape = treat)) +
    geom_point(size = 2) +
    ylab("Alate offspring proportion") +
    scale_color_manual(values = RColorBrewer::brewer.pal(3,"Dark2")[c(1, 2)], guide = FALSE) +
    scale_shape_manual(NULL, values = c(1, 2)) +
    theme_classic() +
    theme(legend.position = c(0.5, 0.9))

df2_p %>%
    ggsave(filename = "~/Desktop/mike_data.pdf", width = 6, height = 3)

df2 %>%
    ggplot(aes(line, al, color = color, shape = treat)) +
    geom_point() +
    scale_color_manual(values = c("green", "red"), guide = FALSE) +
    scale_shape_manual(NULL, values = c(1, 17)) +
    theme_classic() +
    theme(legend.position = c(0.5, 0.9))

