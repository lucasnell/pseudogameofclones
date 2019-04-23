

library(clonewars)
library(readxl)


source(".Rprofile")

comp_df <- read_excel(paste0("~/Dropbox/Aphid Project 2017/Lauren_competition/",
                             "competition_data_entry.xlsx")) %>%
    rename_all(~ tolower(.)) %>%
    filter(!is.na(green_line)) %>%
    mutate(green_line = ifelse(green_line == "Clover-2017-6 (reg+",
                               "Clover-2017-6 (reg+)", green_line),
           red_line = ifelse(red_line == "WI-L4(Ham-)", "WI-L4 (Ham-)", red_line)) %>%
    mutate_at(vars(green_line, red_line), ~ gsub(" \\(Ham\\-\\)", "Ø",
                                                 gsub(" \\(reg\\+\\)| \\(Ham\\+\\)", "",
                                                      .x))) %>%
    mutate(total_red = plant_red + dispersed_red,
           total_green = plant_green + dispersed_green) %>%
    mutate(date = as.Date(paste(year, month, day, sep = "-"))) %>%
    select(-year, -month, -day, -observer, -comments,
           -starts_with("dispersed_"), -starts_with("plant_")) %>%
    select(date, everything()) %>%
    mutate(combo = paste0(green_line, red_line) %>%
               factor() %>% as.integer() %>% factor()) %>%
    # gather("color", "line", green_line, red_line) %>%
    # gather("color2", "total", total_green, total_red) %>%
    # mutate(color = gsub("_line", "", color),
    #        color2 = gsub("total_", "", color2)) %>%
    # filter(color == color2) %>%
    # select(-color2) %>%
    group_by(treatment, combo, rep) %>%
    mutate(date = as.integer(date - min(date))) %>%
    ungroup() %>%
    rename(green = green_line, red = red_line) %>%
    mutate(diff = total_green - total_red)


comp_plot <- function(.temp) {
    comp_df %>%
        filter(grepl(sprintf("^%s", paste(.temp)), treatment)) %>%
        split(interaction(.$combo, .$rep)) %>%
        map_dfr(function(.x) {
            .mean <- mean(c(.x$total_red, .x$total_green))
            .sd <- sd(c(.x$total_red, .x$total_green))
            mutate(.x,
                   total_red = (total_red - .mean) / .sd,
                   total_green = (total_green - .mean) / .sd,
                   diff_z = total_green - total_red)
        }) %>%
        ggplot(aes(date)) +
        # geom_line(aes(y = total_green, linetype = factor(rep))) +
        # geom_line(aes(y = total_red, linetype = factor(rep)), color = "red") +
        geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
        geom_line(aes(y = diff_z, group = factor(rep)), color = "gray50") +
        geom_point(aes(y = diff_z, fill = diff_z), shape = 21, color = "gray40", size = 2) +
        facet_grid(red ~ green) + #, labeller = label_both) +
        # scale_color_gradient2(low = "#d95f02", high = "#1b9e77", mid = "gray80",
        scale_fill_gradient2(low = "firebrick", high = "chartreuse", mid = "gray60",
                             midpoint = 0, guide = FALSE) +
        scale_y_continuous(expression("Scaled green" - "red"), breaks = c(-2, 0, 2)) +
        scale_x_continuous("Date", breaks = seq(0, 30, 10)) +
        theme(strip.text = element_text(size = 10))
}

comp_plot(20) %>%
    ggsave(filename = "~/Desktop/comp_20.pdf", height = 6, width = 6)
comp_plot(27) %>%
    ggsave(filename = "~/Desktop/comp_27.pdf", height = 6, width = 6)







library(clonewars)
source(".Rprofile")

{
    ggplot(data = tibble(x = 0:10), aes(x)) +
        stat_function(fun = sin, size = 1) +
        stat_function(fun = function(x) sin(x - (pi/1.5)), color = "red", size = 1) +
        coord_cartesian(ylim = c(-1,1)) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
} %>%
    ggsave(filename = "~/Desktop/wa.pdf", width = 4.75, height = 2.25)

{
    ggplot(data = tibble(x = 0:20), aes(x)) +
        stat_function(fun = sin, size = 1) +
        stat_function(fun = function(x) sin(x - (pi/2)), color = "red", size = 1) +
        coord_cartesian(ylim = c(-1,2)) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
    } %>%
    ggsave(filename = "~/Desktop/wa_beetles.pdf", width = 4.75, height = 2.25)

{
    ggplot(data = tibble(x = 0:20), aes(x)) +
        stat_function(fun = function(x) {
            z <- x
            z[x <= 1.5 * pi] <- (1 + sin(x[x <= 1.5 * pi])) / 2
            z[x > 1.5 * pi] <- 1.5 * (1 + sin(x[x > 1.5 * pi])) / 2
            z[x > 2.5 * pi] <- 0.5 + (1 + sin(x[x > 2.5 * pi])) / 2
            return(z)
        },
        size = 1) +
        stat_function(fun = function(x) {
            ifelse(x <= 2 * pi, 1, 0.5) * (1 + sin(x - (pi/2))) / 2
        },
        color = "red", size = 1) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) +
        NULL
    } %>%
    ggsave(filename = "~/Desktop/wa_beetles2.pdf", width = 4.75, height = 2.25)



#

















