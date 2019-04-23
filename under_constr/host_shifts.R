

library(clonewars)
library(readxl)

source(".Rprofile")
#
# .i <- function(xlcol) {
#     xlcol <- tolower(xlcol)
#     rcol <- sum(sapply(1:nchar(xlcol),
#                        function (i) {
#                            ind_i <- which(letters == substr(xlcol, i, i))
#                            ind <- length(letters)^(nchar(xlcol) - i) * ind_i
#                            return(ind)
#                        }))
#     return(rcol)
# }
# .i('Z'); .i('AA')
#
# start_date <- as.Date("1899-12-30")
#
# rawdata <- read_excel("~/Desktop/host shift data sheet.xlsx", col_types = "text",
#                       col_names = paste0("V", 1:15))
#
# process_count_cells <- function(.x, .y, .host) {
#     rawdata[.x:(.x+6), .y:(.y+1)] %>%
#         set_names(c("total_red", "total_green")) %>%
#         mutate_all(function(x) ifelse(grepl("plant died", x), NA, x)) %>%
#         mutate_all(as.integer) %>%
#         mutate(host = .host)
# }
# add_cols <- function(.i, .dat) {
#     mutate(.dat[[.i]], combo = combos[.i], date = dates[[.i]],
#            red = lines[[.i]][["red"]], green = lines[[.i]][["green"]])
# }
# # (In strings of red and green lines, red always comes first.)
# process_lines <- function(.x) {
#     strsplit(.x, " and ")[[1]] %>%
#         gsub(pattern = "\\(red\\)|\\(green\\)", replacement = "") %>%
#         gsub(pattern = " - ", replacement = "-") %>%
#         trimws() %>%
#         gsub(pattern = " Ham-", replacement = "Ø") %>%
#         gsub(pattern = " ", replacement = "-") %>%
#         set_names(c("red","green"))
# }
#
#
#
# combos <- map2_chr(rep(seq(1, 78, 11), 2), rep(c(1, .i("I")), each = 8),
#                  ~ unlist(rawdata[.x,.y]))
# lines <- map2_chr(rep(seq(1, 78, 11), 2), rep(c(2, .i("J")), each = 8),
#                    ~ unlist(rawdata[.x,.y])) %>%
#     map(process_lines)
# dates <- map2(rep(seq(4, 81, 11), 2), rep(c(1, .i("I")), each = 8),
#               ~ start_date + as.integer(unlist(rawdata[.x:(.x+6), .y])))
# fava <- map2(rep(seq(4, 81, 11), 2), rep(c(2, .i("J")), each = 8),
#              process_count_cells, .host = "fava") %>%
#              {map_dfr(1:length(.), add_cols, .dat = .)}
# clover <- map2(rep(seq(4, 81, 11), 2), rep(c(4, .i("L")), each = 8),
#                process_count_cells, .host = "clover") %>%
#                {map_dfr(1:length(.), add_cols, .dat = .)}
# alfalfa <- map2(rep(seq(4, 81, 11), 2), rep(c(6, .i("N")), each = 8),
#                 process_count_cells, .host = "alfalfa") %>%
#                 {map_dfr(1:length(.), add_cols, .dat = .)}
#
#
# host_df <- bind_rows(fava, clover, alfalfa) %>%
#     mutate(host = factor(host, levels = c("fava", "alfalfa", "clover")),
#            rep = 1L,
#            green = factor(green),
#            red = factor(red)) %>%
#     filter(!is.na(total_red))
#
# # comp_df
#
# host_df %>%
#     filter(host == "clover", red == "R10", green == "Clover-2017-2")
#
# host_plot <- function(.host) {
#     host_df %>%
#         filter(host == .host) %>%
#         mutate(date = difftime(date, min(date), units = "days") %>% as.integer()) %>%
#         split(interaction(.$combo, .$rep)) %>%
#         map_dfr(function(.x) {
#             .mean <- mean(c(.x$total_red, .x$total_green))
#             .sd <- sd(c(.x$total_red, .x$total_green))
#             if (.sd == 0) .sd <- 1
#             mutate(.x,
#                    total_red = (total_red - .mean) / .sd,
#                    total_green = (total_green - .mean) / .sd,
#                    diff_z = total_green - total_red)
#         }) %>%
#         ggplot(aes(date)) +
#         geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
#         geom_line(aes(y = diff_z, group = factor(rep)), color = "gray50") +
#         geom_point(aes(y = diff_z, fill = diff_z), shape = 21, color = "gray40", size = 2) +
#         facet_grid(green ~ red) +
#         scale_fill_gradient2(low = "firebrick", high = "chartreuse", mid = "gray60",
#                              midpoint = 0, guide = FALSE) +
#         scale_y_continuous(expression("Scaled green" - "red"), breaks = c(-2, 0, 2)) +
#         scale_x_continuous("Date") + # , breaks = seq(0, 30, 10)) +
#         theme(strip.text = element_text(size = 10))
# }
#
# host_plot("alfalfa") %>%
#     ggsave(filename = "~/Desktop/comp_alfala.pdf", height = 6, width = 5)
# host_plot("clover") %>%
#     ggsave(filename = "~/Desktop/comp_clover.pdf", height = 6, width = 5)
# host_plot("fava") %>%
#     ggsave(filename = "~/Desktop/comp_fava.pdf", height = 6, width = 5)
#
#
#
#
# {host_df %>%
#     rename(green_line = green, red_line = red,
#            n_green = total_green, n_red = total_red) %>%
#     mutate(year = format(date, "%Y"),
#            month = format(date, "%m"),
#            day = format(date, "%d")) %>%
#     select(host, green_line, red_line, rep, year, month, day, n_green, n_red) %>%
#     identity()} %>%
#     write_csv(path = paste0("~/Box Sync/2019/host_shift/",
#                             "host_shift_data_entry.csv"))




host_df <- read_excel("under_constr/host_shift/host_shift_data_entry.xlsx") %>%
    mutate(date = as.Date(sprintf("%i-%i-%i", year, month, day)),
           combo = sprintf("%s__%s", green_line, red_line) %>%
               factor() %>% as.integer() %>% factor(),
           rep = as.integer(rep),
           red_line = gsub("√ò", "Ø", red_line)) %>%
    select(host, combo, rep, date, everything(), -year, -month, -day) %>%
    filter(! date %in% as.Date(c("2018-10-23", "2018-10-25")))



host_plot <- function(.host) {
    host_df %>%
        filter(host == .host) %>%
        group_by(combo, rep) %>%
        mutate(date = difftime(date, min(date), units = "days") %>% as.integer()) %>%
        ungroup() %>%
        split(interaction(.$combo, .$rep)) %>%
        map_dfr(function(.x) {
            .mean <- mean(c(.x$n_red, .x$n_green))
            .sd <- sd(c(.x$n_red, .x$n_green))
            if (is.na(.sd) || .sd == 0) .sd <- 1
            mutate(.x,
                   n_red = (n_red - .mean) / .sd,
                   n_green = (n_green - .mean) / .sd,
                   diff_z = n_green - n_red)
        }) %>%
        ggplot(aes(date)) +
        geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
        geom_line(aes(y = diff_z, group = factor(rep)), color = "gray50") +
        geom_point(aes(y = diff_z, fill = diff_z), shape = 21, color = "gray40", size = 2) +
        facet_grid(green_line ~ red_line) +
        scale_fill_gradient2(low = "firebrick", high = "chartreuse", mid = "gray60",
                             midpoint = 0, guide = FALSE) +
        scale_y_continuous(expression("Scaled green" - "red"), breaks = c(-2, 0, 2)) +
        scale_x_continuous("Date") + # , breaks = seq(0, 30, 10)) +
        theme(strip.text = element_text(size = 10))
}


host_plot("fava") %>%
    ggsave(filename = "~/Desktop/comp_fava.pdf", width = 5, height = 6)
host_plot("alfalfa") %>%
    ggsave(filename = "~/Desktop/comp_alfalfa.pdf", width = 5, height = 6)
host_plot("clover") %>%
    ggsave(filename = "~/Desktop/comp_clover.pdf", width = 5, height = 6)
