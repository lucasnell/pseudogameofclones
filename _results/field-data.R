
library(tidyverse)
library(readxl)
library(lubridate)
library(clonewars)
library(viridisLite)

source(".Rprofile")




par_df <- list(
    read_excel(paste0("~/Box Sync/eco-evo_experiments/field-data/",
                      "field-data.xlsx"), sheet = "parasitism",
               na = c("", "NA")) %>%
        select(field, Cycle, DateFormat, year, day, para, paraN) %>%
        rename(cycle = Cycle, para_n = paraN, date = DateFormat) %>%
        mutate(field = paste(field),
               # Field 18 is the same as field 110; the numbers weren't fixed
               # until 2015
               field = ifelse(field == "18", "110", field),
               # These were re-coded to stay as numbers, but I'm converting them
               # back to their original names:
               field = ifelse(field == "506.2", "506N", field),
               field = ifelse(field == "506.1", "506S", field),
               date = as.Date(date)) %>%
        filter(!is.na(para), !is.na(para_n)),
    list.files("~/Box Sync/eco-evo_experiments/field-data/2017--2019", "*.csv",
               full.names = TRUE) %>%
        map_dfr(read_csv, show_col_types = FALSE) %>%
        # I can't find this one in any of the Arlington maps...
        filter(Field != "N1902") %>%
        rename_with(tolower) %>%
        select(field, cycle, date, starts_with("diss")) %>%
        rename_with(function(.x) gsub("^diss_", "", gsub("\\+", "_", .x))) %>%
        mutate(para = g_para + g_p_f + r_para + r_p_f,
               para_n = para + g_unpara + g_fungus + r_unpara + r_fungus,
               para = para / para_n) %>%
        filter(!is.na(para)) %>%
        select(field, cycle, date, para, para_n) %>%
        # Because two date formats are used:
        mutate(date1 = as.Date(date, format = "%m/%d/%y"),
               date2 = as.Date(date, format = "%d-%b-%y"),
               # `structure(...` below is to keep `ifelse` from
               # changing to numeric
               date = structure(ifelse(is.na(date1), date2, date1),
                                class = class(date2))) %>%
        select(-date1, -date2) %>%
        mutate(field = paste(field),
               # Cleaning up some field names
               field = ifelse(field == "349 Clover", "349", field),
               field = ifelse(field == "506SE", "506S", field),
               year = year(date),
               day = yday(date) - 1)  # <- previous years set Jan 1 as day 0
) %>%
    do.call(what = bind_rows) %>%
    # Not sure why, but there's a 2020 data point in the 2019 dataset:
    filter(year != 2020) %>%
    mutate(cycle = floor(cycle),
           harvest = lead(cycle, default = tail(cycle, 1)) - cycle > 0 |
               # This is to force the first cell be `1`
               c(TRUE, rep(FALSE, n()-1)),
           year = factor(year, levels = sort(unique(year))))




# This essentially replicates Nature E&E paper:
par_df %>%
    mutate(cycle_grp = interaction(year, field, cycle, drop = TRUE)) %>%
    split(.$year) %>%
    map_dfr(~ mutate(.x, field_col = factor(field) %>% as.integer())) %>%
    mutate(field_col = factor(field_col)) %>%
    ggplot(aes(day, para)) +
    geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
    # geom_line(aes(color = field_col, group = cycle_grp)) +
    geom_point(aes(color = field_col), alpha = 0.5, size = 1) +
    facet_wrap(~ year, ncol = 1) +
    scale_color_manual(values = viridis(10, begin = 0.2, end = 0.9) %>%
                           .[do.call(c, map(5:1, ~ c(.x, .x + 5)))],
                       guide = "none") +
    # scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    # theme(axis.title.x = element_blank(),
    #       axis.text.x = element_text(size = 10, color = "black",
    #                                  angle = 45, hjust = 0.7)) +
    theme(strip.text = element_text(size = 8)) +
    xlab("Day of year") +
    ylab("Parasitism") +
    NULL



# Looking at effects of different periods of time to split by:
map_dfr(c("3 days", "1 week", "2 weeks", "3 weeks"),
        function(f) {
            par_df %>%
                mutate(period = cut.Date(date, breaks = f)) %>%
                group_by(year, period, field) %>%
                summarize(para = mean(para), .groups = "drop") %>%
                group_by(year) %>%
                mutate(total_f = length(unique(field))) %>%
                group_by(year, period) %>%
                summarize(obs_p = length(unique(field)) / total_f[1],
                          .groups = "drop") %>%
                group_by(year) %>%
                summarize(n_all = sum(obs_p == 1), .groups = "drop") %>%
                mutate(period = gsub(" ", "_", f))
        }) %>%
    pivot_wider(id_cols = c(year, period), names_from = period, values_from = n_all)


# Splitting into averages over 2-week period:

par_df %>%
    mutate(period = cut.Date(date, breaks = "2 weeks")) %>%
    group_by(year, period, field) %>%
    summarize(para = mean(para), .groups = "drop") %>%
    group_by(year) %>%
    mutate(total_f = length(unique(field))) %>%
    group_by(year, period) %>%
    summarize(obs_p = length(unique(field)) / total_f[1],
              .groups = "drop") %>%
    mutate(period = period %>% paste() %>% as.Date(),
           days = yday(period) - 1,
           year = year(period) %>% factor()) %>%
    ggplot(aes(days, obs_f / total_f)) +
    geom_point(alpha = 0.5, size = 1) +
    facet_wrap(~ year, ncol = 1) +
    theme(strip.text = element_text(size = 8)) +
    xlab("Day of year") +
    scale_y_continuous("Proportion of fields", breaks = 0.2 * 0:5)  +
    theme(panel.grid.major.y = element_line(color = "gray80"))


par_df %>%
    mutate(cycle_grp = interaction(year, field, cycle, drop = TRUE)) %>%
    split(.$year) %>%
    map_dfr(~ mutate(.x, field_col = factor(field) %>% as.integer())) %>%
    mutate(field_col = factor(field_col)) %>%
    mutate(week2 = cut.Date(date, breaks = "2 weeks")) %>%
    group_by(year, week2, field) %>%
    summarize(para = mean(para), .groups = "drop") %>%
    mutate(week2 = week2 %>% paste() %>% as.Date(),
           days = yday(week2) - 1) %>%
    split(.$year) %>%
    map_dfr(~ mutate(.x, field_col = factor(field) %>% as.integer())) %>%
    mutate(field_col = factor(field_col)) %>%
    ggplot(aes(days, para)) +
    geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
    # geom_line(aes(color = field_col, group = cycle_grp)) +
    geom_point(aes(color = field_col), alpha = 0.5, size = 1) +
    facet_wrap(~ year, ncol = 1) +
    scale_color_manual(values = viridis(10, begin = 0.2, end = 0.9) %>%
                           .[do.call(c, map(5:1, ~ c(.x, .x + 5)))],
                       guide = "none") +
    # scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    # theme(axis.title.x = element_blank(),
    #       axis.text.x = element_text(size = 10, color = "black",
    #                                  angle = 45, hjust = 0.7)) +
    theme(strip.text = element_text(size = 8)) +
    xlab("Day of year") +
    ylab("Parasitism") +
    NULL





library(sf)
library(gganimate)
library(gifski)
library(transformr)

fields_sf <- st_read(paste0("~/Box Sync/eco-evo_experiments/field-data/",
                         "arlington-fields/Arlington.gpkg")) %>%
    st_transform(st_crs(32616))


#'
#' To combine the field polygons with the parasitism data, I first need to
#' prepare the parasitism data for plotting by observation date, where
#' the dates can be a bit different:
#'
obs_par_df <- par_df %>%
    select(-cycle, -harvest) %>%
    arrange(date) %>%
    mutate(obs = cut.Date(date, breaks = "3 days", labels = FALSE) %>%
               as.numeric())

#'
#' These observation groups had multiple observations of the same fields,
#' so I split them up further.
#'
#' ```
#' 1 2014    371
#' 2 2015    495
#' 3 2016    609
#' 4 2016    623
#' 5 2016    637
#' ```

obs_par_df <- obs_par_df %>%
    group_by(obs) %>%
    mutate(repeats = any(duplicated(field))) %>%
    ungroup() %>%
    split(interaction(.$year, .$obs, drop = TRUE)) %>%
    map_dfr(function(.dd) {
        unq_dates <- length(unique(.dd$date))
        if (.dd$repeats[1]) {
            stopifnot(unq_dates > 1)
            for (i in 2:unq_dates) {
                frac <- (i - 1) / unq_dates
                .dd$obs[.dd$date == unique(.dd$date)[i]] <- .dd$obs[1] + frac
            }
        }
        return(.dd)
    }) %>%
    select(-repeats) %>%
    mutate(obs = factor(obs)) %>%
    # I'm going to want to filter out groups that only sampled few fields:
    group_by(year) %>%
    mutate(n_fields = length(unique(field))) %>%
    group_by(year, obs) %>%
    mutate(obs_n = n()) %>%
    # To simplify plotting by having a single obs number by the obs
    mutate(obs_date = mean(date),
           obs_day = mean(day)) %>%
    ungroup()


#'
#' Plotting data through time where points are colored by obs.
#' Looks like the splitting worked well.
#'
# obs_par_df %>%
#     ggplot(aes(day, para)) +
#     geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
#     geom_point(aes(color = obs), alpha = 0.5, size = 1) +
#     facet_wrap(~ year, ncol = 1) +
#     scale_color_manual(values = rep(c("green", "blue", "red"), 50),
#                        guide = "none") +
#     theme(strip.text = element_text(size = 8)) +
#     xlab("Day of year") +
#     ylab("Parasitism") +
#     NULL


#'
#' Now I can combine them:
#'

fields_par <- obs_par_df %>%
    filter(obs_n >= n_fields / 2) %>%
    # filter(obs_n == n_fields) %>%
    split(1:nrow(.)) %>%
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf %>% filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$date <- .d$date
        .f$obs <- .d$obs
        .f$obs_date <- .d$obs_date
        .f$obs_day <- .d$obs_day
        .f$para <- .d$para
        return(.f)
    }) %>%
    do.call(what = rbind)

fields_par %>% head()
fields_par %>% nrow()




# gifs by year:
for (y in levels(fields_par$year)) {

    .d <- fields_par %>%
        filter(year == y) %>%
        mutate(obs_int = obs %>% fct_drop() %>% as.integer())

    xy_lims <- cbind(
        map(.d$geom, ~ apply(.x[[1]], 2, max)) %>%
            do.call(what = rbind) %>%
            apply(2, max),
        map(.d$geom, ~ apply(.x[[1]], 2, min)) %>%
            do.call(what = rbind) %>%
            apply(2, min))

    gif_p <- .d %>%
        ggplot() +
        geom_sf(aes(fill = para), color = NA) +
        scale_fill_viridis_c(limits = c(0, 0.85), begin = 0.2, end = 0.9) +
        # theme(axis.text = element_blank(), axis.ticks = element_blank()) +
        theme_void() +
        coord_sf(datum = st_crs(32616), xlim = xy_lims[1,], ylim = xy_lims[2,]) +
        ggtitle(paste("Year:", y)) +
        # Here comes the gganimate specific bits
        labs(subtitle = "Observation: {current_frame}") +
        transition_manual(obs_int) +
        NULL

    anim_save(sprintf("~/Desktop/field_gifs/fields_%s.gif", y),
              animate(gif_p, renderer = gifski_renderer(),
                      nframes = max(.d$obs_int), fps = 2))

    cat(sprintf("%s done\n", y))

}; rm(y, .d, xy_lims, gif_p)




