
library(tidyverse)
library(readxl)
library(lubridate)
library(clonewars)
library(viridisLite)
library(patchwork)
library(sf)
library(s2)
library(transformr)
library(ggsn)

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
        # This one had clover for part of it, but where isn't clear
        filter(Field != "349 Clover") %>%
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
               field = ifelse(field == "506SE", "506S", field),
               year = year(date),
               day = yday(date) - 1)  # <- previous years set Jan 1 as day 0
) %>%
    do.call(what = bind_rows) %>%
    # Not sure why, but there's a 2020 data point in the 2019 dataset:
    filter(year != 2020) %>%
    # We need to have at least 10 aphids dissected:
    filter(para_n >= 10) %>%
    mutate(cycle = floor(cycle),
           harvest = lead(cycle, default = tail(cycle, 1)) - cycle > 0 |
               # This is to force the first cell be `1`
               c(TRUE, rep(FALSE, n()-1)),
           year = factor(year, levels = sort(unique(year)))) %>%
    #'
    #' I found that these dates have the same exact numbers for all fields
    #' on dates two days before.
    #' These must be duplicates, so I'm removing them.
    #'
    filter(!date %in% as.Date(c("2011-07-14", "2011-07-21", "2011-08-26",
                                "2012-08-17", "2012-07-11", "2012-06-14")))



# Time series plot ----

#'
#'
#'
#' This essentially replicates (and updates) Nature E&E paper.
#'
#' I'm not using lines here bc I think points illustrate my point better,
#' and because code I used to separate by cycle for years 2011--2016 doesn't
#' appear to work for 2017--2019.
#'

par_ts_p <- par_df %>%
    split(.$year) %>%
    map_dfr(~ mutate(.x, field_col = factor(field) %>% as.integer())) %>%
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) %>%
    ggplot(aes(plot_date, para)) +
    geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
    geom_point(aes(color = field_col), alpha = 0.5, size = 1) +
    facet_wrap(~ year, ncol = 1) +
    scale_color_manual(values = viridis(10, begin = 0.1, end = 0.8) %>%
                           .[do.call(c, map(5:1, ~ c(.x, .x + 5)))],
                       guide = "none") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous("Proportion aphids parasitized", breaks = 0.4*0:2) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black"),
          strip.text = element_text(size = 9)) +
    NULL

# par_ts_p






# splitting into periods ----
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
#' obs_par_df %>%
#'     group_by(obs) %>%
#'     mutate(repeats = any(duplicated(field))) %>%
#'     ungroup() %>%
#'     filter(repeats) %>%
#'     distinct(year, obs)
#'
#' # # A tibble: 5 × 2
#' #   year    obs
#' #   <fct> <dbl>
#' # 1 2014    371
#' # 2 2015    495
#' # 3 2016    609
#' # 4 2016    623
#' # 5 2016    637
#' ```
#'

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
    # Plus I'm going to manually combine these dates because it's clear that
    # half the fields were sampled one day, half the next.
    mutate(obs = ifelse(date %in% as.Date(c("2016-06-13", "2016-06-14")),
                        mean(obs[date %in% as.Date(c("2016-06-13", "2016-06-14"))]),
                        obs)) %>%
    mutate(obs = factor(obs)) %>%
    # I'm going to want to filter out groups that only sampled few fields:
    group_by(year) %>%
    mutate(n_fields = length(unique(field))) %>%
    group_by(year, obs) %>%
    mutate(obs_n = n()) %>%
    # To simplify plotting by having a single obs number by the obs
    mutate(obs_date = round(mean(date)),
           obs_day = mean(day)) %>%
    ungroup()


#'
#' Plotting data through time where points are colored by obs.
#' Looks like the splitting worked well.
#'
# obs_par_df %>%
#     mutate(plot_date = as.Date(day, origin = "2022-01-01")) %>%
#     ggplot(aes(plot_date, para)) +
#     geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
#     geom_point(aes(color = obs), alpha = 0.5, size = 1) +
#     facet_wrap(~ year, ncol = 1) +
#     scale_color_manual(values = rep(c("#1b9e77", "#d95f02", "#7570b3"), 100),
#                        guide = "none") +
#     scale_x_date(date_breaks = "1 month", date_labels = "%b") +
#     scale_y_continuous("Proportion aphids parasitized", breaks = 0.4*0:2) +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(color = "black"),
#           panel.grid.major.x = element_line(color = "gray80"),
#           strip.text = element_text(size = 9)) +
#     NULL



# maps ----


# obs_dates to use for maps:
maps_dates <- c("2011-10-12",
                "2012-06-07",
                "2013-08-09",
                "2014-07-01",
                "2015-06-12",
                "2016-06-14",
                "2017-05-22",
                "2018-07-26",
                "2019-07-01") %>%
    as.Date()



fields_sf <- st_read(paste0("~/Box Sync/eco-evo_experiments/field-data/",
                            "arlington-fields/Arlington.gpkg")) %>%
    st_transform(st_crs(32616)) %>%
    mutate(geom = st_centroid(fields_sf$geom))


obs_fields_par <- obs_par_df %>%
    filter(obs_date %in% maps_dates) %>%
    mutate(obs = obs %>% fct_drop()) %>%
    split(1:nrow(.)) %>%
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf %>% filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$date <- .d$date
        .f$para <- .d$para
        .f$obs <- .d$obs
        .f$obs_date <- .d$obs_date
        return(.f)
    }) %>%
    do.call(what = rbind)


xy_lims <- st_bbox(obs_fields_par) %>% as.list()
xy_lims$xmin <- xy_lims$xmin - 400
xy_lims$xmax <- xy_lims$xmax + 400
xy_lims$ymin <- xy_lims$ymin - 400
xy_lims$ymax <- xy_lims$ymax + 400

fields_par_p_list <- map(
    sort(unique(obs_fields_par$obs_date)),
    function(.o) {
        obs_fields_par %>%
            filter(obs_date == .o) %>%
            ggplot() +
            geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
                      ymin = xy_lims$ymin, ymax = xy_lims$ymax,
                      fill = NA, color = "black") +
            # geom_sf(aes(size = para, fill = para), shape = 21, color = "black") +
            geom_sf(aes(size = para, color = para), shape = 16) +
            scale_fill_viridis_c("parasitism", option = "inferno",
                                 limits = c(0, 0.85), begin = 0.1, end = 0.9,
                                 aesthetics = c("color", "fill"),
                                 breaks = 0.2 * 0:4) +
            scale_size("parasitism", limits = c(0, 0.85), range = c(0.5, 8),
                       breaks = 0.2 * 0:4) +
            guides(fill = guide_legend(),
                   color = guide_legend(),
                   size = guide_legend()) +
            coord_sf(datum = st_crs(32616),
                     xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
                     ylim = as.numeric(xy_lims[c("ymin", "ymax")])) +
            labs(title = format(.o, "%e %b %Y")) +
            theme_void() +
            theme(plot.title = element_text(size = 10)) +
            # theme(legend.position = "none") +
            NULL
    })



# fields_par_p_list[[1]] <-
fields_par_p_list[[1]] +
    scalebar(dist = 0.5, dist_unit = "km", transform = FALSE,
             border.size = 0.5, st.dist = 0.05,
             st.size = 3, height = 0.01, location = "topleft",
             x.min = xy_lims$xmin + 200, x.max = xy_lims$xmax,
             y.min = xy_lims$ymin, y.max = xy_lims$ymax - 200)


fields_par_p <- wrap_plots(fields_par_p_list) +
    plot_layout(guides = "collect")
fields_par_p

fields_par_p_list[[1]] +
    scalebar(dist = 500, dist_unit = "m", transform = FALSE,
             border.size = 0.5,
             st.size = 3, height = 0.01, location = "topleft",
             x.min = xy_lims$xmin, x.max = xy_lims$xmax,
             y.min = xy_lims$ymin, y.max = xy_lims$ymax) +
    north(location = "topleft", symbol = 10,
             x.min = xy_lims$xmin, x.max = xy_lims$xmax,
             y.min = xy_lims$ymin, y.max = xy_lims$ymax)




wi_bounds <- paste0("~/Box Sync/eco-evo_experiments/field-data/",
                    "WI-boundary/Wisconsin_State_Boundary_24K.gpkg") %>%
    st_read() %>%
    st_transform(st_crs(32616))
wi_xy_lims <- st_bbox(wi_bounds) %>% as.list()


wi_bounds %>%
    ggplot() +
    geom_sf(size = 0.25) +
    geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
              ymin = xy_lims$ymin, ymax = xy_lims$ymax, fill = "red", color = NA) +
    coord_sf(datum = st_crs(32616)) +
    scalebar(dist = 50, dist_unit = "km", transform = FALSE,
             border.size = 0.5,
             st.size = 3, height = 0.01, location = "topleft",
             x.min = wi_xy_lims$xmin, x.max = wi_xy_lims$xmax,
             y.min = wi_xy_lims$ymin, y.max = wi_xy_lims$ymax) +
    north(location = "topright", symbol = 10,
          x.min = wi_xy_lims$xmin, x.max = wi_xy_lims$xmax,
          y.min = wi_xy_lims$ymin, y.max = wi_xy_lims$ymax) +
    theme_void() +
    NULL





# ============================================================================*
# ============================================================================*

# MAKING GIFS ----

# ============================================================================*
# ============================================================================*

library(gganimate)
library(gifski)





# Looking at effects of different periods of time to split by:
map_dfr(c("3 days", "1 week", "2 weeks", "3 weeks"),
        function(f) {
            par_df %>%
                arrange(date) %>%
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

#'
#' To combine the field polygons with the parasitism data, I first need to
#' prepare the parasitism data for plotting by 2-week intervals:
#'
gif_par_df <- obs_par_df %>%
    filter(obs_n >= 5)



#'
#' Now I can combine the parasitism date with the field polygons:
#'

gif_fields_par <- gif_par_df %>%
    split(1:nrow(.)) %>%
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf %>% filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$obs <- .d$obs
        .f$date <- .d$date
        .f$obs_date <- .d$obs_date
        .f$para <- .d$para
        return(.f)
    }) %>%
    do.call(what = rbind)

gif_fields_par
nrow(gif_fields_par) == nrow(gif_par_df)








gif_xy_lims <- rbind(
    map_dbl(gif_fields_par$geom, ~ .x[[1]]) %>% range(),
    map_dbl(gif_fields_par$geom, ~ .x[[2]]) %>% range())

gif_all_p <- gif_fields_par %>%
    ggplot() +
    geom_sf(aes(size = para, color = para)) +
    scale_color_viridis_c("parasitism",
                          limits = c(0, 0.85), begin = 0.2, end = 0.9,
                          breaks = 0.2 * 0:4) +
    scale_size("parasitism", limits = c(0, 0.85), range = c(0.5, 8),
               breaks = 0.2 * 0:4) +
    coord_sf(datum = st_crs(32616),
             xlim = gif_xy_lims[1,], ylim = gif_xy_lims[2,]) +
    ggtitle("Year: {year(as.Date(paste(current_frame)))}") +
    theme(plot.title = element_text(size = 10),
          # legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(color = guide_legend(), size = guide_legend()) +
    # Here comes the gganimate specific bits
    labs(subtitle = "Date: {format(as.Date(paste(current_frame)), '%d %b')}") +
    transition_manual(obs_date) +
    NULL



animate(gif_all_p, renderer = gifski_renderer(),
        nframes = length(levels(gif_fields_par$obs)), fps = 10)

anim_save("~/Desktop/field_gifs/fields_2011--2019.gif",
          animate(gif_all_p, renderer = gifski_renderer(),
                  nframes = length(levels(gif_fields_par$obs)), fps = 10))



# gifs by year:
for (y in levels(gif_fields_par$year)) {

    .d <- gif_fields_par %>%
        filter(year == y) %>%
        mutate(obs_int = obs %>% fct_drop() %>% as.integer())

    gif_xy_lims_y <- rbind(map_dbl(.d$geom, ~ .x[[1]]) %>% range(),
                           map_dbl(.d$geom, ~ .x[[2]]) %>% range())

    gif_p <- .d %>%
        ggplot() +
        geom_sf(aes(size = para, color = para)) +
        scale_color_viridis_c("parasitism",
                              limits = c(0, 0.85), begin = 0.2, end = 0.9,
                              breaks = 0.2 * 0:4) +
        scale_size("parasitism", limits = c(0, 0.85), range = c(0.5, 8),
                   breaks = 0.2 * 0:4) +
        coord_sf(datum = st_crs(32616),
                 xlim = gif_xy_lims[1,], ylim = gif_xy_lims[2,]) +
        ggtitle(paste("Year:", y)) +
        theme(plot.title = element_text(size = 10),
              # legend.position = "none",
              axis.text = element_blank(),
              axis.ticks = element_blank()) +
        guides(color = guide_legend(), size = guide_legend()) +
        # Here comes the gganimate specific bits
        labs(subtitle = "Date: {format(as.Date(paste(current_frame)), '%d %b')}") +
        transition_manual(obs_date) +
        NULL

    anim_save(sprintf("~/Desktop/field_gifs/fields_%s.gif", y),
              animate(gif_p, renderer = gifski_renderer(),
                      nframes = max(.d$obs_int), fps = 2))

    cat(sprintf("%s done\n", y))

}; rm(y, .d, gif_xy_lims_y, gif_p)




