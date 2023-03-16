
library(tidyverse)
library(gameofclones)
library(ggtext)             # element_markdown
library(lubridate)          # yday, year
library(viridisLite)        # inferno
library(here)               # here
library(grid)               # grid.newpage, grid.draw
library(sf)                 # st_* (e.g., st_read, st_transform, st_crs)
library(ggmap)              # get_googlemap
library(rnaturalearth)      # ne_countries
library(rnaturalearthdata)  # used for ne_countries


source(".Rprofile")







#'
#' Much of the data here are from https://doi.org/10.6084/m9.figshare.11828865.v1
#'
#' To run the scripts below, download this dataset, then rename
#' `Ives et al. 2020 Data Fig2_1.csv` to `parasitism-2001-2016.csv`
#' and
#' `Ives et al. 2020 Data Fig3A.csv` to `symbionts-2012-2017.csv`
#'
#' Then put both inside the `gameofclones/_results/_data` folder.
#'




#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Parasitism dataset - read and organize ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------


par_df <- list(
    here("_results/_data/parasitism-2001-2016.csv") |>
        read_csv(col_types = cols()) |>
        select(field, Cycle, DateFormat, year, day, para, paraN) |>
        rename(cycle = Cycle, para_n = paraN, date = DateFormat) |>
        mutate(field = paste(field),
               # Field 18 is the same as field 110; the numbers weren't fixed
               # until 2015
               field = ifelse(field == "18", "110", field),
               # These were re-coded to stay as numbers, but I'm converting them
               # back to their original names:
               field = ifelse(field == "506.2", "506N", field),
               field = ifelse(field == "506.1", "506S", field),
               date = as.Date(date)) |>
        filter(!is.na(para), !is.na(para_n)),
    list.files(here("_results/_data"), "parasitism-....\\.csv",
               full.names = TRUE) |>
        map_dfr(read_csv, show_col_types = FALSE) |>
        # I can't find this one in any of the Arlington maps...
        filter(Field != "N1902") |>
        # This one had clover for part of it, but where isn't clear
        filter(Field != "349 Clover") |>
        rename_with(tolower) |>
        select(field, cycle, date, starts_with("diss")) |>
        rename_with(function(.x) gsub("^diss_", "", gsub("\\+", "_", .x))) |>
        mutate(para = g_para + g_p_f + r_para + r_p_f,
               para_n = para + g_unpara + g_fungus + r_unpara + r_fungus,
               para = para / para_n) |>
        filter(!is.na(para)) |>
        select(field, cycle, date, para, para_n) |>
        # Because two date formats are used:
        mutate(date1 = as.Date(date, format = "%m/%d/%y"),
               date2 = as.Date(date, format = "%d-%b-%y"),
               # `structure(...` below is to keep `ifelse` from
               # changing to numeric
               date = structure(ifelse(is.na(date1), date2, date1),
                                class = class(date2))) |>
        select(-date1, -date2) |>
        mutate(field = paste(field),
               field = ifelse(field == "506SE", "506S", field),
               year = year(date),
               day = yday(date) - 1)  # <- previous years set Jan 1 as day 0
) |>
    do.call(what = bind_rows) |>
    # Not sure why, but there's a 2020 data point in the 2019 dataset:
    filter(year != 2020) |>
    # We need to have at least 10 aphids dissected:
    filter(para_n >= 10) |>
    mutate(cycle = floor(cycle),
           harvest = lead(cycle, default = tail(cycle, 1)) - cycle > 0 |
               # This is to force the first cell be `1`
               c(TRUE, rep(FALSE, n()-1)),
           year = factor(year, levels = sort(unique(year)))) |>
    #'
    #' I found that these dates have the same exact numbers for all fields
    #' on dates two days before.
    #' These must be duplicates, so I'm removing them.
    #'
    filter(!date %in% as.Date(c("2011-07-14", "2011-07-21", "2011-08-26",
                                "2012-08-17", "2012-07-11", "2012-06-14"))) |>
    # Add the relative fitness for resistance (r_r / r_s)
    mutate(rr_rs = rel_res_fitness(para))


# Color palettes to use in plots (using binned rr_rs values):
rr_rs_pal <- with(list(pal = inferno, inds = c(80, 60, 20)),
                  list(color = c(pal(100)[inds], "gray70"),
                       fill = c(pal(100)[inds], "white"),
                       fill2 = c(pal(100, alpha = 0.5)[inds], "white")))

# Add factor that breaks rr_rs into bins for plotting:
add_rr_rs_fct <- function(.df) {
    .df |>
    mutate(rr_rs_fct = cut(rr_rs, c(0, 1, 1.05, 1.1, Inf),
                           labels = c("< 1", "1 – 1.05", "1.05 – 1.1", "> 1.1")),
           rr_rs_fct = factor(rr_rs_fct, levels = rev(levels(rr_rs_fct))))
}

#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))

#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Split parasitism into periods ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' To combine the field polygons with the parasitism data, I first need to
#' prepare the parasitism data for plotting by observation date, where
#' the dates can be a bit different:
#'
obs_par_df <- par_df |>
    select(-cycle, -harvest) |>
    arrange(date) |>
    mutate(obs = cut.Date(date, breaks = "3 days", labels = FALSE) |>
               as.numeric())

#'
#' These observation groups had multiple observations of the same fields,
#' so I split them up further.
#'
#' ```
#' obs_par_df |>
#'     group_by(obs) |>
#'     mutate(repeats = any(duplicated(field))) |>
#'     ungroup() |>
#'     filter(repeats) |>
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

obs_par_df <- obs_par_df |>
    group_by(obs) |>
    mutate(repeats = any(duplicated(field))) |>
    ungroup() |>
    mutate(tmpid = interaction(year, obs, drop = TRUE)) |>
    split(~ tmpid) |>
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
    }) |>
    select(-repeats, -tmpid) |>
    # Plus I'm going to manually combine these dates because it's clear that
    # half the fields were sampled one day, half the next.
    mutate(obs = ifelse(date %in% as.Date(c("2016-06-13", "2016-06-14")),
                        mean(obs[date %in% as.Date(c("2016-06-13", "2016-06-14"))]),
                        obs)) |>
    mutate(obs = factor(obs)) |>
    # I'm going to want to filter out groups that only sampled few fields:
    group_by(year) |>
    mutate(n_fields = length(unique(field))) |>
    group_by(year, obs) |>
    mutate(obs_n = n()) |>
    # To simplify plotting by having a single obs number by the obs
    mutate(obs_date = round(mean(date)),
           obs_day = mean(day)) |>
    ungroup()


#'
#' Plotting data through time where points are colored by obs.
#' Looks like the splitting worked well.
#'
# obs_par_df |>
#     mutate(plot_date = as.Date(day, origin = "2022-01-01")) |>
#     ggplot(aes(plot_date, para)) +
#     geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
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





#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Parasitism maps ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------




fields_sf <- st_read(here("_results/_data/Arlington.geojson")) |>
    st_transform(st_crs(3857)) |>
    mutate(geometry = st_centroid(geometry)) |>
    rename(geom = geometry)


obs_fields_par <- obs_par_df |>
    filter(obs_date %in% maps_dates) |>
    # Make sure all fields are present in all dates within year:
    group_by(year, field) |>
    mutate(no = n()) |>
    ungroup() |>
    filter(no == 3) |>
    mutate(obs = obs |> fct_drop(),
           tmpid = 1:n()) |>
    split(~tmpid) |>
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf |> filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$date <- .d$date
        .f$para <- .d$para
        .f$para_n <- .d$para_n
        .f$rr_rs <- .d$rr_rs
        .f$obs <- .d$obs
        .f$obs_date <- .d$obs_date
        return(.f)
    }) |>
    do.call(what = rbind) |>
    mutate(plot_date = factor(paste(obs_date),
                              levels = paste(sort(unique(obs_date))),
                              labels = format(sort(unique(obs_date)),
                                              "%e %b")),
           plot_date_by_yr = (as.integer(plot_date) - 1) %% 3 + 1,
           plot_date_by_yr = factor(plot_date_by_yr)) |>
    add_rr_rs_fct()


xy_lims <- st_bbox(obs_fields_par) |>
    as.list() |> as_tibble() |>
    mutate(xmin = xmin - 400,
           xmax = xmax + 400,
           ymin = ymin - 400,
           ymax = ymax + 400)

fields_par_scale_df <- st_bbox(obs_fields_par) |>
    as.list() |> as_tibble() |>
    mutate(xmax = xmin + 1000,
           ymin = ymax - 100,
           ymax = ymin + 200,
           plot_date = sort(obs_fields_par$plot_date)[1],
           year = "2013", plot_date_by_yr = "1")


fields_par_p <- obs_fields_par |>
    ggplot() +
    geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
              ymin = xy_lims$ymin, ymax = xy_lims$ymax,
              fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(aes(size = para, color = rr_rs_fct, fill = rr_rs_fct), shape = 21) +
    # # ------*
    # # DIY scale bar:
    # geom_rect(data = fields_par_scale_df,
    #              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #              fill = "black", color = NA) +
    # geom_text(data = fields_par_scale_df |>
    #               mutate(x = (xmin + xmax) / 2, y = ymin - 200),
    #           aes(x = x, y = y), label = "1 km",
    #           size = 9 / 2.83465, vjust = 1) +
    # # ------*
    scale_color_manual(NULL, guide = "none",
                       values = rr_rs_pal$color) +
    scale_fill_manual(NULL, guide = "none",
                      values = rr_rs_pal$fill) +
    scale_size("Proportion\nparasitized",
               limits = c(0, 0.85), range = c(0.5, 8),
               breaks = 0.2 * 0:4) +
    guides(size = guide_legend(override.aes = list(shape = 16))) +
    coord_sf(datum = st_crs(3857),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")]),
             clip = "off") +
    facet_wrap(~ plot_date, nrow = 2) +
    theme_void() +
    theme(strip.text = element_text(size = 9, margin = margin(0,0,b=3,t=3)),
          legend.position = "none")


fields_par_p

save_plot(here("_results/_plots/field-data/par-map.pdf"), fields_par_p,
          w = 3.5, h = 2.5)





#' ----------------------
# Parasitism map separate legend
#' ----------------------


fields_par_p_leg <- function() {
    legend <- (fields_par_p + theme(legend.position = "right")) |>
        (function(a.gplot){
            tmp <- ggplot_gtable(ggplot_build(a.gplot))
            leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
            legend <- tmp$grobs[[leg]]
            legend
        })()
    grid.newpage()
    grid.draw(legend)
}

save_plot(here("_results/_plots/field-data/par-map-legend.pdf"), fields_par_p_leg,
          w = 2, h = 3.5)


#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# USA map inset ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------


usa_map_p <- ne_countries(country = "united states of america",
                          scale = "medium", returnclass = "sf") |>
    st_cast("POLYGON") |>
    (\(polies) {
        #' Find the polygon whose western edge is the furthest west without
        #' being in the eastern hemisphere. This should be the contiguous USA.
        bbx = lapply(1:nrow(polies),
                     \(i) st_bbox(polies[i,])[["xmax"]]) |>
            do.call(what = c)
        i <- which(bbx == max(bbx[bbx < 0]))
        return(polies[i,])
    })() |>
    st_transform(st_crs(3857)) |>
    ggplot() +
    geom_sf(size = 0.25) +
    geom_point(data = tibble(x = (xy_lims$xmin + xy_lims$xmax) / 2,
                             y = (xy_lims$ymin + xy_lims$ymax) / 2),
               aes(x, y), color = "black", shape = 20, size  = 4) +
    coord_sf(datum = st_crs(3857)) +
    theme_void()

# usa_map_p

save_plot(here("_results/_plots/field-data/par-map-usa-inset.pdf"), usa_map_p,
          w = 2.5, h = 1.5)

# If you have pdfcrop installed:
# system(paste("pdfcrop", here("_results/_plots/field-data/par-map-usa-inset.pdf")))






#' For WI boundary, first download GeoJSON from here:
#' https://data-wi-dnr.opendata.arcgis.com/datasets/wi-dnr::wisconsin-state-boundary-24k/explore?location=44.724029%2C-89.836300%2C8.61
#' (I accessed this on 15 March 2023)
#'
#' I then reduced the resolution and wrote this to a smaller file using the
#' following code in R:
#' ```
#' library(sf)
#' library(here)
#' wi_bounds_hires <- "Wisconsin_State_Boundary_24K.geojson" |>
#'     st_read() |>
#'     st_transform(st_crs(3857))
#' wi_bounds <- wi_bounds_hires |> st_simplify(preserveTopology = FALSE, dTolerance = 1000)
#' st_write(wi_bounds, here("_results/_data/WI-boundary.geojson"))
#' ```


wi_bounds  <- here("_results/_data/WI-boundary.geojson") |>
    st_read()


wi_inset <- wi_bounds |>
    ggplot() +
    geom_sf(size = 0.25) +
    geom_point(data = tibble(x = (xy_lims$xmin + xy_lims$xmax) / 2,
                             y = (xy_lims$ymin + xy_lims$ymax) / 2),
               aes(x, y), color = "black", shape = 20, size  = 4) +
    coord_sf(datum = st_crs(3857)) +
    theme_void()

# wi_inset


save_plot(here("_results/_plots/field-data/par-map-wi-inset.pdf"), wi_inset,
          w = 2, h = 2)


#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Satellite map ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

arl_google <- get_googlemap(c(-89.35399, 43.30917), zoom = 13, maptype = "satellite")

arl_goog_p <- ggmap(arl_google) +
    # geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
    #           ymin = xy_lims$ymin, ymax = xy_lims$ymax,
    #           fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(datum = st_crs(3857),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")]),
             clip = "on") +
    theme_void()


save_plot(here("_results/_plots/field-data/par-map-google-inset.pdf"), arl_goog_p,
          w = 1.5, h = 1.5)




