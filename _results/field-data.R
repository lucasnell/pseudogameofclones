
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
library(parallel)

source(".Rprofile")


options(mc.cores = max(1, detectCores()-2))



#'
#' Calculate r_r / r_s (relative fitness for resistance vs susceptible aphids).
#' This is almost directly from code in Ives et al. (2020).
#'
#' This function will be very slow if used inside an `apply` function.
#' Instead, `para` should be a vector.
#'
rel_res_fitness <- function(para, p_res = 0.48) {

    stopifnot(is.numeric(para) && is.null(dim(para)))
    stopifnot(is.numeric(p_res) && length(p_res) == 1)

    if (file.exists("_results/rr_rs_lookup.rds")) {

        rr_rs_df <- readRDS("_results/rr_rs_lookup.rds")

    } else {

        # We first create lookup table based on Ives et al. (2020) code
        growthcost = 6.5266912e-01
        resist = 2.7982036e-01
        Leslie_s <- matrix(c(0.5, 0, 0, 0, 2.55, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5,
                             0, 0, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0.5, 0.8),
                           byrow = TRUE, nrow=5)
        Leslie_r <- Leslie_s
        Leslie_r[1,5] <- growthcost * Leslie_r[1,5]
        relatt <- c(0.1200, 0.2700, 0.3900, 0.1600, 0.0600)
        kk <- 0.3475
        yLeslie <- matrix(c(0.5, 0, 0, 0, 0,
                            0.5, 0.5, 0, 0, 0,
                            0, 0.5, 0.5, 0, 0,
                            0, 0, 0.5, 0.667, 0,
                            0, 0, 0, 0.333, 0.95),
                          byrow = TRUE, nrow = 5)
        one_calc <- function(a){
            A <- (1 + exp(a) * relatt / kk)^(-kk)
            rs <- max(abs(eigen(diag(A) %*% Leslie_s)$values))
            rr <- max(abs(eigen(diag(1 - resist * (1-A)) %*% Leslie_r)$values))
            B11 <- diag(A) %*% Leslie_s
            B21 <- matrix(0,nrow=5, ncol=5)
            B21[1,] <- 1-A
            B12 <- matrix(0,nrow=5, ncol=5)
            B22 <- yLeslie
            B <- rbind(cbind(B11, B12), cbind(B21, B22))
            #colSums(B)
            SADA <- eigen(B)$vectors[,1]
            ps <- Re(sum(SADA[8:9])/sum(SADA[c(4:5,8:9)]))
            if (is.nan(ps)) ps <- 1
            B11 <- diag(1 - resist*(1-A)) %*% Leslie_r
            B21 <- matrix(0,nrow=5, ncol=5)
            B21[1,] <- resist*(1-A)
            B12 <- matrix(0,nrow=5, ncol=5)
            B22 <- yLeslie
            B <- rbind(cbind(B11, B12), cbind(B21, B22))
            SADA <- eigen(B)$vectors[,1]
            pr <- Re(sum(SADA[8:9])/sum(SADA[c(4:5,8:9)]))
            rr_rs <- rr / rs
            return(c(a, rs, rr, ps, pr, rr_rs))
        }

        # Takes ~15 sec on my machine.
        rr_rs_list <- mclapply(round(seq(-20, 5, 1e-4), 4), one_calc)
        rr_rs_df <- as.data.frame(do.call(rbind, rr_rs_list))
        colnames(rr_rs_df) <- c("a", "rs", "rr", "ps", "pr", "rr_rs")
        # saveRDS(rr_rs_df, "rr_rs_lookup.rds")
    }

    # Now we look up values of relative fitnesses based on the attack rate
    # associated with the parasitism reported (and the assumed proportion of
    # aphids that are currently resistant)
    para_lookup <- (1 - p_res) * rr_rs_df$ps + p_res * rr_rs_df$pr
    rr_rs <- sapply(para, function(p) {
        nearest <- which(abs(para_lookup - p) < 1.5e-5)

        # dealing with know discontinuity when p_res = 0.48 (default):
        if (length(nearest) == 0 &&
            p > para_lookup[rr_rs_df$a == 1.0581] &&
            p < para_lookup[rr_rs_df$a == 1.0582]) {
            nearest <- which(rr_rs_df$a == 1.0581 | rr_rs_df$a == 1.0582)
        }
        if (length(nearest) > 0) {
            return(mean(rr_rs_df$rr_rs[nearest], na.rm=TRUE))
        } else stop("cannot find fit for p")
    })

    return(rr_rs)

}



# Main dataset:

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
                                "2012-08-17", "2012-07-11", "2012-06-14"))) %>%
    # Add the relative fitness for resistance (r_r / r_s)
    mutate(rr_rs = rel_res_fitness(para))


# Color palette to use in plots (using binned rr_rs values):
rr_rs_pal <- c(viridis(3, begin = 0.1, end = 0.8, direction = -1), "gray70")
# Add factor that breaks rr_rs into bins for plotting:
add_rr_rs_fct <- function(.df) {
    .df %>%
    mutate(rr_rs_fct = cut(rr_rs, c(0, 1, 1.05, 1.1, Inf),
                           labels = c("< 1", "1 – 1.05", "1.05 – 1.1", "> 1.1")),
           rr_rs_fct = factor(rr_rs_fct, levels = rev(levels(rr_rs_fct))))
}

#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))



#' #'
#' #' Thresholds are the "... level of parasitism above which selection favours
#' #' resistant Hamiltonella–APSE3 clones" for...
#' #' - `min` minimum observed Hamiltonella frequency (0.02)
#' #' - `mean` mean observed Hamiltonella frequency (0.48)
#' #' - `max` max observed Hamiltonella frequency (0.88)
#' #'
#' #' For simplicity, I'm only using `mean` in plots below.
#' #'
#' #' These are from Ives et al. (2020)
#' #'
#' para_thresh = list(min = 0.13, mean = 0.21, max = 0.30)
#'
#' # Color palette to use in plots:
#' thresh_pal <- c("gray70", viridis(3, begin = 0.1, end = 0.8)) %>% rev()
#'
#'
#' #'
#' #' To add whether parasitism goes beyond each threshold.
#' #'
#' add_thresh <- function(.df) {
#'     .df %>%
#'         mutate(thresh_grp = case_when(para < para_thresh$min ~ 1,
#'                                       para < para_thresh$mean ~ 2,
#'                                       para < para_thresh$max ~ 3,
#'                                       TRUE ~ 4) %>%
#'                    factor(levels = 4:1,
#'                           labels = c("very often", "often", "rare", "very rare")))
#'                           # labels = c("p[res] <= 0.02", "p[res] <= 0.48",
#'                           #            "p[res] <= 0.88", "p[res] > 0.88")))
#' }




# Time series plot ----

#'
#' This essentially replicates (and updates) Nature E&E paper.
#'
#' I'm not using lines here bc I think points illustrate my point better,
#' and because code I used to separate by cycle for years 2011--2016 doesn't
#' appear to work for 2017--2019.
#'
#' To show the
#'

par_ts_p <- par_df %>%
    split(.$year) %>%
    map_dfr(~ mutate(.x, field_col = factor(field) %>% as.integer())) %>%
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) %>%
    add_rr_rs_fct() %>%
    ggplot(aes(plot_date, para)) +
    geom_hline(yintercept = 0, color = "gray70", size = 0.5) +
    # geom_hline(yintercept = para_thresh$min, color = "gray70",
    #            size = 0.5, linetype = 3) +
    # geom_hline(yintercept = para_thresh$max, color = "gray70",
    #            size = 0.5, linetype = 3) +
    # geom_hline(yintercept = para_thresh$mean, color = "gray70", size = 0.5) +
    geom_segment(data = tibble(plot_date = yday(maps_dates) %>%
                                   as.Date(origin = "2021-12-31"),
                               year = year(maps_dates) %>% factor(),
                               para = -0.1),
               aes(xend = plot_date, yend = para + 0.05),
               size = 0.5, linejoin = "mitre",
               arrow = arrow(length = unit(0.1, "lines"), type = "closed")) +
    # geom_point(aes(color = field_col), alpha = 0.5, size = 1) +
    # geom_point(aes(color = thresh_grp), alpha = 0.5, size = 1) +
    geom_point(aes(color = rr_rs_fct), alpha = 0.5, size = 1) +
    facet_wrap(~ year, ncol = 3) +
    # scale_color_manual(values = viridis(10, begin = 0.1, end = 0.8) %>%
    #                        .[do.call(c, map(5:1, ~ c(.x, .x + 5)))],
    #                    guide = "none") +
    # scale_color_manual("Selection for\nresistance", values = thresh_pal) +  # ,
    #                    # labels = str2expression) +
    scale_color_manual(paste0("Relative fitness for\nresistant aphids\n",
                              "(r\U1D63 / r\U209B)"),
                      values = rr_rs_pal) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous("Parasitism", breaks = 0.4*0:2) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          legend.title = element_text(hjust = 0),
          strip.text = element_text(size = 9)) +
    NULL

par_ts_p








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


fields_sf <- st_read(paste0("~/Box Sync/eco-evo_experiments/field-data/",
                            "arlington-fields/Arlington.gpkg")) %>%
    st_transform(st_crs(32616)) %>%
    mutate(geom = st_centroid(geom))

# [1] "2013-07-01" "2014-07-01" "2015-06-16" "2016-06-16" "2018-07-26" "2019-07-01"

# 1 2015  2015-06-12    10      5      3      1      1
# 2 2016  2016-06-08     6      4      1      0      1
# 3 2016  2016-06-14    10      3      0      0      7

obs_par_df %>%
    # filter(year %in% c(2016),
    #        obs_day > yday(ymd("2022-05-30"))-1, obs_day < yday(ymd("2022-07-15"))-1) %>%
    filter(year == 2013) %>%
    add_rr_rs_fct() %>%
    group_by(year, obs_date, obs_day) %>%
    summarize(obs_n = n(),
              rr_rs0 = sum(rr_rs_fct == levels(rr_rs_fct)[4]),
              rr_rs1 = sum(rr_rs_fct == levels(rr_rs_fct)[3]),
              rr_rs2 = sum(rr_rs_fct == levels(rr_rs_fct)[2]),
              rr_rs3 = sum(rr_rs_fct == levels(rr_rs_fct)[1]),
              .groups = "drop") %>%
    # filter(rr_rs0 > 0, (rr_rs2 + rr_rs3) > 0) %>%
    # filter(obs_n == max(obs_n)) %>%
    print(n = 50)




# obs_par_df %>%
#     filter(year == 2013) %>%
#     filter(obs_day %in% c(174, 177)) %>%
#     # .[["field"]] %>% unique() %>% length() %>%
#     # group_by(field) %>% summarize(no = n(), min = min(para), max = max(para)) %>%
#     identity()
#
# # Looking at pairs of dates:
# fdf <- obs_par_df %>%
#     filter(year == 2013) %>%
#     filter(!obs_day %in% c(181, 223, 230)) %>%
#     group_by(obs_date, obs_day) %>%
#     summarize(fields = list(field), .groups = "drop")
#
#
# combn(nrow(fdf), 2) %>%
#     t() %>%
#     as.data.frame() %>%
#     setNames(c("d1", "d2")) %>%
#     as_tibble() %>%
#     mutate(n_unq = map2_int(d1, d2, function(.x, .y) {
#         c(fdf$fields[[.x]], fdf$fields[[.y]]) %>%
#             unique() %>%
#             length()
#     }),
#     d1 = map_chr(d1, ~ paste(fdf$obs_date)[.x]) %>% as.Date(),
#     d2 = map_chr(d2, ~ paste(fdf$obs_date)[.x]) %>% as.Date()) %>%
#     filter((d2 - d1) < 7) %>%
#     filter(n_unq > 5) %>%
#     print(n = 50)
#
#
# obs_par_df %>%
#     filter(obs_date == as.Date(c("2014-06-06"))) %>%
#     .[["field"]] %>%
#     {paste0('"', ., '"', collapse = ", ")} %>%
#     cat("\n")




obs_fields_par <- obs_par_df %>%
    filter(obs_date %in% maps_dates) %>%
    # Make sure all fields are present in all dates within year:
    group_by(year, field) %>%
    mutate(no = n()) %>%
    ungroup() %>%
    filter(no == 3) %>%
    mutate(obs = obs %>% fct_drop()) %>%
    split(1:nrow(.)) %>%
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf %>% filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$date <- .d$date
        .f$para <- .d$para
        .f$para_n <- .d$para_n
        .f$rr_rs <- .d$rr_rs
        .f$obs <- .d$obs
        .f$obs_date <- .d$obs_date
        return(.f)
    }) %>%
    do.call(what = rbind) %>%
    mutate(plot_date = factor(paste(obs_date),
                              levels = paste(sort(unique(obs_date))),
                              labels = format(sort(unique(obs_date)),
                                              "%Y\n%e %b")))


xy_lims <- st_bbox(obs_fields_par) %>% as.list()
xy_lims$xmin <- xy_lims$xmin - 400
xy_lims$xmax <- xy_lims$xmax + 400
xy_lims$ymin <- xy_lims$ymin - 400
xy_lims$ymax <- xy_lims$ymax + 400

fields_par_p <- obs_fields_par %>%
    # add_thresh() %>%
    add_rr_rs_fct() %>%
    ggplot() +
    geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
              ymin = xy_lims$ymin, ymax = xy_lims$ymax,
              fill = NA, color = "black", size = 0.5) +
    # geom_sf(aes(size = para, color = para), shape = 16) +
    geom_sf(aes(size = para, color = rr_rs_fct), shape = 16) +
    # scale_fill_viridis_c("Parasitism", option = "inferno",
    #                      limits = c(0, 0.85), begin = 0.1, end = 0.9,
    #                      aesthetics = c("color", "fill"),
    #                      breaks = 0.2 * 0:4) +
    # scale_fill_manual(values = thresh_pal,
    #                   aesthetics = c("color", "fill"),
    #                   guide = "none") +
    scale_color_manual(paste("'Relative fitness\nfor resistance\n('",
                             "* r[r] / r[s] * ')'") %>%
                           str2expression(),
                       guide = "none",
                       values = rr_rs_pal) +
    scale_size("Parasitism", limits = c(0, 0.85), range = c(0.5, 8),
               breaks = 0.2 * 0:4) +
    guides(# color = guide_legend(override.aes = list(alpha = 1, size = 3)),
           size = guide_legend()) +
    coord_sf(datum = st_crs(32616),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")])) +
    facet_wrap(~ plot_date, nrow = 2) +
    theme_void() +
    theme(plot.title = element_text(size = 8, margin = margin(0,0,0,t=3))) +
    theme(strip.text = element_text(size = 8, margin = margin(0,0,0,t=3))) +
    scalebar(data = obs_fields_par,
             dist = 1, dist_unit = "km", transform = FALSE,
             border.size = 0.5, st.dist = 0.06,
             st.size = 3, height = 0.025, location = "topleft",
             x.min = xy_lims$xmin + 400, x.max = xy_lims$xmax,
             y.min = xy_lims$ymin, y.max = xy_lims$ymax - 600,
             facet.var = "plot_date",
             facet.lev = levels(obs_fields_par$plot_date)[1]) +
    NULL





mosaic_p <- par_ts_p + fields_par_p +
    plot_annotation(tag_levels = "a") +
    plot_layout(nrow = 2, heights = c(1, 1.5)) &
    theme(plot.tag = element_text(size = 14, face = "bold"))


{
    cairo_pdf("~/Desktop/mosaic.pdf", width = 7, height = 7)
    plot(mosaic_p)
    dev.off()
}









wi_bounds <- paste0("~/Box Sync/eco-evo_experiments/field-data/",
                    "WI-boundary/Wisconsin_State_Boundary_24K.gpkg") %>%
    st_read() %>%
    st_transform(st_crs(32616))
wi_xy_lims <- st_bbox(wi_bounds) %>% as.list()


wi_inset <- wi_bounds %>%
    ggplot() +
    geom_sf(size = 0.25) +
    # geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
    #           ymin = xy_lims$ymin, ymax = xy_lims$ymax,
    #           fill = "black", color = NA) +
    geom_point(x = (xy_lims$xmin + xy_lims$xmax) / 2,
               y = (xy_lims$ymin + xy_lims$ymax) / 2,
               color = "black", shape = 20, size  = 4) +
    coord_sf(datum = st_crs(32616)) +
    # scalebar(dist = 100, dist_unit = "km", transform = FALSE,
    #          border.size = 0.5, st.dist = 0.05,
    #          st.size = 3, height = 0.02, location = "topleft",
    #          x.min = wi_xy_lims$xmin, x.max = wi_xy_lims$xmax,
    #          y.min = wi_xy_lims$ymin, y.max = wi_xy_lims$ymax + 50e3) +
    north(location = "bottomleft", symbol = 10, scale = 0.2,
          x.min = wi_xy_lims$xmin, x.max = wi_xy_lims$xmax,
          y.min = wi_xy_lims$ymin, y.max = wi_xy_lims$ymax) +
    theme_void() +
    NULL

# wi_inset

# ggsave("~/Desktop/wi_inset.pdf", wi_inset, width = 3, height = 3)


fig1bc <- fields_par_p

# wi_inset + fields_par_p +
#     plot_layout(design = c(area(1, 1, 2, 2),
#                            area(1, 2, 10, 10)))



fig1a + (fig1bc) +
    plot_layout(nrow = 1, widths = c(1, 1.5)) +
    plot_annotation(tag_levels = list(c("A", "B")))







