
library(tidyverse)
library(ggtext)
library(lubridate)
library(gameofclones)
library(viridisLite)
library(patchwork)
library(sf)
library(s2)
library(transformr)
library(ggsn)
library(parallel)
library(here)

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

    # File of saved r_r / r_s lookup table (may not exist):
    rr_rs_rds <- here("_results/_data/rr_rs_lookup.rds")

    if (file.exists(rr_rs_rds)) {

        rr_rs_df <- readRDS(rr_rs_rds)

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
        saveRDS(rr_rs_df, rr_rs_rds)
    }

    # Now we look up values of relative fitnesses based on the attack rate
    # associated with the parasitism reported (and the assumed proportion of
    # aphids that are currently resistant)
    para_lookup <- (1 - p_res) * rr_rs_df$ps + p_res * rr_rs_df$pr
    rr_rs <- sapply(para, function(p) {
        nearest <- which(abs(para_lookup - p) < 1.5e-5)

        # dealing with known discontinuity when p_res = 0.48 (default):
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





# Time series plot ----

#'
#' This essentially replicates (and updates) Nature E&E paper.
#'
#' I'm not using lines here bc I think points illustrate my point better,
#' and because code I used to separate by cycle for years 2011--2016 doesn't
#' appear to work for 2017--2019.
#'



# This is used for both time series and map
rr_rs_legend_title <- str_c("Relative fitness for<br>resistant aphids<br>",
                          "(*r<sub>r</sub>* / *r<sub>s</sub>*)")


par_ts_p <- par_df |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer())) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) |>
    add_rr_rs_fct() |>
    ggplot(aes(plot_date, para)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
    geom_segment(data = tibble(plot_date = yday(maps_dates) |>
                                   as.Date(origin = "2021-12-31"),
                               year = year(maps_dates) |> factor(),
                               para = -0.15),
               aes(xend = plot_date, yend = -0.05),
               linewidth = 0.5, linejoin = "mitre",
               arrow = arrow(length = unit(0.1, "lines"), type = "closed")) +
    geom_point(aes(color = rr_rs_fct, fill = rr_rs_fct),
               size = 1, shape = 21) +
    facet_wrap(~ year, ncol = 3) +
    scale_color_manual(rr_rs_legend_title,
                      values = rr_rs_pal$color) +
    scale_fill_manual(rr_rs_legend_title,
                      values = rr_rs_pal$fill2) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous("Parasitism", breaks = 0.4*0:2) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          legend.title = element_markdown(hjust = 0),
          strip.text = element_text(size = 9)) +
    NULL

# par_ts_p








# splitting into periods ----
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


obs_par_df |>
    filter(obs_n >= 5) |>
    group_by(obs_date) |>
    summarize(para = sd(para)) |>
    summarize(median = median(para),
              mean = mean(para))


obs_par_df |>
    group_by(field, year) |>
    summarize(para = sd(para), .groups = "drop") |>
    summarize(median = median(para),
              mean = mean(para))


# maps ----


fields_sf <- st_read(here("_results/_data/Arlington.geojson")) |>
    st_transform(st_crs(32616)) |>
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
           plot_date_by_yr = factor(plot_date_by_yr))


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
           year = "2013", plot_date_by_yr = "1")
fields_par_lab_df <- obs_fields_par |>
    distinct(year, plot_date_by_yr, plot_date) |>
    mutate(x = xy_lims$xmin + 100,
           x = (xy_lims$xmin + xy_lims$xmax) / 2,
           y = xy_lims$ymax - 100,
           y = xy_lims$ymax + 100)

fields_par_p <- obs_fields_par |>
    add_rr_rs_fct() |>
    ggplot() +
    geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
              ymin = xy_lims$ymin, ymax = xy_lims$ymax,
              fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(aes(size = para, color = rr_rs_fct, fill = rr_rs_fct), shape = 21) +
    # geom_text(data = fields_par_lab_df, aes(x, y, label = plot_date),
    #           size = 9 / 2.83465, hjust = 0.5, vjust = 0) +
    # # ------*
    # # DIY scale bar:
    # geom_rect(data = fields_par_scale_df,
    #              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #              fill = "black", color = NA) +
    # geom_text(data = fields_par_scale_df |>
    #               mutate(x = (xmin + xmax) / 2, y = ymin - 200),
    #           aes(x = x, y = y), label = "1 km",
    #           size = 9 / 2.83465, vjust = 1) +
    # ------*
    scale_color_manual(rr_rs_legend_title, guide = "none",
                       values = rr_rs_pal$color) +
    scale_fill_manual(rr_rs_legend_title, guide = "none",
                      values = rr_rs_pal$fill) +
    scale_size("Parasitism", limits = c(0, 0.85), range = c(0.5, 8),
               breaks = 0.2 * 0:4) +
    guides(size = guide_legend(override.aes = list(shape = 16))) +
    coord_sf(datum = st_crs(32616),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")]),
             clip = "off") +
    facet_wrap(~ plot_date, nrow = 2) +
    # facet_grid(year ~ plot_date_by_yr, switch = "y") +
    theme_void() +
    theme(strip.text = element_text(size = 9, margin = margin(0,0,b=3,t=3))) +
    # theme(strip.text.x = element_blank(),
    #       strip.text.y.left = element_text(size = 12, angle = 0,
    #                                        margin = margin(0,0,0,r=6)),
    #       panel.spacing.y = unit(9, "pt")) +
    #       # plot.margin = margin(0,0,0,t=9)) +
    scalebar(data = obs_fields_par,
             dist = 1, dist_unit = "km", transform = FALSE,
             border.size = 0.5, st.dist = 0.06,
             st.size = 3, height = 0.025, location = "topleft",
             x.min = xy_lims$xmin + 400, x.max = xy_lims$xmax,
             y.min = xy_lims$ymin, y.max = xy_lims$ymax - 600,
             facet.var = "plot_date",
             facet.lev = levels(obs_fields_par$plot_date)[1]) +
             # facet.var = c("plot_date_by_yr", "year"),
             # facet.lev = c("1", "2013")) +
    NULL







mosaic_p <- function() {
    p <- par_ts_p + fields_par_p +
        plot_annotation(tag_levels = "A") +
        plot_layout(nrow = 2, heights = c(1, 1.5)) &
        theme(plot.tag = element_text(size = 14, face = "bold"))
    plot(p)
    grid.text("2013", x = unit(0.08, "npc"), y = unit(0.385, "npc"),
              just = c("right", "center"))
    grid.text("2015", x = unit(0.08, "npc"), y = unit(0.135, "npc"),
              just = c("right", "center"))
}

# mosaic_p()



# save_plot(here("_results/_plots/mosaic.pdf"), mosaic_p, 7, 7)





# --------------------------------------------------------*
# Hamiltonella frequencies ----
# These will be listed in text.
# --------------------------------------------------------*

#'
#' Frequencies of Hamiltonella by field and time
#'
ham_df <- here("_results/_data/symbionts-2012-2017.csv") |>
    read_csv(col_types = cols()) |>
    select(year, late, date, field, Hamiltonella) |>
    group_by(year, late, date, field) |>
    summarize(ham = mean(Hamiltonella), n = n(), .groups = "drop")

#' Average number of aphids assayed:
ham_df$n |> mean()

#' Some summaries of Hamiltonella frequencies:
ham_df$ham |> mean()
ham_df$ham |> sd()
ham_df$ham |> range()

