
library(tidyverse)
library(gameofclones)
library(ggtext)         # element_markdown
library(lubridate)      # yday, year, ymd
library(viridisLite)    # inferno
library(readxl)         # readxl
library(here)           # here
library(grid)           # grid.newpage, grid.draw



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
# Shared objects ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------


# Color palettes to use in plots (using binned rr_rs values):
rr_rs_pal <- with(list(pal = inferno, inds = c(80, 60, 20)),
                  list(color = c(pal(100)[inds], "gray40"),
                       fill = c(pal(100)[inds], "white"),
                       fill2 = c(pal(100, alpha = 0.5)[inds], "white")))

# Add factor that breaks rr_rs into bins for plotting:
add_rr_rs_fct <- function(.df) {
    .df |>
        mutate(rr_rs_fct = cut(rr_rs, c(0, 1, 1.05, 1.1, Inf),
                               labels = c("< 1", "1 – 1.05", "1.05 – 1.1", "> 1.1")),
               rr_rs_fct = factor(rr_rs_fct, levels = rev(levels(rr_rs_fct))))
}

#' Function to create time series plots and keep axes and font sizes consistent.
ts_p_maker <- function(.df, .y, .ylab, ...) {
    ggplot(.df, aes(plot_date, {{ .y }})) +
        geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
        list(...) +
        facet_wrap(~ year, ncol = 3, drop = FALSE) +
        scale_x_date(breaks = ymd(sprintf("2022-%02i-01", 5:11)),
                     labels = c("","Jun","","Aug", "","Oct",""),
                     limits = as.Date(c("2022-04-27", "2022-10-30"))) +
        scale_y_continuous(.ylab, limits = c(-0.2, 0.9), breaks = 0.4*0:2) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(color = "black", size = 8),
              axis.title.y = element_markdown(hjust = 0.5),
              legend.title = element_markdown(hjust = 0),
              legend.position = "none",
              strip.text = element_text(size = 9))
}

#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))





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




#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Hamiltonella dataset - read and organize ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------



older_ham_df <- here("_results/_data/symbionts-2018-2019.csv") |>
    read_csv(col_types = cols()) |>
    base::`[`()



newer_ham_df <- here("_results/_data/symbionts-2012-2017.csv") |>
    read_csv(col_types = cols()) |>
    mutate(season = case_when(is.na(date) ~ "fall",
                              late == 1 ~ "fall",
                              late == 0 ~ "spring")) |>
    select(year, season, date, field, clone, Hamiltonella) |>
    rename(ham = Hamiltonella) |>
    # Filling in this date manually based on info from Kerry Oliver:
    mutate(date = ifelse(year == 2012, "9/4/12", date),
           date = as.Date(date, format = "%m/%d/%y"))


ham_df <- bind_rows(newer_ham_df, older_ham_df) |>
    group_by(year, season, date, field) |>
    summarize(ham = mean(ham), n = n(), .groups = "drop") |>
    #' filter for fields where at least 10 aphids were surveyed:
    filter(n >= 10) |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer())) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(yday(date), origin = "2022-01-01"),
           field_id = interaction(field, year, drop = TRUE),
           year = factor(year, levels = 2011:2019))


#' If you want just some basic stats:
ham_df |>
    summarize(n = mean(n), mean = mean(ham), sd = sd(ham),
              min = min(ham), max = max(ham))




#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Make plots ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' The parasitism time series plot essentially replicates (and updates)
#' Ives et al. (2020, Nature E&E) paper.
#'
#' I'm not using lines here bc I think points illustrate my point better,
#' and because code I used to separate by cycle for years 2011--2016 doesn't
#' appear to work for 2017--2019.
#'


rr_rs_legend_title <- str_c("Relative fitness for<br>resistant aphids<br>",
                            "(*r<sub>r</sub>* / *r<sub>s</sub>*)")


par_ts_p <- par_df |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer())) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) |>
    add_rr_rs_fct() |>
    ts_p_maker(para, "Proportion parasitized",
               geom_segment(data = tibble(plot_date = yday(maps_dates) |>
                                              as.Date(origin = "2021-12-31"),
                                          year = year(maps_dates) |> factor(),
                                          para = -0.15),
                            aes(xend = plot_date, yend = -0.05),
                            linewidth = 0.5, linejoin = "mitre",
                            arrow = arrow(length = unit(0.1, "lines"), type = "closed")),
               geom_point(aes(color = rr_rs_fct, fill = rr_rs_fct),
                          size = 1, shape = 21, stroke = 0.5)) +
    scale_color_manual(rr_rs_legend_title,
                       values = rr_rs_pal$color) +
    scale_fill_manual(rr_rs_legend_title,
                      values = rr_rs_pal$fill2) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

# par_ts_p

save_plot(here("_results/_plots/field-data/par-time.pdf"), par_ts_p,
          w = 4, h = 2.5)

par_ts_p_leg <- function() {
    legend <- (par_ts_p + theme(legend.position = "right")) |>
        (function(a.gplot){
            tmp <- ggplot_gtable(ggplot_build(a.gplot))
            leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
            legend <- tmp$grobs[[leg]]
            legend
        })()
    grid.newpage()
    grid.draw(legend)
}

save_plot(here("_results/_plots/field-data/par-time-legend.pdf"), par_ts_p_leg,
          w = 2, h = 2.5)



#' Same thing but without coloring by fitness and with lines connecting fields.
#' The stuff with `max_fc` is to keep the colors from all being similar when
#' not many fields were sampled in a given year
par_ts_wlines_p <- par_df |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer(),
                     max_fc = max(field_col))) |>
    mutate(max_max_fc = max(max_fc)) |>
    split(~ year) |>
    map_dfr(\(x) {
        mfc <- x[["max_fc"]][[1]]
        mmfc <- x[["max_max_fc"]][[1]]
        x |>
            mutate(field_col = map_int(field_col,
                                        \(i) {
                                            ss <- seq(1, mmfc, length.out = mfc)
                                            as.integer(round(ss))[i]
                                        }))
    }) |>
    select(-ends_with("max_fc")) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) |>
    add_rr_rs_fct() |>
    ts_p_maker(para, "Proportion parasitized",
               geom_line(aes(color = field_col), linewidth = 0.5),
               geom_point(aes(color = field_col, fill = field_col),
                          size = 0.5, shape = 21)) +
    scale_color_viridis_d(guide = "none") +
    scale_fill_viridis_d(guide = "none")

save_plot(here("_results/_plots/par-time-by-field.pdf"), par_ts_wlines_p,
          w = 6.5, h = 5)





ham_ts_p <- ham_df |>
    ts_p_maker(ham, "Proportion *H. defensa*",
               geom_line(aes(group = field_id), color = "gray60"),
               geom_jitter(size = 1, shape = 1, width = 3, height = 0))

# ham_ts_p

save_plot(here("_results/_plots/field-data/ham-time.pdf"), ham_ts_p,
          w = 4, h = 2.5, seed = 380247925)

