
# Shared packages
suppressPackageStartupMessages({
    library(tidyverse)
    library(pseudogameofclones)
    library(future.apply)
    library(progressr)
    library(patchwork)
    library(viridisLite)
    library(ggtext)
})

# number of threads (used not just for futures):
.n_threads <- max(1L, parallel::detectCores() - 2L)

plan(multisession, workers = .n_threads)
handlers(global = TRUE)
handlers("progress")



#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


# color palette for plant types:
type_pal <- c(viridis(3, begin = 0.1, end = 0.9), "gray80")[c(2,1,4,3)] |>
    set_names(c("virus", "pseudo", "none", "both"))


