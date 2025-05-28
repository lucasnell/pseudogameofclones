
#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}
