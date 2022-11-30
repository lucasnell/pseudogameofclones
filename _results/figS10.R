
library(gameofclones)
library(scales)
library(grid)

# colors for resistant, susceptible, and parasitoid wasps, respectively
col_pal <- list(r = "#CCCC00", s = "blue", w = "red")

# File name for figure:
file_out <- "_results/plots/fig_S10.pdf"


rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
}



#'
#' Define aphid line information.
#' Both lines start with 32 adult aphids.
#'

# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Resistant line: high resistance, low parasitized-aphid survival rate,
#                 low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")

# shared end ----


delta.list <- exp(.5*(-10:10))

# baseline
wasp.disp.list <- sort(c(.002*(0:25), .001,.003))

disp <- .2
a <- 2.32
#a <- 2.32/4

n_fields <- 4
save_every <- 1
wasp_density_0 <- .01*rep(1,n_fields)
wasp_density_0[1] <- 0

harvesting.length <- 28
day.interval <- harvesting.length/n_fields
# max_t needs to be a multiple of the harvesting length
max_t <- harvesting.length * 3000
tol <- 1e-3

extinct_N <- 1e-10

s <- 0.025
min.surv <- s
max.surv <- s
n.events <- round(max_t/day.interval)

perturb <- data.frame(when=rep(day.interval*(1:n.events), each=3), where=rep(1:n_fields, each=3), who=c("resistant","susceptible","mummies"), how=runif(3*n.events*n_fields, min=min.surv, max=max.surv))
# kill all mummies
perturb$how[perturb$who == "mummies"] <- 0

perturb <- perturb[1:(3*n.events),]

pick.wasp.disp.list <- 16:17

cairo_pdf(file_out, width = 7, height = 7)
	for(pick.wasp.disp in pick.wasp.disp.list){

		wasp.disp <- wasp.disp.list[pick.wasp.disp]
		if(pick.wasp.disp == pick.wasp.disp.list[1]) {
			par(mfrow=c(2,1), mai=c(.1,1,.8,.1))
			xaxt <- "n"
			plab <- "A"
			plab_y <- 0.9
		} else {
			par(mai=c(.8,1,.1,.1))
			xaxt <- NULL
			plab <- "B"
			plab_y <- 0.5
		}

		sim <- sim_experiments(clonal_lines = c(line_s, line_r),
								n_fields = n_fields,
								wasp_density_0 = wasp_density_0,
								alate_field_disp_p = disp,
								wasp_disp_p = wasp.disp,
								perturb = perturb,
								a = a,
								K = 12500,
								#h = 0,
								#extinct_N = 1e-10,
								extinct_N = extinct_N,
		                        max_t = max_t, save_every = save_every)
		sim <- rm_tibs(sim)

		w <- sim$aphids
		ww <- sim$wasps

		count <- 16
		w <- w[w$time > max_t - count * harvesting.length,]
		ww <- ww[ww$time > max_t - count * harvesting.length,]


		w$time0 <- w$time - min(w$time)
		ww$time0 <- ww$time - min(ww$time)

		.lty = 1
		.lwd = 1.5
		.alpha = 0.25

		i.field <- 2
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col=col_pal$r, log="y", xlab="Days", ylim = c(.00001,10000), xaxt = xaxt, yaxt = "n", lwd = 2, ylab = "Abundance")
		aty <- axTicks(2)
		y_labels <- sapply(aty,function(i) {
		    li <- log10(i)
		    if (li %% 3 == 0) return(as.expression(bquote(10^ .(li))))
		    expression("")
		})
		axis(2,at=aty,labels=y_labels, las=1)
		grid.text(plab, x = unit(0, "npc"), y = unit(plab_y, "npc"),
		          just = c("left", "top"),
		          gp = gpar(fontsize = 18, fontface = "bold"))
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col=col_pal$s, lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col=col_pal$w, lwd = 2)
		for(i.field in c(1,3,4)){
			lines(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], col=alpha(col_pal$r, .alpha), lwd = .lwd, lty=.lty)
			lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col=alpha(col_pal$s, .alpha), lwd = .lwd, lty=.lty)
			lines(wasps ~ time0, data=ww[ww$field == i.field,], col=alpha(col_pal$w, .alpha), lwd = .lwd, lty=.lty)
		}
	}
dev.off()



#'
#' If pdfcrop and ghostscript are installed, use `pdfcrop` to trim whitespace.
#' This isn't necessary to replicate the main plot, but it helps in
#' incorporating the plot into the final document.
#'
if (nzchar(Sys.which("pdfcrop")) && nzchar(tools::find_gs_cmd())) {
    system2("pdfcrop", shQuote(c(file_out, file_out)), stdout = FALSE)
}
