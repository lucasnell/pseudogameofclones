
library(numDeriv)
library(lattice)
library(plot3D)
#' Install `clonewars` package from GitHub if necessary
if (!require("gameofclones")) {
    if (!require("remotes")) install.packages("remotes")
    remotes::install_github("lucasnell/gameofclones", force = TRUE)
    library(gameofclones)
}



# cairo_pdf(filename = sprintf("_results/_plots/%s", fn), width = w, height = h, ...)
# dev.off()


# Used to construct the jacobian
comm <- function(x, sim, new_starts = NULL, max_t = 1){
	if(is.null(new_starts)) new_starts <- sim$all_info[[1]]
	new_starts$N <- x
	restart.sim <- restart_experiment(sim, new_starts = new_starts,
	                                  max_t = max_t)$all_info[[1]]
	return(restart.sim$N)
}

conf_bounds <- function(x, y.lower, y.upper, col="lightgray"){
		polygon(c(x, rev(x)), c(y.upper, rev(y.lower)), col=col, border=NA)
}


clone_traj <- function(sim, delta, category = "resistant", max_t = 400,
                       perturb = NULL){

	new.starts <- sim$all_info[[1]]
	pick <- !is.na(new.starts$line) & (new.starts$line == category)
	new.starts[pick, "N"] <- delta * new.starts[pick, "N"]

	new.sim <- restart_experiment(sim, new_starts = new.starts, max_t = max_t,
	                              perturb = perturb)
	new.sim <- rm_tibs(new.sim)
	d <- new.sim$aphids
	d$line.type <- paste0(d$line,".", d$type)
	d$line.type[d$line.type == "NA.mummy"] <- "mummy"
	d <- d[,-c(1,5,6)]
	d <- spread(d, "line.type", "N")
	d <- cbind(delta = delta, d, wasps = new.sim$wasps[,4])
	return(d)
}

clone_converge <- function(sim, delta, category = "resistant", max_t = 1e4,
                           tol = 1e-10, perturb = NULL){

	new.starts <- sim$all_info[[1]]
	pick <- !is.na(new.starts$line) & (new.starts$line == category)
	new.starts[pick, "N"] <- delta * new.starts[pick, "N"]

	new.sim <- restart_experiment(sim, new_starts = new.starts, max_t = max_t,
	                              perturb = perturb)
	ns.df <- new.sim$all_info[[1]]

	S.resistant <- sum(ns.df$N[!is.na(ns.df$line) & ns.df$line == "resistant"])
	S.susceptible <- sum(ns.df$N[!is.na(ns.df$line) &
	                                 ns.df$line == "susceptible"])
	test <- (S.resistant/(S.resistant + S.susceptible) > tol)
	return(test)
}

clone_wasp_converge <- function(sim, delta, max_t = 1e4, tol = 1e-8,
                                perturb = NULL, harvesting.length,
                                day.interval){

	new.starts <- sim$all_info[[1]]
	pick <- !is.na(new.starts$line) & (new.starts$line == "resistant")
	new.starts[pick, "N"] <- delta * new.starts[pick, "N"]

	new.sim <- restart_experiment(sim, new_starts = new.starts,
	                              max_t = max_t, perturb = perturb)
	d <- new.sim$wasps
	d <- d[d$time >= .9 * max_t,]
	d <- d[d$time %% day.interval == 0,]

	d <- aggregate(wasps ~ time, data=d, FUN = sum)
	test.wasp <- max(d$wasps) > tol

	d <- new.sim$aphids
	d <- d[d$time >= .9 * max_t,]
	d <- d[d$time %% day.interval == 0,]

	d <- aggregate(N ~ time + line, data=d, FUN = sum)
	test.resistant <- max(d$N[d$line == "resistant"] > tol)
	test.susceptible <- max(d$N[d$line == "susceptible"] > tol)

	return(all(test.wasp, test.resistant, test.susceptible))
}

clone_plot <- function(sim, delta.list, labels = "", max_t = 1e3,
                       max_t.converge = 1e4, perturb = NULL){
	d <- sim$aphids
    # d <- d[d$time <= max_t,]
	d$line.type <- paste0(d$line,".", d$type)
	d$line.type[d$line.type == "NA.mummy"] <- "mummy"
	d <- d[,-c(1,5,6)]
	#d <- pivot_wider(d, names_from = line.type, values_from = N)
	d <- spread(d, "line.type", "N")
	d <- cbind(delta = 0, d, wasps = sim$wasps[,4], converge = TRUE)

	d <- d[d$time <= max_t,]

	for(i.delta in delta.list){
		dd <- cbind(clone_traj(sim = sim, delta = i.delta,
		                       category = "resistant", max_t = max_t,
		                       perturb = perturb),
		            converge = clone_converge(sim = sim, delta = i.delta,
		                                      category = "resistant",
		                                      max_t = max_t.converge,
		                                      perturb = perturb))
		d <- rbind(d, dd)
	}

	d$resistant <- d$resistant.alate + d$resistant.apterous +
	    d$resistant.parasitized
	d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
	    d$susceptible.parasitized
	d$prop <- d$resistant/(d$resistant + d$susceptible)

	d$col <- "gray"
	d$col[d$delta == 1] <- "black"

	par(mfrow=c(3,1), mai = c(.6,.7,.1,.1))
	plot(prop ~ time, data=d[d$delta == 1 & d$field == 1,], typ="l",
	     ylim=c(0,1), ylab = "Resistant (field 1)", xlab = "Time")
	for (i.delta in delta.list) {
	    lines(prop ~ time, data=d[d$delta == i.delta & d$field == 1,],
	          col=col, lty=2-as.numeric(converge))
    }
	text(x=max_t/2, y=.9, labels=labels, cex=2)
	lines(prop ~ time, data=d[d$delta == 1 & d$field == 1,])

	plot(prop ~ time, data=d[d$delta == 1 & d$field == 2,], typ="l",
	     ylim=c(0,1), ylab = "Resistant (field 2)", xlab = "Time")
	for (i.delta in delta.list) {
	    lines(prop ~ time, data=d[d$delta == i.delta & d$field == 2,],
	          col=col, lty=2-as.numeric(converge))
    }
	lines(prop ~ time, data=d[d$delta == 1 & d$field == 2,])

	plot(wasps ~ time, data=d[d$delta == 1 & d$field == 1,], typ="l",
	     ylim = c(0, max(d$wasps)), ylab = "Wasps (field 1)", xlab = "Time")
	for(i.delta in delta.list) {
	    lines(wasps ~ time, data=d[d$delta == i.delta & d$field == 1,],
	          col=col, lty=2-as.numeric(converge))
    }
	lines(wasps ~ time, data=d[d$delta == 1 & d$field == 1,])
}

clone_eig <- function(sim){
	x <- sim$all_info[[1]]$N
	C <- jacobian(comm, x, sim=sim)
	eig <- eigen(C)
	value <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
	magnitude <- Re(abs(eig$value[abs(eig$values) == max(abs(eig$values))][1]))
	return(list(value = value, magnitude = magnitude))
}

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

    ## Get RGB values for named color
    rgb.val <- col2rgb(color)

    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)

    ## Save the color
    invisible(t.col)
}
## END



#'
#' WARNING: Some of the code below contains tibbles.
#' To remove them, you can use this function on the output
#' from `sim_experiments` or `restart_experiment`:
#'
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
#'
#' #'
#' #' The resulting `aphid` objects have the following fields:
#' #'
#' #' - `name`       : string indicate the name of the line ("susceptible" or
#' #'                  "resistant")
#' #' - `density_0`  : 29x2 matrix with starting abundances of each stage,
#' #'                  for apterous aphids (col 1) and alates (col 2)
#' #' - `attack_surv`: length-2 vector of wasp attack survivals for singly and
#' #'                  multiply parasitized aphids
#' #' - `leslie`     : 29x29x3 array with the 3 leslie matrices for apterous
#' #'                  aphids, alates, and parasitized aphids, respectively.
#' #'                  The latter is largely ignored as in the original model,
#' #'                  but is included for compatibility with the C++ code.
#' #'
#'
#'
#' #'
#' #' The defaults for the `sim_experiments` function are designed to
#' #' approximately replicate the experimental results.
#' #' The two things you'll definitely need to provide it are the clonal lines
#' #' we've defined and the maximum time to run the simulations.
#' #' The default for the latter is 250 days, so that won't work for reaching
#' #' equilibrium.
#' #'
#' #' This takes ~ 5 sec on my computer:
#' #'
#' alate_field_disp_p <- 0.2
#' sims <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = alate_field_disp_p,
#' 						extinct_N = 1e-6,
#'                         max_t = 1e5, save_every = 1e3)
#' sims <- rm_tibs(sims)
#' # sims
#'
#' # par(mfrow=c(2,1))
#' # plot(N ~ time, data=sims$aphids[sims$aphids$field == 1 & sims$aphids$line == "resistant" & sims$aphids$type == "apterous",], typ="l", col="blue", log="y", ylim=c(1,1500))
#' # lines(N ~ time, data=sims$aphids[sims$aphids$field == 1 & sims$aphids$line == "susceptible" & sims$aphids$type == "apterous",])
#' # lines(wasps ~ time, data=sims$wasps[sims$wasps$field == 1,], col="red")
#' #
#' # plot(N ~ time, data=sims$aphids[sims$aphids$field == 2 & sims$aphids$line == "susceptible" & sims$aphids$type == "apterous",], typ="l", log="y", ylim=c(1,1500))
#' # lines(N ~ time, data=sims$aphids[sims$aphids$field == 2 & sims$aphids$line == "resistant" & sims$aphids$type == "apterous",], col="blue")
#' # lines(wasps ~ time, data=sims$wasps[sims$wasps$field == 2,], col="red")
#'
#' ######### strange behavior ----
#' alate_field_disp_p <- 0.23058
#' sims <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = alate_field_disp_p,
#' 						extinct_N = 1e-6,
#'                         max_t = 2*1e3, save_every = 1e1)
#' sims <- rm_tibs(sims)
#' sims
#'
#' par(mfrow=c(2,2))
#' plot(N ~ time, data=sims$aphids[sims$aphids$field == 1 & sims$aphids$line == "resistant" & sims$aphids$type == "apterous",], typ="l", col="blue", log="y", ylim=c(1,1500))
#' lines(N ~ time, data=sims$aphids[sims$aphids$field == 1 & sims$aphids$line == "susceptible" & sims$aphids$type == "apterous",])
#' lines(wasps ~ time, data=sims$wasps[sims$wasps$field == 1,], col="red")
#'
#' plot(N ~ time, data=sims$aphids[sims$aphids$field == 1 & sims$aphids$line == "resistant" & sims$aphids$type == "apterous",], typ="l", col="blue", ylim=c(700,800), xlim=c(100,700))
#'
#' plot(N ~ time, data=sims$aphids[sims$aphids$field == 2 & sims$aphids$line == "susceptible" & sims$aphids$type == "apterous",], typ="l", log="y", ylim=c(1,1500))
#' lines(N ~ time, data=sims$aphids[sims$aphids$field == 2 & sims$aphids$line == "resistant" & sims$aphids$type == "apterous",], col="blue")
#' lines(wasps ~ time, data=sims$wasps[sims$wasps$field == 2,], col="red")
#'
#' plot(N ~ time, data=sims$aphids[sims$aphids$field == 2 & sims$aphids$line == "susceptible" & sims$aphids$type == "apterous",], typ="l", ylim=c(1200,1400), xlim=c(100,700))
#'
#'
#' #'
#' #' The resulting `cloneSims` object contains the following fields
#' #' (the last two don't show up when printing `sims` bc they shouldn't be
#' #'  messed with):
#' #'
#' #' - `aphids`       : A `tibble` of aphid abundances through time and by cage.
#' #'                    These data are not stage structured.
#' #'                    Note that in the simulations, the `field` column is
#' #'                    equivalent to the experimental cage.
#' #' - `wasps`        : A `tibble` of wasp abundances through time and by cage.
#' #' - `all_info`     : A list of length 1 containing a data frame with the
#' #'                    abundances at the end of the simulations.
#' #'                    These data are stage structured and should be used
#' #'                    as the starting point for perturbations if using
#' #'                    the `restart_experiment` function described below.
#' #' - `all_info_xptr`: A pointer to the C++ object used to run the simulations
#' #'                    that is used to restart the simulations for doing
#' #'                    perturbations. This should not be touched.
#' #' - `call`         : A long list containing the call information for the
#' #'                    simulations. This can and should be ignored.
#' #'
#'
#' #'
#' #' To perturb and restart experimental simulations.
#' #'
#' #' Data frame of stage-structured ending abundances:
#' ss_df <- sims$all_info[[1]]
#' # ss_df
#'
#' #' The data frame should have the following columns:
#' #'
#' #' - `field`: This is equivalent to experimental cage.
#' #'            Field 1 has wasps, field 2 doesn't.
#' #' - `plant`: This is a holdover from when I was simulating individual plants
#' #'            withering. It can be ignored.
#' #' - `line` : The name of the associated aphid line. This is an empty string
#' #'            for wasps and mummies.
#' #' - `type` : The organism type. This can be `"wasp"`, `"mummy"`, `"apterous"`,
#' #'            `"alate"`, or `"parasitized"`.
#' #' - `stage`: The stage (in days). There should be 29 stages for alate and
#' #'            apterous aphids, 7 for parasitized aphids, 4 for mummies, and
#' #'            1 for adult wasps.
#' #' - `N`    : Abundance.
#' #'
#'
#' #'
#' #' How to perturb abundances and restart the simulations.
#' #' Three key points for perturbations:
#' #' 1) You shouldn't change the `all_info` field directly from the original
#' #'    simulation object (`sims` in this case).
#' #'    Always copy to a new object before editing.
#' #' 2) You should only change the `N` column and leave everything else the same.
#' #' 3) Only try perturbation on output from `sim_experiments`, not on output
#' #'    from `restart_experiments`. If you need to do multiple perturbations,
#' #'    the `perturb` argument to the `sim_experiments` function is a better
#' #'    option.
#' #'
#' #'
#' #' Here are some examples:
#' #'
#'
#' # Perturb wasps in the wasp cage:
#' delta_ss_df <- ss_df
#' inds <- ss_df$field == 2 & ss_df$type == "wasp"
#' delta_ss_df[inds, "N"] <- ss_df[inds, "N"] * (1 + 1e-5)
#'
#' # Perturb 10-day-old resistant alates in the no-wasp cage:
#' delta_ss_df <- ss_df
#' inds <- ss_df$field == 2 & ss_df$type == "alate" & ss_df$line == "resistant" & ss_df$stage == 10
#' delta_ss_df[inds, "N"] <- ss_df[inds, "N"] * (1 + 1e-5)
#'
#' # Do new simulation for 1 time step:
#' new_sims <- restart_experiment(sims, new_starts = delta_ss_df, max_t = 1)
#' new_sims <- rm_tibs(new_sims)
#' # new_sims
#'
#' # new_sims$aphids
#'
#' # Full stage-structured data from new sims:
#' new_ss_df <- new_sims$all_info[[1]]
#' # new_ss_df
#'
#'
#' ####################################################################*
#' ####################################################################*
#' # plot of stage structure through time ----
#' delta <- 1e-3
#'
#' ss_df <- sims$all_info[[1]]
#' delta_ss_df <- ss_df
#'
#' pick <- delta_ss_df$N > -1
#' pick <- delta_ss_df$line == "resistant" & delta_ss_df$type == "apterous"
#' pick <- delta_ss_df$field == 1 & delta_ss_df$line == "resistant" & delta_ss_df$type == "apterous"
#' # pick <- delta_ss_df$field == 1 & delta_ss_df$line == "resistant"
#' # pick <- delta_ss_df$field == 2 & delta_ss_df$line == "susceptible"
#' # pick <- delta_ss_df$field == 2 & delta_ss_df$type == "wasp"
#'
#' delta_ss_df$N[pick] <- delta_ss_df$N[pick] + delta * delta_ss_df$N[pick]
#' ss_df_29 <- restart_experiment(sims, new_starts = delta_ss_df, max_t = 29)$all_info[[1]]
#'
#' Tmax <- 6*29 + 1
#' X <- array(dim = c(Tmax, dim(delta_ss_df)[1]))
#' X[1,] <- delta_ss_df$N
#'
#' for(t in 2:Tmax){
#' 	sims.restart <- restart_experiment(sims, new_starts = delta_ss_df, max_t = 1)
#' 	delta_ss_df <- sims.restart$all_info[[1]]
#' 	X[t,] <- delta_ss_df$N
#' }
#' par(mfrow=c(4,1))
#' matplot(X[,pick] - matrix(1,Tmax, 1) %*% ss_df$N[pick], typ="l")
#' matplot(X[,pick] / (matrix(1,Tmax, 1) %*% ss_df$N[pick]), typ="l")
#'
#' # pick.short <- pick
#' # pick.short[pick][-c(1,2)] <- FALSE
#' # matplot(X[,pick.short]/(matrix(1,Tmax, 1) %*% ss_df$N[pick.short]), typ="l")
#'
#' plot(rowSums(X[,pick])/sum(ss_df$N[pick]), typ="l")
#' arrows(29,1+delta,29,1)
#' arrows(2*29,1+delta,2*29,1)
#' arrows(3*29,1+delta,3*29,1)
#'
#' # Another way to view stage structure through time.
#' # Just added this feature for the `restart_experiments` function
#' delta_ss_df <- ss_df
#' delta_ss_df$N[pick] <- ss_df$N[pick] + delta * ss_df$N[pick]
#' sims.alt <- restart_experiment(sims, new_starts = delta_ss_df, max_t = Tmax,
#'                              stage_ts_out = TRUE)
#'
#' X <- array(dim = c(Tmax, dim(delta_ss_df)[1]))
#' X[1,] <- delta_ss_df$N
#' for(t in 2:Tmax){
#' 	X[t,] <- sims.alt$stage_ts[[t]]$N
#' }
#' matplot(X[,pick]/(matrix(1,Tmax, 1) %*% ss_df$N[pick]), typ="l")
#'
#'
#' ####################################################################*
#' ####################################################################*
#' # JACOBIAN
#'
#' alate_field_disp_p <- 0.23058
#' max_t <- 1000
#' sim <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = alate_field_disp_p,
#' 						extinct_N = 1e-6,
#'                         max_t = max_t, save_every = 1e0)
#' sim <- rm_tibs(sim)
#' # sim$all_info[[1]]
#'
#'
#' par(mfrow=c(2,1))
#' plot(N ~ time, data=sim$aphids[sim$aphids$field == 1 & sim$aphids$line == "resistant" & sim$aphids$type == "apterous",], typ="l", col="blue", log="y", ylim=c(.1,1500))
#' lines(N ~ time, data=sim$aphids[sim$aphids$field == 1 & sim$aphids$line == "susceptible" & sim$aphids$type == "apterous",])
#' lines(wasps ~ time, data=sim$wasps[sim$wasps$field == 1,], col="red")
#'
#' plot(N ~ time, data=sim$aphids[sim$aphids$field == 2 & sim$aphids$line == "susceptible" & sim$aphids$type == "apterous",], typ="l", log="y", ylim=c(.1,1500))
#' lines(N ~ time, data=sim$aphids[sim$aphids$field == 2 & sim$aphids$line == "resistant" & sim$aphids$type == "apterous",], col="blue")
#' lines(wasps ~ time, data=sim$wasps[sim$wasps$field == 2,], col="red")
#'
#' # demonstration that the Jacobian works ----
#' max_t <- 400
#' sim <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = alate_field_disp_p,
#' 						extinct_N = 1e-6,
#'                         max_t = max_t, save_every = 1e0)
#' x <- sim$all_info[[1]]$N
#' C <- jacobian(comm, x, sim=sim)
#' eig <- eigen(C)
#' sort(abs(eig$values), decreasing=TRUE)
#'
#'
#' ####################################################################*
#' ####################################################################*
#' # aphid dispersal
#'
#' disp.list <- .01*(1:30)
#' max_t <- 400
#' sim <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = disp.list[15],
#' 						extinct_N = 1e-10,
#'                         max_t = max_t, save_every = 1e0)
#' sim <- rm_tibs(sim)
#' # sim$all_info[[1]]
#'
#' par(mfrow=c(2,1))
#' plot(N ~ time, data=sim$aphids[sim$aphids$field == 1 & sim$aphids$line == "resistant" & sim$aphids$type == "apterous",], typ="l", col="blue", log="y", ylim=c(.1,1500))
#' lines(N ~ time, data=sim$aphids[sim$aphids$field == 1 & sim$aphids$line == "susceptible" & sim$aphids$type == "apterous",])
#' lines(wasps ~ time, data=sim$wasps[sim$wasps$field == 1,], col="red")
#'
#' plot(N ~ time, data=sim$aphids[sim$aphids$field == 2 & sim$aphids$line == "susceptible" & sim$aphids$type == "apterous",], typ="l", log="y", ylim=c(.1,1500))
#' lines(N ~ time, data=sim$aphids[sim$aphids$field == 2 & sim$aphids$line == "resistant" & sim$aphids$type == "apterous",], col="blue")
#' lines(wasps ~ time, data=sim$wasps[sim$wasps$field == 2,], col="red")
#'
#' max_t <- 1e4
#' sim <- sim_experiments(clonal_lines = c(line_s, line_r),
#' 						alate_field_disp_p = disp.list[15],
#' 						extinct_N = 1e-10,
#'                         max_t = max_t, save_every = 1e0)
#'
#' x <- sim$all_info[[1]]$N
#' C <- jacobian(comm, x, sim=sim)
#' eig <- eigen(C)
#' eig$value[abs(eig$values) == max(abs(eig$values))][1]
#' abs(eig$value[abs(eig$values) == max(abs(eig$values))][1])
#'
#' ########################*
#' # iterate over disp.list
#' sim.list <- list()
#' sim.list[[15]] <- sim
#' df <- data.frame(count = 1:length(disp.list), var = disp.list, eig = NA)
#'
#' aphids <- sim$aphids[sim$aphid$time == max(sim$aphid$time),]
#' df$aphids.sus.0[df$var==disp.list[15]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "susceptible"])
#' df$aphids.res.0[df$var==disp.list[15]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "resistant"])
#' df$aphids.sus.1[df$var==disp.list[15]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "susceptible"])
#' df$aphids.res.1[df$var==disp.list[15]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "resistant"])
#'
#' wasps <- sim$wasps[sim$wasp$time == max(sim$wasp$time),]
#' df$wasps.0[df$var==disp.list[15]] <- wasps$wasps[1]
#' df$wasps.1[df$var==disp.list[15]] <- wasps$wasps[2]
#'
#' df$eig[df$var==disp.list[15]] <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
#' df$abseig[df$var==disp.list[15]] <- abs(eig$value[abs(eig$values) == max(abs(eig$values))][1])
#'
#' new.starts <- sim$all_info[[1]]
#' for(count in 16:30){
#'
#' 	new.sim <- restart_experiment(sim, new_starts = new.starts, alate_field_disp_p = disp.list[count], max_t = 1e4)
#' 	sim.list[[count]] <- new.sim
#'
#' 	new.starts$N <- new.sim$all_info[[1]]$N
#' 	new.x <- new.sim$all_info[[1]]$N
#' 	# NOTE: this violates the restriction that restart_experiment() isn't started with an object produced by restart_experiment()
#' 	C <- jacobian(comm, new.x, sim=new.sim)
#' 	eig <- eigen(C)
#'
#' 	aphids <- new.sim$aphids[new.sim$aphid$time == max(new.sim$aphid$time),]
#' 	df$aphids.sus.0[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "susceptible"])
#' 	df$aphids.res.0[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "resistant"])
#' 	df$aphids.sus.1[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "susceptible"])
#' 	df$aphids.res.1[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "resistant"])
#'
#' 	wasps <- new.sim$wasps[new.sim$wasp$time == max(new.sim$wasp$time),]
#' 	df$wasps.0[df$var==disp.list[count]] <- wasps$wasps[1]
#' 	df$wasps.1[df$var==disp.list[count]] <- wasps$wasps[2]
#'
#' 	df$eig[df$var==disp.list[count]] <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
#' 	df$abseig[df$var==disp.list[count]] <- abs(eig$value[abs(eig$values) == max(abs(eig$values))][1])
#'
#' 	show(df[count,])
#'
#' 	if(max(abs(eig$values)) > 1) break
#' }
#'
#' new.starts <- sim$all_info[[1]]
#' for(count in 14:1){
#'
#' 	new.sim <- restart_experiment(sim, new_starts = new.starts, alate_field_disp_p = disp.list[count], max_t = 1e4)
#' 	sim.list[[count]] <- new.sim
#'
#' 	new.starts$N <- new.sim$all_info[[1]]$N
#' 	new.x <- new.sim$all_info[[1]]$N
#' 	# NOTE: this violates the restriction that restart_experiment() isn't started with an object produced by restart_experiment()
#' 	C <- jacobian(comm, new.x, sim=new.sim)
#' 	eig <- eigen(C)
#'
#' 	aphids <- new.sim$aphids[new.sim$aphid$time == max(new.sim$aphid$time),]
#' 	df$aphids.sus.0[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "susceptible"])
#' 	df$aphids.res.0[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 1 & aphids$line == "resistant"])
#' 	df$aphids.sus.1[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "susceptible"])
#' 	df$aphids.res.1[df$var==disp.list[count]] <- mean(aphids$N[aphids$field == 2 & aphids$line == "resistant"])
#'
#' 	wasps <- new.sim$wasps[new.sim$wasp$time == max(new.sim$wasp$time),]
#' 	df$wasps.0[df$var==disp.list[count]] <- wasps$wasps[1]
#' 	df$wasps.1[df$var==disp.list[count]] <- wasps$wasps[2]
#'
#' 	df$eig[df$var==disp.list[count]] <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
#' 	df$abseig[df$var==disp.list[count]] <- abs(eig$value[abs(eig$values) == max(abs(eig$values))][1])
#'
#' 	show(df[count,])
#'
#' 	if(max(abs(eig$values)) > 1) break
#' }
#'
#' df <- df[!is.na(df$eig),]
#' par(mfrow=c(2,1))
#' yvals <- c(df$aphids.res.0, df$aphids.sus.0, df$wasps.0)
#' plot(aphids.res.0 ~ var, data=df, typ="l", col="blue", log="y", ylim=c(1+min(yvals),max(yvals)))
#' lines(aphids.sus.0 ~ var, data=df)
#' lines(wasps.0 ~ var, data=df, col="red")
#'
#' yvals <- c(df$aphids.res.1, df$aphids.sus.1, df$wasps.1)
#' plot(aphids.res.1 ~ var, data=df, typ="l", col="blue", log="y", ylim=c(1+min(yvals),max(yvals)))
#' lines(aphids.sus.1 ~ var, data=df)
#' lines(wasps.1 ~ var, data=df, col="red")
#'
#' # proportions
#' par(mfrow=c(2,1))
#' plot(df$var, df$aphids.res.0/(df$aphids.res.0 + df$aphids.sus.0), typ="l", ylim=c(0,1))
#' plot(df$var, df$aphids.res.1/(df$aphids.res.1 + df$aphids.sus.1), typ="l", ylim=c(0,1))
#'
#' saveRDS(sim.list, file="sim.list.RDS")
#'
#' df

####################################################################*
####################################################################*
# domains of attraction

########################*
# 3D plot
# library(plot3D)

max_t <- 1e3
disp_p <- .258
delta.list <- c(50, 40,30, 2, .85, .8)
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
						alate_field_disp_p = .15,
						extinct_N = 1e-10,
                        max_t = max_t, save_every = 1e2)
sim <- restart_experiment(sim,
					alate_field_disp_p = disp_p,
					new_starts = sim$all_info[[1]]
					)


d <- clone_traj(sim, delta = delta.list[1], max_t = max_t)
for(i.delta in delta.list[2:length(delta.list)]) d <- rbind(d, clone_traj(sim, delta = i.delta, max_t = max_t))
d$resistant <- d$resistant.alate + d$resistant.apterous + d$resistant.parasitized
d$susceptible <- d$susceptible.alate + d$susceptible.apterous + d$susceptible.parasitized

d <- aggregate(cbind(resistant, susceptible, wasps) ~ time + delta, data=d, FUN = sum)
d <- d[d$time > 100,]
d$prop <- d$resistant/(d$susceptible + d$resistant)

# S9 ----
pdf("Nell fig S9.pdf", width = 10, height = 4)
	par(mfrow=c(1,3))

	dd <- d[d$delta == delta.list[1],]
	lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps, colkey = FALSE, theta = 0, phi = 0, xlab = "Resistant", ylab = "Susceptible", zlab = "Wasps")
	#points3D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], dd$wasps[nrow(dd)], colkey = FALSE, add = TRUE)
	for(i.delta in delta.list[2:length(delta.list)]) {
		dd <- d[d$delta == i.delta,]
		lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps, colkey = FALSE, add = TRUE)
		#points3D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], dd$wasps[nrow(dd)], colkey = FALSE, add = TRUE)
	}

	dd <- d[d$delta == delta.list[1],]
	lines2D(dd$resistant, dd$susceptible, colvar = dd$wasps, colkey = FALSE, xlab = "Resistant", ylab = "Susceptible")
	points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
	for(i.delta in delta.list[2:length(delta.list)]) {
		dd <- d[d$delta == i.delta,]
		lines2D(dd$resistant, dd$susceptible, add = TRUE, colvar = dd$wasps, colkey = FALSE)
		points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
	}

	dd <- d[d$delta == delta.list[1],]
	lines2D(dd$time, dd$prop, colvar = dd$wasps, colkey = FALSE, xlab = "Time", ylab = "Proportion resistant")
	for(i.delta in delta.list[2:length(delta.list)]) {
		dd <- d[d$delta == i.delta,]
		lines2D(dd$time, dd$prop, add = TRUE, colvar = dd$wasps, colkey = FALSE)
	}
dev.off()

########################*
# domain of attraction
max_t <- 1e4
tol <- 1e-3

# baseline
disp.list <- c(.01*(1:25))
disp.list <- c(.0001*(1:9), disp.list, .258 + .0001*(1:6), .26, .28)

pick.disp <- 25
delta.list <- exp(.25*(-20:20))
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
						alate_field_disp_p = disp.list[pick.disp],
						extinct_N = 1e-10,
                        max_t = max_t, save_every = 1e0)
#sim <- restart_experiment(sim, new_starts = sim$all_info[[1]], max_t = max_t, alate_field_disp_p = disp.list[pick.disp])

# S8 ----
pdf("Nell fig S8.pdf", width = 5, height = 7)
	clone_plot(sim, delta.list, max_t = 500)
dev.off()

clone_eig(sim)

df <- data.frame(disp = rep(disp.list, each = length(delta.list)), delta = rep(delta.list, times = length(disp.list)), peak.resistant = NA, peak.susceptible = NA, peak.wasps = NA, trough.resistant = NA, trough.susceptible = NA, trough.wasps = NA, prop.delta = NA, converge = NA, peak.prop = NA, trough.prop = NA)
counter <- 0

new.sim <- sim
new.starts <- new.sim$all_info[[1]]
for(count.disp in pick.disp:length(disp.list)){
# for(count.disp in (pick.disp-1):1){
	i.disp <- disp.list[count.disp]
	show(i.disp)

	old.sim <- new.sim
	new.sim <- restart_experiment(old.sim, new_starts = new.starts, alate_field_disp_p = i.disp, max_t = max_t)
	new.starts <- new.sim$all_info[[1]]

	S.resistant <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "resistant"])
	S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "susceptible"])
	test <- (S.resistant/(S.resistant + S.susceptible) > tol)
	if(!test) {
		for(i.delta in rev(delta.list)) {
			test.starts <- old.sim$all_info[[1]]
			pick <- (!is.na(test.starts$line)) & (test.starts$line == "resistant")
			test.starts[pick, "N"] <- i.delta * test.starts[pick, "N"]

			test.sim <- restart_experiment(old.sim, new_starts = test.starts, alate_field_disp_p = i.disp, max_t = max_t)
			test.starts <- test.sim$all_info[[1]]
			S.resistant <- sum(test.starts$N[(!is.na(test.starts$line)) & (test.starts$line == "resistant")])
			S.susceptible <- sum(test.starts$N[(!is.na(test.starts$line)) & (test.starts$line == "susceptible")])

			show(c(i.disp, i.delta, S.resistant/(S.resistant + S.susceptible)))

			new.test <- (S.resistant/(S.resistant + S.susceptible) > tol)
			if(new.test) break
		}
		new.sim <- test.sim
		new.starts <- test.sim$all_info[[1]]
	}
	S.resistant <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "resistant"])
	S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "susceptible"])

	# pdf(paste0("sim.",i.disp,".pdf"), height = 6, width = 4)
	# 	clone_plot(sim = new.sim, delta.list = delta.list, labels = "", max_t = 500)
	# dev.off()

	# aphids
	d <- new.sim$aphids
	d <- d[d$type != "mummy",]
	d$line.type <- paste0(d$line,".", d$type)
	d <- d[,-c(1,5,6)]
	d <- spread(d, "line.type", "N")
	d <- d[d$time >= 9000,]

	d$resistant <- d$resistant.alate + d$resistant.apterous + d$resistant.parasitized
	d$susceptible <- d$susceptible.alate + d$susceptible.apterous + d$susceptible.parasitized

	d <- aggregate(cbind(resistant, susceptible) ~ time, data=d, FUN = sum)
	d$prop <- d$resistant/(d$resistant + d$susceptible)

	df$peak.prop[df$disp == i.disp] <- max(d$prop)
	df$trough.prop[df$disp == i.disp] <- min(d$prop)
	df$peak.resistant[df$disp == i.disp] <- max(d$resistant)
	df$peak.susceptible[df$disp == i.disp] <- max(d$susceptible)
	df$trough.resistant[df$disp == i.disp] <- min(d$resistant)
	df$trough.susceptible[df$disp == i.disp] <- min(d$susceptible)

	# wasps
	d <- new.sim$wasps
	d <- d[d$time >= 9000,]

	d <- aggregate(wasps ~ time, data=d, FUN = sum)
	df$peak.wasps[df$disp == i.disp] <- max(d$wasps)
	df$trough.wasps[df$disp == i.disp] <- min(d$wasps)

	eig <- clone_eig(new.sim)$value
	df$eig[df$disp == i.disp] <- eig

	for(i.delta in delta.list) {
		df$converge[df$disp == i.disp & df$delta == i.delta] <- clone_converge(new.sim, i.delta, max_t = max_t)
		df$prop.delta[df$disp == i.disp & df$delta == i.delta] <- i.delta*sum(S.resistant)/(i.delta*sum(S.resistant) + sum(S.susceptible))
	}
	show(df[df$disp == i.disp,])
}
df.save <- df

df <- rbind(df, df.save)
write.table(df, "df_Fig4ab_26Oct22.csv", sep=",", row.names=FALSE)

########################*
df <- read.csv("df_Fig4ab_26Oct22.csv")


########################*
# crude figure
df$col.converge <- df$converge+1
df$lty.eig <- 2 - (Im(df$eig) == 0)


par(mfrow=c(1,1))
plot(prop.delta ~ disp, data = df, col = col.converge, pch=3)

points(prop.delta ~ disp, data = df[df$delta == 1 & df$lty.eig == 2,], lty = lty.eig)
points(peak.prop ~ disp, data = df[df$delta == 1,], lty = lty.eig)
points(trough.prop ~ disp, data = df[df$delta == 1,], lty = lty.eig)

########################*
# pretty figure
w <- data.frame(disp = unique(df$disp), stable.equil1 = NA, stable.equil2 = NA, unstable.equil = NA, upper = NA)
for(i.disp in rev(w$disp)){
	ww <- df[df$disp == i.disp,]

	w$stable.equil1[w$disp == i.disp] <- ww$peak.prop[1]
	w$stable.equil2[w$disp == i.disp] <- ww$trough.prop[1]

	if(any(ww$converge[!is.na(ww$converge)])){
		www <- ww[ww$converge == TRUE,]
		w$unstable.equil[w$disp == i.disp] <- www$prop.delta[1]
		w$upper[w$disp == i.disp] <- www$prop.delta[nrow(www)]
	}
	w$peak.resistant[w$disp == i.disp] <- ww$peak.resistant[1]
	w$peak.susceptible[w$disp == i.disp] <- ww$peak.susceptible[1]
	w$peak.wasps[w$disp == i.disp] <- ww$peak.wasps[1]
	w$trough.resistant[w$disp == i.disp] <- ww$trough.resistant[1]
	w$trough.susceptible[w$disp == i.disp] <- ww$trough.susceptible[1]
	w$trough.wasps[w$disp == i.disp] <- ww$trough.wasps[1]

}
w <- w[w$disp >= 0.0005,]

#w <- w[-nrow(w),]


w$trough.resistant <- log10(w$trough.resistant + .0001) #/max(w$peak.resistant, na.rm = TRUE)
w$peak.resistant <- log10(w$peak.resistant + .0001) #/max(w$peak.resistant, na.rm = TRUE)
w$trough.susceptible <- log10(w$trough.susceptible + .0001) #/max(w$peak.susceptible, na.rm = TRUE)
w$peak.susceptible <- log10(w$peak.susceptible + .0001) #/max(w$peak.susceptible, na.rm = TRUE)
w$trough.wasps <- log10(w$trough.wasps + .0001) #/max(w$peak.wasps, na.rm = TRUE)
w$peak.wasps <- log10(w$peak.wasps + .0001) #/max(w$peak.wasps, na.rm = TRUE)

w.spline <- smooth.spline(x = w$disp, y = w$unstable.equil, df = 6)
w$unstable.equil.spline <- predict(w.spline, w$disp)$y

w$stable.equil1[w$disp >= 0.2586] <- NA
w$stable.equil2[w$disp >= 0.2586] <- NA
w$peak.resistant[w$disp > 0.2586] <- NA
w$trough.resistant[w$disp > 0.2586] <- NA
w$unstable.equil.spline[w$disp >= 0.2586] <- NA


#
ww <- w[w$disp <= 0.2585,]

ww <- rbind(ww, rep(NA, ncol(ww)), rep(NA, ncol(ww)))
ww$disp[nrow(ww)-1] <- w$disp[w$disp == 0.2585]
ww$upper[nrow(ww)-1] <- 0
ww$disp[nrow(ww)] <- .28
ww$upper[nrow(ww)] <- 0

# 4ab ----
pdf("Fig4ab 26Oct22.pdf", height = 8, width = 6)
	par(mfrow = c(2,1), mai=c(1,1,.1,.1))
	plot(peak.resistant ~ disp, data = w, typ="l", ylab = "Abundance", xlab = "Aphid dispersal", col="#CCCC00", lwd = 2, xlim = c(0,.28), ylim = c(-4,4), yaxt = "n")
	mtext(side = 2, at=2*(-2:2), text=c(expression(10^-4), expression(10^-2), expression(10^0), expression(10^2), expression(10^4)), las=2, adj=1.2)
	lines(trough.resistant ~ disp, data = w, col="#CCCC00", lwd = 2)
	lines(peak.susceptible ~ disp, data = w, col="blue", lwd = 2)
	lines(trough.susceptible ~ disp, data = w, col="blue", lwd = 2)
	lines(peak.wasps ~ disp, data = w, col="red", lwd = 2)
	lines(trough.wasps ~ disp, data = w, col="red", lwd = 2)

	w <- w[w$disp <= 0.2585,]
	plot(stable.equil1 ~ disp, data = w, typ="l", xlim = c(0,.28), ylim = c(0,1), ylab = "Proportion resistant", xlab = "Aphid dispersal")
	#conf_bounds(x = w$disp, y.lower = w$unstable.equil.spline, y.upper = w$upper, col="lightgray")

	conf_bounds(x = ww$disp, y.lower = ww$upper, y.upper = rep(1,nrow(ww)), col="lightgray")
	conf_bounds(x = w$disp, y.upper = w$unstable.equil.spline, y.lower = rep(0,nrow(w)), col="lightgray")

	lines(stable.equil1 ~ disp, data = w, lwd = 2, col="green")
	lines(stable.equil2 ~ disp, data = w, lwd = 2, col="green")
	points(stable.equil1 ~ disp, data = w[w$disp > .04 & w$disp <= .15,], col="green")
	points(stable.equil2 ~ disp, data = w[w$disp > .04 & w$disp <= .15,], col="green")
	lines(unstable.equil.spline ~ disp, data = w, col = "black", lwd = 2)

dev.off()



####################################################################*
# n fields
####################################################################*

# for wasp.disp = 0
sim0 <- sim

delta.list <- exp(.5*(-10:10))

# baseline
wasp.disp.list <- sort(c(.002*(0:25), .001,.003))
pick.wasp.disp <- 10
wasp.disp <- wasp.disp.list[pick.wasp.disp]
# wasp.disp <- 0

disp <- .2
a <- 2.32
#a <- 2.32/4

n_fields <- 28
save_every <- 1
wasp_density_0 <- .1*rep(1,n_fields)
# wasp_density_0[1] <- 0

harvesting.length <- 28
day.interval <- harvesting.length/n_fields
# max_t needs to be a multiple of the harvesting length
max_t <- harvesting.length*500
tol <- 1e-3

extinct_N <- 1e-10

s <- 0.025
min.surv <- s
max.surv <- s
n.events <- round(max_t/day.interval)

perturb <- data.frame(when=rep(day.interval*(1:n.events), each=3), where=rep(1:n_fields, each=3), who=c("resistant","susceptible","mummies"), how=runif(3*n.events*n_fields, min=min.surv, max=max.surv))
# kill all mummies
perturb$how[perturb$who == "mummies"] <- 0

# kill all resistant
# perturb$how[perturb$who == "resistant"] <- 0 * perturb$how[perturb$who == "resistant"]

perturb <- perturb[1:(3*n.events),]

sim <- sim_experiments(clonal_lines = c(line_s, line_r),
						n_fields = n_fields,
						wasp_density_0 = wasp_density_0,
						alate_field_disp_p = disp,
						wasp_disp_m0 = wasp.disp,
						perturb = perturb,
						a = a,
						K = 12500,
						#h = 0,
						#extinct_N = 1e-10,
						extinct_N = extinct_N,
                        max_t = max_t, save_every = save_every)
sim <- rm_tibs(sim)

pdf("n field trajectories.pdf", width = 6, height = 6)
	par(mfrow=c(2,1), mai=c(.4,.4,.6,.1))
	# w <- sim$aphids[sim$aphids$ time %% harvesting.length == 0,]
	# ww <- sim$wasps[sim$wasps$ time %% harvesting.length == 0,]
	w <- sim$aphids[sim$aphids$time %% harvesting.length == 0 & sim$aphids$time > harvesting.length*400,]
	ww <- sim$wasps[sim$wasps$time %% harvesting.length == 0 & sim$wasps$time > harvesting.length*400,]

	w$time0 <- w$time - min(w$time)
	ww$time0 <- ww$time - min(ww$time)

	i.field <- 1
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="", ylim = c(.001,10000), xaxt = "n", lwd = 2)
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)

	i.field <- 2
		par(mai=c(.9,.4,.1,.1))
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="Days", ylim = c(.001,10000), lwd = 2)
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)

dev.off()

########################*
# phase portrait
df <- data.frame(wasp.disp = rep(wasp.disp.list, each = length(delta.list)), delta = rep(delta.list, times = length(wasp.disp.list)), prop.delta = NA, converge = NA, peak.resistant1 = NA, peak.resistant2 = NA, trough.resistant1 = NA, trough.resistant2 = NA, peak.susceptible1 = NA, peak.susceptible2 = NA, trough.susceptible1 = NA, trough.susceptible2 = NA, peak.wasps1 = NA, peak.wasps2 = NA, trough.wasps1 = NA, trough.wasps2 = NA, peak.prop1 = NA, peak.prop2 = NA, trough.prop1 = NA, trough.prop2 = NA)

for(case in 1:2){
	if(case == 1) {
		count.list <- pick.wasp.disp:length(wasp.disp.list)
	}else{
		sim <- sim0
		count.list <- (pick.wasp.disp-1):1
	}
	new.sim <- sim
	new.starts <- new.sim$all_info[[1]]
	for(count.wasp.disp in count.list){
		i.wasp.disp <- wasp.disp.list[count.wasp.disp]
		show(i.wasp.disp)

		old.sim <- new.sim
		new.sim <- restart_experiment(old.sim, new_starts = new.starts, max_t = max_t, perturb = perturb, wasp_disp_m0 = i.wasp.disp)
		new.starts <- new.sim$all_info[[1]]
		#new.sim <- restart_experiment(sim, max_t = max_t, perturb = perturb, wasp_disp_m0 = i.wasp.disp)

		S.resistant <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "resistant"])
		S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) & new.starts$line == "susceptible"])

		# pdf(paste0("sim.multifield.",i.wasp.disp,".pdf"), height = 6, width = 4)
			# clone_plot(sim = new.sim, delta.list = delta.list, labels = "", max_t = harvesting.length * 40, perturb = perturb)
		# dev.off()

		# aphids
		d <- new.sim$aphids
		d <- d[d$type != "mummy",]
		d$line.type <- paste0(d$line,".", d$type)
		d <- d[,-c(1,5,6)]
		d <- spread(d, "line.type", "N")
		d <- d[d$time >= harvesting.length*400,]
		d <- d[d$time %% day.interval == 0,]

		d$resistant <- d$resistant.alate + d$resistant.apterous + d$resistant.parasitized
		d$susceptible <- d$susceptible.alate + d$susceptible.apterous + d$susceptible.parasitized

		d <- aggregate(cbind(resistant, susceptible) ~ time + field, data=d, FUN = sum)
		d$prop <- d$resistant/(d$resistant + d$susceptible)

		# wasps
		dw <- new.sim$wasps
		dw <- dw[dw$time >= harvesting.length*400,]
		dw <- dw[dw$time %% day.interval == 0,]

		dw <- aggregate(wasps ~ time + field, data=dw, FUN = sum)

		d <- merge(d, dw)

		max.matrix <- matrix(NA,n_fields,4)
		min.matrix <- matrix(NA,n_fields,4)
		for(i.field in 1:n_fields){
			max.matrix[i.field,1] <- max(d$resistant[d$field == i.field])
			min.matrix[i.field,1] <- min(d$resistant[d$field == i.field])

			max.matrix[i.field,2] <- max(d$susceptible[d$field == i.field])
			min.matrix[i.field,2] <- min(d$susceptible[d$field == i.field])

			max.matrix[i.field,3] <- max(d$wasps[d$field == i.field])
			min.matrix[i.field,3] <- min(d$wasps[d$field == i.field])

			max.matrix[i.field,4] <- max(d$prop[d$field == i.field])
			min.matrix[i.field,4] <- min(d$prop[d$field == i.field])
		}

		df$peak.resistant1[df$wasp.disp == i.wasp.disp] <- max(max.matrix[,1])
		df$peak.resistant2[df$wasp.disp == i.wasp.disp] <- min(max.matrix[,1])

		df$peak.susceptible1[df$wasp.disp == i.wasp.disp] <- max(max.matrix[,2])
		df$peak.susceptible2[df$wasp.disp == i.wasp.disp] <- min(max.matrix[,2])

		df$peak.wasps1[df$wasp.disp == i.wasp.disp] <- max(max.matrix[,3])
		df$peak.wasps2[df$wasp.disp == i.wasp.disp] <- min(max.matrix[,3])

		df$peak.prop1[df$wasp.disp == i.wasp.disp] <- max(max.matrix[,4])
		df$peak.prop2[df$wasp.disp == i.wasp.disp] <- min(max.matrix[,4])

		df$trough.resistant1[df$wasp.disp == i.wasp.disp] <- max(min.matrix[,1])
		df$trough.resistant2[df$wasp.disp == i.wasp.disp] <- min(min.matrix[,1])

		df$trough.susceptible1[df$wasp.disp == i.wasp.disp] <- max(min.matrix[,2])
		df$trough.susceptible2[df$wasp.disp == i.wasp.disp] <- min(min.matrix[,2])

		df$trough.wasps1[df$wasp.disp == i.wasp.disp] <- max(min.matrix[,3])
		df$trough.wasps2[df$wasp.disp == i.wasp.disp] <- min(min.matrix[,3])

		df$trough.prop1[df$wasp.disp == i.wasp.disp] <- max(min.matrix[,4])
		df$trough.prop2[df$wasp.disp == i.wasp.disp] <- min(min.matrix[,4])

		for(i.delta in delta.list) {
			df$converge[df$wasp.disp == i.wasp.disp & df$delta == i.delta] <- clone_wasp_converge(new.sim, i.delta, max_t = max_t, perturb = perturb, tol = extinct_N, day.interval = day.interval)
			df$prop.delta[df$wasp.disp == i.wasp.disp & df$delta == i.delta] <- i.delta*sum(S.resistant)/(i.delta*sum(S.resistant) + sum(S.susceptible))
		}
		show(df[df$wasp.disp == i.wasp.disp,])

		par(mfrow=c(7,1), mai=c(.3,.5,.1,.1))
		w <- new.sim$aphids[new.sim$aphids$time %% harvesting.length == 0 & new.sim$aphids$time > harvesting.length*400,]
		ww <- new.sim$wasps[new.sim$wasps$time %% harvesting.length == 0 & new.sim$wasps$time > harvesting.length*400,]
		for(i.field in 1:7){
			plot(N ~ time, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="blue", log="y", xlab="", ylim = c(.01,10000))
			lines(N ~ time, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",])
			lines(wasps ~ time, data=ww[ww$field == i.field,], col="red")
		}

	}
}

write.table(df, "df_Fig4cd_26Oct22.csv", sep=",", row.names=FALSE)


########################*
df <- read.csv("df_Fig4cd_26Oct22.csv")


########################*
# pretty figure

# df <- df[df$wasp.disp <= .18,]
# df <- df[df$wasp.disp >=.01,]
w <- data.frame(wasp.disp = unique(df$wasp.disp))
for(i.resist in c(0,1)) for(i.wasp.disp in w$wasp.disp){

	ww <- df[df$wasp.disp == i.wasp.disp,]

	w$peak.resistant1[w$wasp.disp == i.wasp.disp] <- ww$peak.resistant1[1]
	w$peak.susceptible1[w$wasp.disp == i.wasp.disp] <- ww$peak.susceptible1[1]
	w$peak.wasps1[w$wasp.disp == i.wasp.disp] <- ww$peak.wasps1[1]
	w$trough.resistant1[w$wasp.disp == i.wasp.disp] <- ww$trough.resistant1[1]
	w$trough.susceptible1[w$wasp.disp == i.wasp.disp] <- ww$trough.susceptible1[1]
	w$trough.wasps1[w$wasp.disp == i.wasp.disp] <- ww$trough.wasps1[1]

	w$peak.resistant2[w$wasp.disp == i.wasp.disp] <- ww$peak.resistant2[1]
	w$peak.susceptible2[w$wasp.disp == i.wasp.disp] <- ww$peak.susceptible2[1]
	w$peak.wasps2[w$wasp.disp == i.wasp.disp] <- ww$peak.wasps2[1]
	w$trough.resistant2[w$wasp.disp == i.wasp.disp] <- ww$trough.resistant2[1]
	w$trough.susceptible2[w$wasp.disp == i.wasp.disp] <- ww$trough.susceptible2[1]
	w$trough.wasps2[w$wasp.disp == i.wasp.disp] <- ww$trough.wasps2[1]

	w$peak.prop1[w$wasp.disp == i.wasp.disp] <- ww$peak.prop1[1]
	w$peak.prop2[w$wasp.disp == i.wasp.disp] <- ww$peak.prop2[1]
	w$trough.prop1[w$wasp.disp == i.wasp.disp] <- ww$trough.prop1[1]
	w$trough.prop2[w$wasp.disp == i.wasp.disp] <- ww$trough.prop2[1]

	www <- ww[ww$converge == TRUE,]
	if(nrow(www) > 0){
		w$lower[w$wasp.disp == i.wasp.disp] <- www$prop.delta[1]
		w$upper[w$wasp.disp == i.wasp.disp] <- www$prop.delta[nrow(www)]
	}
}
w

for(i in 1:3){
	# w[,i + 1 + 3*(0:3)] <- log10(w[,i + 1 + 3*(0:3)]/max(w[,i + 1 + 3*(0:3)], na.rm = TRUE)+.0001)
	w[,i + 1 + 3*(0:3)] <- log10(w[,i + 1 + 3*(0:3)]+.0001)
}
w[,2 + 3*(0:3)]

w$upper[is.na(w$upper)] <- 0
w$upper[23] <- .85


# w.spline <- smooth.spline(x = w$s, y = w$unstable.equil, df = 6)
# w$unstable.equil.spline <- predict(w.spline, w$s)$y

# 4cd ----
pdf("Fig4cd 26Oct22.pdf", height = 8, width = 6)
	par(mfrow = c(2,1), mai=c(1,1,.1,.1))

	plot(peak.resistant1 ~ wasp.disp, data = w, typ="l", ylab = "Abundance", xlab = "Wasp dispersal", col="#CCCC00", lty = 2-w$resistance.flag, ylim = c(-4,4), yaxt = "n")
	mtext(side = 2, at=2*(-2:2), text=c(expression(10^-4), expression(10^-2), expression(10^0), expression(10^2), expression(10^4)), las=2, adj=1.2)

	col <- t_col("#CCCC00", percent = 70)
	conf_bounds(x = w$wasp.disp, y.lower = w$peak.resistant2, y.upper = w$peak.resistant1, col=col)
	conf_bounds(x = w$wasp.disp, y.lower = w$trough.resistant2, y.upper = w$trough.resistant1, col=col)
	lines(peak.resistant1 ~ wasp.disp, data = w, col="#999900")
	lines(peak.resistant2 ~ wasp.disp, data = w, col="#999900")
	lines(trough.resistant1 ~ wasp.disp, data = w, col="#999900", lty = 2)
	lines(trough.resistant2 ~ wasp.disp, data = w, col="#999900", lty = 2)

	col <- t_col("blue", percent = 70)
	conf_bounds(x = w$wasp.disp, y.lower = w$peak.susceptible2, y.upper = w$peak.susceptible1, col=col)
	conf_bounds(x = w$wasp.disp, y.lower = w$trough.susceptible2, y.upper = w$trough.susceptible1, col=col)
	lines(peak.susceptible1 ~ wasp.disp, data = w, col="blue")
	lines(peak.susceptible2 ~ wasp.disp, data = w, col="blue")
	lines(trough.susceptible1 ~ wasp.disp, data = w, col="blue", lty = 2)
	lines(trough.susceptible2 ~ wasp.disp, data = w, col="blue", lty = 2)

	col <- t_col("red", percent = 70)
	conf_bounds(x = w$wasp.disp, y.lower = w$peak.wasps2, y.upper = w$peak.wasps1, col=col)
	conf_bounds(x = w$wasp.disp, y.lower = w$trough.wasps2, y.upper = w$trough.wasps1, col=col)
	lines(peak.wasps1 ~ wasp.disp, data = w, col="red")
	lines(peak.wasps2 ~ wasp.disp, data = w, col="red")
	lines(trough.wasps1 ~ wasp.disp, data = w, col="red", lty = 2)
	lines(trough.wasps2 ~ wasp.disp, data = w, col="red", lty = 2)

	plot(peak.prop1 ~ wasp.disp, data = w, typ="l", ylim = c(0,1), ylab = "Proportion resistant", xlab = "Wasp dispersal")

	conf_bounds(x = w$wasp.disp[!is.na(w$upper)], y.lower = w$upper[!is.na(w$upper)], y.upper = rep(1,sum(!is.na(w$upper))), col="lightgray")
	conf_bounds(x = w$wasp.disp[!is.na(w$lower)], y.upper = w$lower[!is.na(w$lower)], y.lower = rep(0,sum(!is.na(w$lower))), col="lightgray")
	conf_bounds(x = w$wasp.disp, y.lower = w$peak.prop2, y.upper = w$peak.prop1, col="lightgreen")
	conf_bounds(x = w$wasp.disp, y.lower = w$trough.prop2, y.upper = w$trough.prop1, col="lightgreen")

	lines(peak.prop1 ~ wasp.disp, data = w, col="darkgreen")
	lines(peak.prop2 ~ wasp.disp, data = w, col="darkgreen")
	lines(trough.prop1 ~ wasp.disp, data = w, lty = 2, col="darkgreen")
	lines(trough.prop2 ~ wasp.disp, data = w, lty = 2, col="darkgreen")


dev.off()

####################################################################*
# Figuring out alternative states in n fields
####################################################################*
# line_s <- clonal_line("susceptible",
#                       density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
#                       surv_juv_apterous = "high",
#                       surv_adult_apterous = "high",
#                       repro_apterous = "high")

delta.list <- exp(.5*(-10:10))

# baseline
wasp.disp.list <- sort(c(.002*(0:25), .001,.003))
pick.wasp.disp <- 10
wasp.disp <- wasp.disp.list[pick.wasp.disp]

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

# S10 ----
pdf("Nell fig S10.pdf", width = 7, height = 7)
	for(remove in c(TRUE,FALSE)){

		if(remove) {
			par(mfrow=c(2,1), mai=c(.1,1,.8,.1))
			xaxt <- "n"
		} else {
			par(mai=c(.8,1,.1,.1))
			xaxt <- NULL
		}

		if(remove) density_0 <- cbind(c(0,0,0,0,0), rep(0, 5)) else density_0 <- cbind(c(0,0,0,0,32), rep(0, 5))
		line_r_pert <- clonal_line("resistant",
	                      density_0 = density_0,
	                      resistant = TRUE,
	                      surv_paras = 0.57,
	                      surv_juv_apterous = "low",
	                      surv_adult_apterous = "low",
	                      repro_apterous = "low")

		sim <- sim_experiments(clonal_lines = c(line_s, line_r_pert),
								n_fields = n_fields,
								wasp_density_0 = wasp_density_0,
								alate_field_disp_p = disp,
								wasp_disp_m0 = wasp.disp,
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

		i.field <- 2
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="#CCCC00", log="y", xlab="Days", ylim = c(.00001,10000), xaxt = xaxt, lwd = 2, ylab = "Abundance")
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)
		for(i.field in c(1,3,4)){
			lines(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], col="#CCCC00", lwd = 2, lty=3)
			lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2, lty=3)
			lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2, lty=3)
		}
	}
dev.off()


par(mfrow=c(4,1), mai=c(.1,.1,.1,.1))
for(i.field in 1:4){
	plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="", ylim = c(.00001,10000), xaxt = "n", lwd = 2)
	lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
	lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)
}


pdf("n field trajectories.pdf", width = 6, height = 6)
	par(mfrow=c(2,1), mai=c(.4,.4,.6,.1))
	# w <- sim$aphids[sim$aphids$ time %% harvesting.length == 0,]
	# ww <- sim$wasps[sim$wasps$ time %% harvesting.length == 0,]
	w <- sim$aphids[sim$aphids$time %% harvesting.length == 0 & sim$aphids$time > harvesting.length*400,]
	ww <- sim$wasps[sim$wasps$time %% harvesting.length == 0 & sim$wasps$time > harvesting.length*400,]

	w$time0 <- w$time - min(w$time)
	ww$time0 <- ww$time - min(ww$time)

	i.field <- 1
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="", ylim = c(.001,10000), xaxt = "n", lwd = 2)
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)

	i.field <- 2
		par(mai=c(.9,.4,.1,.1))
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="Days", ylim = c(.001,10000), lwd = 2)
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)

dev.off()


####################################################################*
# Figuring out why high dispersal causes resistant clone to go extinct
####################################################################*
# line_s <- clonal_line("susceptible",
#                       density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
#                       surv_juv_apterous = "high",
#                       surv_adult_apterous = "high",
#                       repro_apterous = "high")
#
# line_r <- clonal_line("resistant",
#                       density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
#                       resistant = TRUE,
#                       surv_paras = 0.57,
#                       surv_juv_apterous = "low",
#                       surv_adult_apterous = "low",
#                       repro_apterous = "low")

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

wasp.disp.list[pick.wasp.disp.list]

pick.wasp.disp.list

# S11 ----

pdf("Nell fig S11.pdf", width = 7, height = 7)
	for(pick.wasp.disp in pick.wasp.disp.list){

		wasp.disp <- wasp.disp.list[pick.wasp.disp]
		if(pick.wasp.disp == pick.wasp.disp.list[1]) {
			par(mfrow=c(2,1), mai=c(.1,1,.8,.1))
			xaxt <- "n"
		} else {
			par(mai=c(.8,1,.1,.1))
			xaxt <- NULL
		}

		sim <- sim_experiments(clonal_lines = c(line_s, line_r),
								n_fields = n_fields,
								wasp_density_0 = wasp_density_0,
								alate_field_disp_p = disp,
								wasp_disp_m0 = wasp.disp,
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

		i.field <- 2
		plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="#CCCC00", log="y", xlab="Days", ylim = c(.00001,10000), xaxt = xaxt, lwd = 2, ylab = "Abundance")
		lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
		lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)
		for(i.field in c(1,3,4)){
			lines(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], col="#CCCC00", lwd = 2, lty=3)
			lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2, lty=3)
			lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2, lty=3)
		}
	}
dev.off()


par(mfrow=c(4,1), mai=c(.1,.1,.1,.1))
for(i.field in 1:4){
	plot(N ~ time0, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col="green", log="y", xlab="", ylim = c(.00001,10000), xaxt = "n", lwd = 2)
	lines(N ~ time0, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",], col="blue", lwd = 2)
	lines(wasps ~ time0, data=ww[ww$field == i.field,], col="red", lwd = 2)
}



