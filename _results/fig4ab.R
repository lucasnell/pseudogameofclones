
library(numDeriv)
library(gameofclones)

rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
}
comm <- function(x, sim, new_starts = NULL, max_t = 1){
    if(is.null(new_starts)) new_starts <- sim$all_info[[1]]
    new_starts$N <- x
    restart.sim <- restart_experiment(sim, new_starts = new_starts,
                                      max_t = max_t)$all_info[[1]]
    return(restart.sim$N)
}
clone_eig <- function(sim){
    x <- sim$all_info[[1]]$N
    C <- jacobian(comm, x, sim=sim)
    eig <- eigen(C)
    value <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
    magnitude <- Re(abs(eig$value[abs(eig$values) == max(abs(eig$values))][1]))
    return(list(value = value, magnitude = magnitude))
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
conf_bounds <- function(x, y.lower, y.upper, col="lightgray"){
    polygon(c(x, rev(x)), c(y.upper, rev(y.lower)), col=col, border=NA)
}

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

max_t <- 1e4
tol <- 1e-3

# baseline
disp.list <- c(.01*(1:25))
disp.list <- c(.0001*(1:9), disp.list, .258 + .0001*(1:6), .26, .28)

pick.disp <- 15
delta.list <- exp(.25*(-20:20))
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
                       alate_field_disp_p = disp.list[pick.disp],
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e0)
sim0 <- sim

# clone_eig(sim)

# between here and `shared end` is from figS8 ----

# Re-run simulations?
rerun_sims <- !file.exists("_results/fig4ab.csv")



# NOTE: The starting conditions for each new disp value are based on the previous value of disp. This is necessary, because the stationary point is only locally stable. The most stable stationary point is near the middle of values of disp.list, so the 15th element was selected to start, and the range of disp.list is completed going in both directions from 15.

if (rerun_sims) {

    df <- data.frame(disp = rep(disp.list, each = length(delta.list)),
                     delta = rep(delta.list, times = length(disp.list)),
                     peak.resistant = NA, peak.susceptible = NA, peak.wasps = NA,
                     trough.resistant = NA, trough.susceptible = NA,
                     trough.wasps = NA, prop.delta = NA, converge = NA,
                     peak.prop = NA, trough.prop = NA)
    counter <- 0
    verbose <- FALSE

    for(case in 1:2){
        if(case == 1) {
            count.list <- pick.disp:length(disp.list)
        }else{
            sim <- sim0
            count.list <- (pick.disp-1):1
        }
        new.sim <- sim
        new.starts <- new.sim$all_info[[1]]
        for(count.disp in count.list){

            i.disp <- disp.list[count.disp]
            if (verbose) show(i.disp)

            old.sim <- new.sim
            new.sim <- restart_experiment(old.sim, new_starts = new.starts,
                                          alate_field_disp_p = i.disp, max_t = max_t)
            new.starts <- new.sim$all_info[[1]]

            S.resistant <- sum(new.starts$N[!is.na(new.starts$line) &
                                                new.starts$line == "resistant"])
            S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) &
                                                  new.starts$line == "susceptible"])
            test <- (S.resistant/(S.resistant + S.susceptible) > tol)
            if (!test) {
                for(i.delta in rev(delta.list)) {
                    test.starts <- old.sim$all_info[[1]]
                    pick <- (!is.na(test.starts$line)) &
                        (test.starts$line == "resistant")
                    test.starts[pick, "N"] <- i.delta * test.starts[pick, "N"]

                    test.sim <- restart_experiment(old.sim, new_starts = test.starts,
                                                   alate_field_disp_p = i.disp,
                                                   max_t = max_t)
                    test.starts <- test.sim$all_info[[1]]
                    S.resistant <- sum(test.starts$N[(!is.na(test.starts$line)) &
                                                         (test.starts$line == "resistant")])
                    S.susceptible <- sum(test.starts$N[(!is.na(test.starts$line)) &
                                                           (test.starts$line == "susceptible")])

                    if (verbose) show(c(i.disp, i.delta, S.resistant/
                                            (S.resistant + S.susceptible)))

                    new.test <- (S.resistant/(S.resistant + S.susceptible) > tol)
                    if(new.test) break
                }
                new.sim <- test.sim
                new.starts <- test.sim$all_info[[1]]
            }
            S.resistant <- sum(new.starts$N[!is.na(new.starts$line) &
                                                new.starts$line == "resistant"])
            S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) &
                                                  new.starts$line == "susceptible"])

            # aphids
            d <- new.sim$aphids
            d <- d[d$type != "mummy",]
            d$line.type <- paste0(d$line,".", d$type)
            d <- d[,-c(1,5,6)]
            d <- spread(d, "line.type", "N")
            d <- d[d$time >= 9000,]

            d$resistant <- d$resistant.alate + d$resistant.apterous +
                d$resistant.parasitized
            d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
                d$susceptible.parasitized

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
                df$converge[df$disp == i.disp & df$delta == i.delta] <-
                    clone_converge(new.sim, i.delta, max_t = max_t)
                df$prop.delta[df$disp == i.disp & df$delta == i.delta] <-
                    i.delta * sum(S.resistant) /
                    (i.delta * sum(S.resistant) + sum(S.susceptible))
            }
            if (verbose) show(df[df$disp == i.disp,])
        }
    }
    write.csv(df, "_results/Fig4ab.csv", row.names = FALSE)

} else {

    df <- read.csv("_results/Fig4ab.csv")

}




########################*
# pretty figure
########################*

# pre-processing
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

w$trough.resistant <- log10(w$trough.resistant + .0001)
w$peak.resistant <- log10(w$peak.resistant + .0001)
w$trough.susceptible <- log10(w$trough.susceptible + .0001)
w$peak.susceptible <- log10(w$peak.susceptible + .0001)
w$trough.wasps <- log10(w$trough.wasps + .0001)
w$peak.wasps <- log10(w$peak.wasps + .0001)

ww <- w[!is.na(w$unstable.equil),]

ww.spline <- smooth.spline(x = ww$disp, y = ww$unstable.equil, df = 6)
ww$unstable.equil.spline <- predict(ww.spline, ww$disp)$y

www <- w[w$disp <= 0.2585,]

# This prevents line from continuing to the end for resistant aphids:
w$peak.resistant[tail(which(w$peak.resistant == log10(.0001)), -1)] <- NA
w$trough.resistant[tail(which(w$trough.resistant == log10(.0001)), -1)] <- NA




# ============================================================================*
# figure itself ----
# ============================================================================*

cairo_pdf(filename = "_results/plots/fig_4ab.pdf", height = 8, width = 6)

par(mfrow = c(2,1), mai=c(1,1,.1,.1))
plot(peak.resistant ~ disp, data = w, typ="l", ylab = "Abundance", xlab = "Aphid dispersal", col="#CCCC00", lwd = 2, xlim = c(0,.28), ylim = c(-4,4), yaxt = "n")
mtext(side = 2, at=2*(-2:2), text=c(expression(10^-4), expression(10^-2), expression(10^0), expression(10^2), expression(10^4)), las=2, adj=1.2)
lines(trough.resistant ~ disp, data = w, col="#CCCC00", lwd = 2)
lines(peak.susceptible ~ disp, data = w, col="blue", lwd = 2)
lines(trough.susceptible ~ disp, data = w, col="blue", lwd = 2)
lines(peak.wasps ~ disp, data = w, col="red", lwd = 2)
lines(trough.wasps ~ disp, data = w, col="red", lwd = 2)

plot(stable.equil1 ~ disp, data = www, typ="l", xlim = c(0,.28), ylim = c(0,1), ylab = "Proportion resistant", xlab = "Aphid dispersal")

conf_bounds(x = ww$disp, y.lower = ww$upper, y.upper = rep(1,nrow(ww)), col="lightgray")
conf_bounds(x = ww$disp, y.upper = ww$unstable.equil.spline, y.lower = rep(0,nrow(ww)), col="lightgray")
conf_bounds(x = c(0.2584, .28), y.lower = c(0,0), y.upper = c(1,1), col="lightgray")

lines(stable.equil1 ~ disp, data = www, lwd = 2, col="green")
lines(stable.equil2 ~ disp, data = www, lwd = 2, col="green")
points(stable.equil1 ~ disp, data = www[www$disp > .04 & www$disp <= .15,], col="green")
points(stable.equil2 ~ disp, data = www[www$disp > .04 & www$disp <= .15,], col="green")
lines(unstable.equil.spline ~ disp, data = ww, col = "black", lwd = 2)

dev.off()


