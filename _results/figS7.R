
library(gameofclones)
library(grid)
library(here)


# File name for figure:
file_out <- here("_results/_plots/trajectories-experiment.pdf")



rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
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
clone_plot <- function(sim, delta.list, labels = "", max_t = 1e3,
                       max_t.converge = 1e4, perturb = NULL){

    d <- sim$aphids
    d$line.type <- paste0(d$line,".", d$type)
    d$line.type[d$line.type == "NA.mummy"] <- "mummy"
    d <- d[,-c(1,5,6)]
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
    grid.text("A", x = unit(0, "npc"), y = unit(1, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 18, fontface = "bold"))
    for (i.delta in delta.list) {
        lines(prop ~ time, data=d[d$delta == i.delta & d$field == 1,],
              col=col, lty=2-as.numeric(converge))
    }
    text(x=max_t/2, y=.9, labels=labels, cex=2)
    lines(prop ~ time, data=d[d$delta == 1 & d$field == 1,])

    plot(prop ~ time, data=d[d$delta == 1 & d$field == 2,], typ="l",
         ylim=c(0,1), ylab = "Resistant (field 2)", xlab = "Time")
    grid.text("B", x = unit(0, "npc"), y = unit(0.66, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 18, fontface = "bold"))
    for (i.delta in delta.list) {
        lines(prop ~ time, data=d[d$delta == i.delta & d$field == 2,],
              col=col, lty=2-as.numeric(converge))
    }
    lines(prop ~ time, data=d[d$delta == 1 & d$field == 2,])

    plot(wasps ~ time, data=d[d$delta == 1 & d$field == 1,], typ="l",
         ylim = c(0, max(d$wasps)), ylab = "Wasps (field 1)", xlab = "Time")
    grid.text("C", x = unit(0, "npc"), y = unit(0.33, "npc"),
              just = c("left", "top"),
              gp = gpar(fontsize = 18, fontface = "bold"))
    for(i.delta in delta.list) {
        lines(wasps ~ time, data=d[d$delta == i.delta & d$field == 1,],
              col=col, lty=2-as.numeric(converge))
    }
    lines(wasps ~ time, data=d[d$delta == 1 & d$field == 1,])
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



# figure itself ----
cairo_pdf(filename = file_out, width = 5, height = 7)
clone_plot(sim, delta.list, max_t = 500)
dev.off()
