
library(numDeriv)
library(lattice)
library(plot3D)
library(gameofclones)

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
rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
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
for(i.delta in delta.list[2:length(delta.list)]) {
    d <- rbind(d, clone_traj(sim, delta = i.delta, max_t = max_t))
}
d$resistant <- d$resistant.alate + d$resistant.apterous +
    d$resistant.parasitized
d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
    d$susceptible.parasitized

d <- aggregate(cbind(resistant, susceptible, wasps) ~ time + delta, data=d,
               FUN = sum)
d <- d[d$time > 100,]
d$prop <- d$resistant/(d$susceptible + d$resistant)


cairo_pdf(filename = "_results/plots/fig_S9.pdf", width = 10, height = 4)
par(mfrow=c(1,3))

dd <- d[d$delta == delta.list[1],]
lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps,
        colkey = FALSE, theta = 0, phi = 0, xlab = "Resistant",
        ylab = "Susceptible", zlab = "Wasps")
for(i.delta in delta.list[2:length(delta.list)]) {
    dd <- d[d$delta == i.delta,]
    lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps,
            colkey = FALSE, add = TRUE)
}

dd <- d[d$delta == delta.list[1],]
lines2D(dd$resistant, dd$susceptible, colvar = dd$wasps, colkey = FALSE,
        xlab = "Resistant", ylab = "Susceptible")
points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
for(i.delta in delta.list[2:length(delta.list)]) {
    dd <- d[d$delta == i.delta,]
    lines2D(dd$resistant, dd$susceptible, add = TRUE, colvar = dd$wasps,
            colkey = FALSE)
    points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
}

dd <- d[d$delta == delta.list[1],]
lines2D(dd$time, dd$prop, colvar = dd$wasps, colkey = FALSE, xlab = "Time",
        ylab = "Proportion resistant")
for(i.delta in delta.list[2:length(delta.list)]) {
    dd <- d[d$delta == i.delta,]
    lines2D(dd$time, dd$prop, add = TRUE, colvar = dd$wasps, colkey = FALSE)
}
dev.off()
