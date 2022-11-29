library(gameofclones)

#'
#' WARNING: Some of the code below contains tibbles.
#' To remove these abominations, you can use this function on the output
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







#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#'
#' BELOW IS THE CODE TONY SENT:
#'
#' ------------------------------------------------------------------------
#' ------------------------------------------------------------------------
#'
#'
#'
#' I’m getting close to having figures for your ms.
#' I want to do a version with multiple fields and wasp dispersal, so I’m
#' using the perturb option in sim_experiments().
#' The perturbations don’t seem to carry into restart_experiments(),
#' unless I’m doing something wrong. Here is example code
#'



disp <- 0
wasp_disp <- 0
#a <- 2.32
a <- .5

n_fields <- 7
wasp_density_0 <- rep(.3, times=n_fields)

harvesting.length <- 28
day.interval <- harvesting.length/n_fields
# max_t needs to be a multiple of the harvesting length
max_t <- harvesting.length*30

min.surv <- .05
max.surv <- .05
n.events <- round(max_t/day.interval)

perturb <- data.frame(when=rep(day.interval*(1:n.events), each=3), where=rep(1:n_fields, each=3), who=c("resistant","susceptible","mummies"), how=runif(3*n.events*n_fields, min=min.surv, max=max.surv))
perturb$how[perturb$who == "mummies"] <- 0

perturb <- perturb[1:(3*n.events),]

sim <- sim_experiments(clonal_lines = c(line_s, line_r),
                       n_fields = n_fields,
                       wasp_density_0 = wasp_density_0,
                       alate_field_disp_p = disp,
                       wasp_disp_p = wasp_disp,
                       perturb = perturb,
                       a = a,
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e0)
sim <- rm_tibs(sim)


new.sim <- restart_experiment(sim, max_t = max_t, perturb = perturb)
new.new.sim <- restart_experiment(new.sim, max_t = max_t, perturb = perturb)


par(mfrow=c(n_fields,1), mai=c(.3,.5,.1,.1))
for(i.field in 1:n_fields){
    plot(N ~ time, data=sim$aphids[sim$aphids$field == i.field & sim$aphids$line == "resistant" & sim$aphids$type == "apterous",], typ="l", col="blue", log="y", xlab="")
    lines(N ~ time, data=sim$aphids[sim$aphids$field == i.field & sim$aphids$line == "susceptible" & sim$aphids$type == "apterous",])
    lines(wasps ~ time, data=sim$wasps[sim$wasps$field == i.field,], col="red")
}

######################
# try restart
new.sim <- restart_experiment(sim, new_starts = sim$all_info[[1]],
                              max_t = max_t, perturb = perturb)

par(mfrow=c(n_fields,1))
for(i.field in 1:n_fields){
    plot(N ~ time, data=new.sim$aphids[new.sim$aphids$field == i.field & new.sim$aphids$line == "resistant" & new.sim$aphids$type == "apterous",], typ="l", col="blue", log="y")
    lines(N ~ time, data=new.sim$aphids[new.sim$aphids$field == i.field & new.sim$aphids$line == "susceptible" & new.sim$aphids$type == "apterous",])
    lines(wasps ~ time, data=new.sim$wasps[new.sim$wasps$field == i.field,], col="red")
}

