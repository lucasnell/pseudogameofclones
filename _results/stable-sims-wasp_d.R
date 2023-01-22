library(gameofclones)
library(scales)
library(viridisLite)
library(here)


# colors for resistant, susceptible, and parasitoid wasps, respectively
col_pal <- list(r = viridis(100)[50],
                s = viridis(100)[95],
                w = viridis(100)[1])

# Directory where plots produced here will be added:
plot_dir_out <- here("_results/_plots/stable-sims")
if (!dir.exists(plot_dir_out)) dir.create(plot_dir_out, recursive = TRUE)
# Names of files produced here:
plots_out <- list(N = paste0(plot_dir_out, "/stable-sims-wasp_d-abundance.pdf"),
                  P = paste0(plot_dir_out, "/stable-sims-wasp_d-resistance.pdf"))
# Name of temporary results file produced here:
tmp_results <- here("_results/_data/stable-sims-wasp_d.csv")



rm_tibs <- function(.sims) {
    for (n in c("aphids", "wasps")) {
        .sims[[n]] <- as.data.frame(.sims[[n]])
    }
    return(.sims)
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

disp.list <- c(.01*(1:25))
disp.list <- c(.0001*(1:9), disp.list, .258 + .0001*(1:6), .26, .28)

pick.disp <- 25
delta.list <- exp(.25*(-20:20))
sim0 <- sim_experiments(clonal_lines = c(line_s, line_r),
                       alate_field_disp_p = disp.list[pick.disp],
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e0)
sim0 <- rm_tibs(sim0)
#

delta.list <- exp(.5*(-10:10))

# baseline
wasp.disp.list <- sort(c(.002*(0:25), .001,.003))
pick.wasp.disp <- 10
wasp.disp <- wasp.disp.list[pick.wasp.disp]

disp <- .2
a <- 2.32

n_fields <- 28
save_every <- 1
wasp_density_0 <- .1*rep(1,n_fields)

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
sim0 <- sim

verbose <- FALSE

# Re-run simulations?
rerun_sims <- !file.exists(tmp_results)



if (rerun_sims) {
    # Takes ~18 min
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
            if (verbose) show(i.wasp.disp)

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
            if (verbose) {
                show(df[df$wasp.disp == i.wasp.disp,])

                par(mfrow=c(7,1), mai=c(.3,.5,.1,.1))
                w <- new.sim$aphids[new.sim$aphids$time %% harvesting.length == 0 & new.sim$aphids$time > harvesting.length*400,]
                ww <- new.sim$wasps[new.sim$wasps$time %% harvesting.length == 0 & new.sim$wasps$time > harvesting.length*400,]
                for(i.field in 1:7){
                    plot(N ~ time, data=w[w$field == i.field & w$line == "resistant" & w$type == "apterous",], typ="l", col=col_pal$s, log="y", xlab="", ylim = c(.01,10000))
                    lines(N ~ time, data=w[w$field == i.field & w$line == "susceptible" & w$type == "apterous",])
                    lines(wasps ~ time, data=ww[ww$field == i.field,], col=col_pal$w)
                }

            }


        }
    }

    write.csv(df, tmp_results, row.names = FALSE)

} else {

    ########################*
    df <- read.csv(tmp_results)

}




# ============================================================================*
# pre-processing for figures ----
# ============================================================================*


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

for(i in 1:3){
    w[,i + 1 + 3*(0:3)] <- log10(w[,i + 1 + 3*(0:3)]+.0001)
}

w$upper[is.na(w$upper)] <- 0
w$upper[23] <- .85





# ============================================================================*
# figures themselves ----
# ============================================================================*


# For equilibrium abundances ~ wasp dispersal

cairo_pdf(filename = plots_out$N, height = 4, width = 6)
{
    par(mai=c(0.9, 0.9, 0.1, 0.1))

    plot(peak.resistant1 ~ wasp.disp, data = w, typ="l", ylab = "Abundance",
         xlab = "Wasp dispersal", col=col_pal$r, lty = 2,
         ylim = c(-4,4), yaxt = "n")
    aty <- axTicks(2)
    y_labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))))
    axis(2,at=aty,labels=y_labels, las=1)

    col <- alpha(col_pal$r, 0.30)
    conf_bounds(x = w$wasp.disp, y.lower = w$peak.resistant2,
                y.upper = w$peak.resistant1, col=col)
    conf_bounds(x = w$wasp.disp, y.lower = w$trough.resistant2,
                y.upper = w$trough.resistant1, col=col)
    lines(peak.resistant1 ~ wasp.disp, data = w, col=col_pal$r, lwd = 2)
    lines(peak.resistant2 ~ wasp.disp, data = w, col=col_pal$r, lwd = 2)
    lines(trough.resistant1 ~ wasp.disp, data = w, col=col_pal$r, lty = 2, lwd = 2)
    lines(trough.resistant2 ~ wasp.disp, data = w, col=col_pal$r, lty = 2, lwd = 2)

    col <- alpha(col_pal$s, 0.30)
    conf_bounds(x = w$wasp.disp, y.lower = w$peak.susceptible2,
                y.upper = w$peak.susceptible1, col=col)
    conf_bounds(x = w$wasp.disp, y.lower = w$trough.susceptible2,
                y.upper = w$trough.susceptible1, col=col)
    lines(peak.susceptible1 ~ wasp.disp, data = w, col=col_pal$s, lwd = 2)
    lines(peak.susceptible2 ~ wasp.disp, data = w, col=col_pal$s, lwd = 2)
    lines(trough.susceptible1 ~ wasp.disp, data = w, col=col_pal$s,
          lty = 2, lwd = 2)
    lines(trough.susceptible2 ~ wasp.disp, data = w, col=col_pal$s,
          lty = 2, lwd = 2)

    col <- alpha(col_pal$w, 0.30)
    conf_bounds(x = w$wasp.disp, y.lower = w$peak.wasps2,
                y.upper = w$peak.wasps1, col=col)
    conf_bounds(x = w$wasp.disp, y.lower = w$trough.wasps2,
                y.upper = w$trough.wasps1, col=col)
    lines(peak.wasps1 ~ wasp.disp, data = w, col=col_pal$w, lwd = 2)
    lines(peak.wasps2 ~ wasp.disp, data = w, col=col_pal$w, lwd = 2)
    lines(trough.wasps1 ~ wasp.disp, data = w, col=col_pal$w, lty = 2, lwd = 2)
    lines(trough.wasps2 ~ wasp.disp, data = w, col=col_pal$w, lty = 2, lwd = 2)

    text(x = max(w$wasp.disp), y = 0.9, labels = "parasitoid", col=col_pal$w,
         adj = c(1, 1), font = 2)
    text(x = 0.039, y = 2.5, labels = "resistant", col=col_pal$r,
         adj = c(0, 0.5), font = 2)
    text(x = max(w$wasp.disp), y = 3.9, labels = "susceptible", col=col_pal$s,
         adj = c(1, 1), font = 2)
}
dev.off()







# For equilibrium proportion resistance ~ wasp dispersal

cairo_pdf(filename = plots_out$P, height = 3.6, width = 5.5)
{
    par(mai=c(0.1, 0.5, 0.5, 0.1))

    plot(peak.prop1 ~ wasp.disp, data = w, typ="l", ylim = c(0,1), xaxt = "n",
         ylab = "", xlab = "")
         # ylab = "Proportion resistant", xlab = "Wasp dispersal")
    axis(3)

    conf_bounds(x = w$wasp.disp[!is.na(w$upper)],
                y.lower = w$upper[!is.na(w$upper)],
                y.upper = rep(1,sum(!is.na(w$upper))), col="lightgray")
    conf_bounds(x = w$wasp.disp[!is.na(w$lower)],
                y.upper = w$lower[!is.na(w$lower)],
                y.lower = rep(0,sum(!is.na(w$lower))), col="lightgray")
    conf_bounds(x = w$wasp.disp, y.lower = w$peak.prop2, y.upper = w$peak.prop1,
                col=alpha("dodgerblue", 0.30))
    conf_bounds(x = w$wasp.disp, y.lower = w$trough.prop2,
                y.upper = w$trough.prop1, col=alpha("dodgerblue", 0.30))

    lines(peak.prop1 ~ wasp.disp, data = w, col="dodgerblue3", lwd=2)
    lines(peak.prop2 ~ wasp.disp, data = w, col="dodgerblue3", lwd=2)
    lines(trough.prop1 ~ wasp.disp, data = w, lty = 2, col="dodgerblue3", lwd=2)
    lines(trough.prop2 ~ wasp.disp, data = w, lty = 2, col="dodgerblue3", lwd=2)
}
dev.off()
