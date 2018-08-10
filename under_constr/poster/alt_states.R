
source("under_constr/poster/_preamble.R")



k1 = k2 = 100

a12 = a21 = 2
# Slopes for isocline lines:
m1 = -1 / a12
m2 = -a21

# Non-stable equilibrium coordinates
ns_equil <- c(x = {k1 - a12*k2} / {a12*(m2-m1)},
              y = m2 * ({k1 - a12*k2} / {a12*(m2-m1)}) + k2)
# Stable equilibrium coordinates
equil <- list(sp1 = c(x = k1, y = 0), sp2 = c(x = 0, y = k2))

# ----------
# Creating data frame for red and blue lines (the "zero isoclines")
# ----------
line_df <- data_frame(x = seq(0, 100, length.out = 11)) %>%
    mutate(sp1 = m1 * x + (k1 / a12), sp2 = m2 * x + k2) %>%
    gather(spp, y, sp1:sp2) %>%
    mutate(spp = factor(spp)) %>%
    filter(y >= 0)

# ----------
# Simulating starting points and where they'd go:
# ----------

sim_lv <- function(N01, N02, nt,
                   r = cbind(0.1, 0.1),
                   k = cbind(k1, k2),
                   a = cbind(a12, a21)) {

    N <- matrix(0, nrow = nt + 1, ncol = 2)
    N[1,] <- c(N01, N02)
    for (t in 1:nt) {
        N[t+1,] <- N[t,] + r * N[t,] * (k - {N[t,] + a * rev(N[t,])}) / k
    }

    N <- N %>%
        as_data_frame() %>%
        set_names(c("x", "y")) %>%
        identity()

    return(N)
}

# Creating geoms to add to plot

geom_lotka <- function(N01, N02, nt, arrow_n, ...) {
    lv <- sim_lv(N01 = N01, N02 = N02, nt = nt, ...)
    # geom_out <- lapply(arrow_n,
    #                    function(n_) {
    #                        geom_path(data = eval(lv[1:n_,]),
    #                                  size = 1, linejoin = "mitre",
    #                                  arrow = grid::arrow(type = "closed", angle = 20))
    #                    })
    geom_out <- list()
    geom_out <- c(geom_out,
                  list(geom_path(data = eval(lv), size = 1, linetype = 2),
                       geom_point(data = data_frame(x = N01, y = N02), size = 4,
                                  color = "black", shape = 4)))
    return(geom_out)
}


# ----------
# Creating ggplot object of alternative states
# ----------
as_plot <- ggplot(line_df, aes(x, y)) +
    # Red and blue lines
    geom_line(aes(color = spp), size = 2) +
    # 1:1 line:
    geom_path(data = data_frame(x = c(0, k2 + 5), y = c(0, k2 + 5)),
              size = 1, color = "gray80", linetype = 3) +
    # Hacky axes
    geom_path(data = data_frame(x = c(0, k1 + 10), y = c(0,0)), size = 1) +
    geom_path(data = data_frame(x = c(0,0), y = c(0, k2 + 10)), size = 1) +
    # Trajectories
    geom_lotka(90, 80, 200, c(5, 70)) +
    geom_lotka(90, 45, 200, 10) +
    # geom_lotka(90, 30, 200, 10) +
    geom_lotka(80, 90, 200, c(5, 70)) +
    geom_lotka(45, 90, 100, 5) +
    # geom_lotka(30, 90, 200, 10) +
    geom_lotka(10, 5,  100, 20) +
    geom_lotka(5,  10, 100, 20) +
    # Equilibrium points, first unstable, then stable
    geom_point(data = data_frame(x = ns_equil['x'], y = ns_equil['y']),
               shape = 21, color= 'black', fill = 'white', size = 8) +
    geom_point(data = data_frame(x = c(k1, 0), y = c(0, k2),
                                 species = factor(c('sp1', 'sp2'))),
               aes(color = species), size = 8, shape = 16) +
    # Aesthetics
    scale_color_manual(values = with(palette, c(default_primary, accent)),
                       guide = FALSE)


# Now simplifying it to have no axes:
as_plot <- as_plot +
    theme(axis.line.x = element_blank(),
          axis.line.y.left = element_blank())

as_plot

ggsave(filename = "figs/alt_states.pdf", plot = as_plot,
       width = 16, height = 16, units = "cm", bg = "white", useDingbats = FALSE)


n_ <- 150
as_traj_plot <- sim_lv(90, 45, n_) %>%
    bind_rows() %>%
    mutate(t = 0:n_) %>%
    gather("spp", "N", x:y, factor_key = TRUE) %>%
    ggplot(aes(t, N, color = spp)) +
    geom_line(size = 2) +
    scale_y_continuous(limits = c(0, k1 + 30)) +
    scale_color_manual(values = with(palette, c(default_primary, accent)),
                       guide = FALSE)


ggsave(filename = "figs/alt_states_traj.pdf", plot = as_traj_plot,
       width = 8, height = 8, units = "cm", bg = "white", useDingbats = FALSE)


