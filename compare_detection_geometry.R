library(tidyverse)
library(cmdstanr)
source("explore_detection_curvature.R", local = TRUE)

set.seed(20260410)

fit <- as_cmdstan_fit(list.files("tst_folder_sig", full.names = TRUE))
draw_dat <- rebuild_utah_data()

lp_draws <- fit$draws("lp__", format = "draws_matrix")[, 1]
draw_id <- which.max(lp_draws)

ei_draw_tbl <- extract_draw_ei(
  fit = fit,
  draw_id = draw_id,
  TT = draw_dat$TT,
  N_C = draw_dat$N_C
) %>%
  mutate(n_incident_draw = ei_prop * draw_dat$pop_size[county_index])

valid_obs_base <- draw_dat$obs_tbl %>%
  left_join(
    tibble(
      county_index = seq_len(draw_dat$N_C),
      first = draw_dat$first,
      last = draw_dat$last
    ),
    by = "county_index"
  ) %>%
  filter(time >= first, time <= last) %>%
  left_join(ei_draw_tbl, by = c("time", "county_index")) %>%
  mutate(n_incident_draw = pmax(n_incident_draw, 1e-8))

simulate_obs <- function(n_incident, p_true) {
  rbinom(length(n_incident), size = pmax(round(n_incident), 0), prob = p_true)
}

make_surface <- function(ii_obs, p_true, label) {
  alpha_grid <- expand_grid(
    p = seq(0.05, 0.80, by = 0.005),
    alpha = seq(0.45, 1.80, by = 0.01)
  ) %>%
    mutate(
      log_post = purrr::pmap_dbl(
        list(p, alpha),
        function(p, alpha) {
          loglik_obs(
            p = p,
            ii_obs = ii_obs,
            n_incident = alpha * valid_obs_base$n_incident_draw
          )
        }
      ),
      rel_log_post = log_post - max(log_post),
      scenario = label,
      p_true = p_true
    )

  ridge_tbl <- alpha_grid %>%
    group_by(p) %>%
    slice_max(order_by = log_post, n = 1, with_ties = FALSE) %>%
    ungroup()

  alpha_slice <- alpha_grid %>%
    filter(abs(p - p_true) == min(abs(p - p_true))) %>%
    mutate(
      alpha_centered = alpha - 1,
      scenario = label
    )

  width_2ll <- approx(
    x = alpha_slice$rel_log_post,
    y = alpha_slice$alpha,
    xout = -2,
    rule = 2
  )$y

  curvature_alpha <- with(alpha_slice, {
    idx <- which.min(abs(alpha - 1))
    if (idx %in% c(1, nrow(alpha_slice))) {
      NA_real_
    } else {
      h <- alpha[idx + 1] - alpha[idx]
      -(rel_log_post[idx + 1] - 2 * rel_log_post[idx] + rel_log_post[idx - 1]) / (h^2)
    }
  })

  list(
    surface = alpha_grid,
    ridge = ridge_tbl,
    alpha_slice = alpha_slice,
    summary = tibble(
      scenario = label,
      p_true = p_true,
      alpha_at_ridge_near_truth = ridge_tbl$alpha[which.min(abs(ridge_tbl$p - p_true))],
      alpha_upper_at_rel_logpost_minus2 = width_2ll,
      curvature_alpha_at_truth = curvature_alpha,
      approx_halfwidth_alpha_for_rel_logpost_minus1 = sqrt(2 / curvature_alpha)
    )
  )
}

scenarios <- tibble(
  p_true = c(0.15, 0.25, 0.50, 0.75, 0.85),
  label = c(
    "True p = 0.15",
    "True p = 0.25",
    "True p = 0.50",
    "True p = 0.75",
    "True p = 0.85"
  )
)

results <- purrr::map2(
  scenarios$p_true,
  scenarios$label,
  function(p_true, label) {
    ii_obs <- simulate_obs(valid_obs_base$n_incident_draw, p_true)
    make_surface(ii_obs = ii_obs, p_true = p_true, label = label)
  }
)

surface_df <- purrr::map_dfr(results, "surface")
ridge_df <- purrr::map_dfr(results, "ridge")
alpha_slice_df <- purrr::map_dfr(results, "alpha_slice")
summary_df <- purrr::map_dfr(results, "summary")

surface_plot <- ggplot(surface_df, aes(x = p, y = alpha, fill = rel_log_post)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(
    aes(z = rel_log_post),
    color = "white",
    alpha = 0.35,
    bins = 8
  ) +
  geom_line(data = ridge_df, aes(x = p, y = alpha), color = "#b22222", linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5) +
  geom_point(
    data = scenarios %>% mutate(alpha = 1),
    aes(x = p_true, y = alpha),
    inherit.aes = FALSE,
    color = "black",
    size = 1.6
  ) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Detection probability (p)",
    y = "Global multiplier on latent EI",
    fill = "Relative\nlog posterior",
    title = "Observation-model geometry across simulated detection settings",
    subtitle = "Each panel uses the same latent EI baseline from one selected draw; black point marks the generating value (p_true, alpha = 1)"
  ) +
  theme_bw()

alpha_slice_plot <- ggplot(alpha_slice_df, aes(x = alpha, y = rel_log_post, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5) +
  labs(
    x = "alpha at fixed p = p_true",
    y = "Relative log posterior",
    color = NULL,
    title = "Cross-sections through alpha at the generating detection probability",
    subtitle = "Narrower curves indicate stronger local information about latent incident infections"
  ) +
  theme_bw()

ggsave(
  "detection_geometry_surfaces.png",
  surface_plot,
  width = 14,
  height = 6.8,
  dpi = 300
)

ggsave(
  "detection_geometry_cross_sections.png",
  alpha_slice_plot,
  width = 10,
  height = 5.6,
  dpi = 300
)

write_csv(summary_df, "detection_geometry_comparison_summary.csv")
print(summary_df)
