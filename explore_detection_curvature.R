library(tidyverse)
library(cmdstanr)

rebuild_utah_data <- function() {
  case <- read_csv("./data/time_series_covid19_confirmed_US.csv", show_col_types = FALSE) %>%
    filter(`Province_State` == "Utah") %>%
    dplyr::select(-c(1:5, 7:11)) %>%
    pivot_longer(-Admin2, names_to = "date", values_to = "cases") %>%
    mutate(date = lubridate::mdy(date)) %>%
    group_by(Admin2, date) %>%
    summarize(cases = sum(as.numeric(cases)), .groups = "drop") %>%
    mutate(date = 7 * (as.numeric(date - min(date)) %/% 7) + min(date)) %>%
    group_by(date, Admin2) %>%
    summarize(cases = sum(cases), .groups = "drop") %>%
    filter(date < lubridate::ymd("2021-06-02")) %>%
    group_by(Admin2) %>%
    mutate(new_cases = cases - lag(cases)) %>%
    arrange(Admin2, date) %>%
    filter(!is.na(new_cases)) %>%
    ungroup()

  dat_final <- case %>% filter(Admin2 %in% readRDS("final_counties.rds"))
  pop <- read_csv("./data/utah_counties_pop_coord.csv", show_col_types = FALSE) %>%
    arrange(desc(Population_2020))

  d1 <- dat_final %>%
    filter(date < lubridate::ymd("2020/8/19")) %>%
    rename(County = "Admin2") %>%
    left_join(pop, by = "County") %>%
    arrange(desc(Population_2020), date) %>%
    dplyr::select(-Population_2020, -Latitude, -Longitude, -cases) %>%
    pivot_wider(names_from = County, values_from = new_cases) %>%
    dplyr::select(-date)

  d1[6, 9] <- 1

  counties <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12)
  ii <- as.matrix(d1)[, counties]
  TT <- nrow(ii) + 1
  N_C <- length(counties)
  pop_size <- pop$Population_2020[counties]

  first <- apply(matrix(as.numeric(ii != 0), ncol = N_C), 2, function(x) which(x == 1)[1])
  last <- purrr::map_int(seq_len(N_C), function(col) {
    w <- which(ii[, col] == 0)
    w <- w[w > first[col]][1] - 1
    ifelse(is.na(w), TT - 1, w)
  })

  county_lookup <- tibble(
    county = colnames(ii),
    county_index = seq_len(N_C)
  )

  obs_tbl <- as_tibble(ii, .name_repair = "minimal") %>%
    mutate(time = row_number()) %>%
    pivot_longer(-time, names_to = "county", values_to = "ii_obs") %>%
    left_join(county_lookup, by = "county")

  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size,
    first = first,
    last = last,
    obs_tbl = obs_tbl
  )
}

extract_draw_ei <- function(fit, draw_id, TT, N_C) {
  draw_vec <- fit$draws("ei_t", format = "draws_matrix")[draw_id, ]
  tibble(
    variable = colnames(fit$draws("ei_t", format = "draws_matrix")),
    ei_prop = as.numeric(draw_vec)
  ) %>%
    tidyr::extract(
      variable,
      into = c("time", "county_index"),
      regex = "ei_t\\[(\\d+),(\\d+)\\]",
      convert = TRUE
    ) %>%
    filter(time <= TT - 1, county_index <= N_C)
}

loglik_obs <- function(p, ii_obs, n_incident) {
  if (p <= 0 || p >= 1) {
    return(-Inf)
  }

  sd_obs <- sqrt(pmax(n_incident * p * (1 - p), 1e-10))
  sum(dnorm(ii_obs, mean = p * n_incident, sd = sd_obs, log = TRUE))
}

if (sys.nframe() == 0) {
  fit <- as_cmdstan_fit(list.files("tst_folder_sig", full.names = TRUE))
  draw_summary <- fit$summary(c("p", "lp__"))
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

  valid_obs <- draw_dat$obs_tbl %>%
    left_join(
      tibble(
        county_index = seq_len(draw_dat$N_C),
        first = draw_dat$first,
        last = draw_dat$last
      ),
      by = "county_index"
    ) %>%
    filter(time >= first, time <= last) %>%
    left_join(ei_draw_tbl, by = c("time", "county_index"))

  p_draw <- fit$draws("p", format = "draws_matrix")[draw_id, 1]

  p_grid <- tibble(p = seq(0.02, 0.95, by = 0.0025)) %>%
    mutate(
      log_post = purrr::map_dbl(
        p,
        ~ loglik_obs(.x, ii_obs = valid_obs$ii_obs, n_incident = valid_obs$n_incident_draw)
      ),
      rel_log_post = log_post - max(log_post),
      density = exp(rel_log_post),
      panel = "Conditional posterior slice in p"
    )

  alpha_grid <- expand_grid(
    p = seq(0.02, 0.95, by = 0.01),
    alpha = seq(0.35, 2.5, by = 0.02)
  ) %>%
    mutate(
      log_post = purrr::pmap_dbl(
        list(p, alpha),
        function(p, alpha) {
          loglik_obs(
            p = p,
            ii_obs = valid_obs$ii_obs,
            n_incident = alpha * valid_obs$n_incident_draw
          )
        }
      ),
      rel_log_post = log_post - max(log_post)
    )

  ridge_tbl <- alpha_grid %>%
    group_by(p) %>%
    slice_max(order_by = log_post, n = 1, with_ties = FALSE) %>%
    ungroup()

  mode_row <- p_grid %>%
    slice_max(order_by = log_post, n = 1, with_ties = FALSE)

  slice_plot <- ggplot(p_grid, aes(x = p, y = rel_log_post)) +
    geom_line(linewidth = 0.8, color = "#1f4e79") +
    geom_vline(xintercept = p_draw, linetype = 2, color = "#b22222") +
    geom_vline(xintercept = mode_row$p, linetype = 3, color = "#555555") +
    annotate(
      "text",
      x = p_draw,
      y = min(p_grid$rel_log_post) * 0.1,
      label = paste0("draw p = ", round(p_draw, 3)),
      hjust = -0.05,
      vjust = 1,
      size = 3
    ) +
    labs(
      x = "Detection probability (p)",
      y = "Relative log posterior",
      title = "Conditional posterior slice with all other parameters fixed",
      subtitle = "One posterior draw from the Utah fit; posterior is proportional to the full observation likelihood because p has a uniform prior on (0, 1)"
    ) +
    theme_bw()

  surface_plot <- ggplot(alpha_grid, aes(x = p, y = alpha, fill = rel_log_post)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(
      data = alpha_grid,
      aes(x = p, y = alpha, z = rel_log_post),
      inherit.aes = FALSE,
      color = "white",
      alpha = 0.45,
      bins = 10
    ) +
    geom_hline(yintercept = 1, linetype = 2, color = "black", linewidth = 0.6) +
    geom_line(data = ridge_tbl, aes(x = p, y = alpha), inherit.aes = FALSE, color = "#b22222", linewidth = 0.8) +
    geom_point(
      data = tibble(p = p_draw, alpha = 1),
      aes(x = p, y = alpha),
      inherit.aes = FALSE,
      color = "black",
      size = 2.2
    ) +
    annotate(
      "text",
      x = min(alpha_grid$p) + 0.02,
      y = 1.03,
      label = "alpha = 1 (selected draw EI)",
      hjust = 0,
      vjust = -0.2,
      size = 3
    ) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    labs(
      x = "Detection probability (p)",
      y = "Global multiplier on latent EI",
      fill = "Relative\nlog posterior",
      title = "Observation-model ridge over p and latent incidence scale",
      subtitle = "alpha rescales all latent incident infections from the chosen draw; alpha = 1 is the unmodified EI trajectory from the selected posterior draw"
    ) +
    theme_bw()

  ggsave(
    "detection_p_slice.png",
    plot = slice_plot,
    width = 8,
    height = 4.6,
    dpi = 300
  )

  ggsave(
    "detection_p_alpha_surface.png",
    plot = surface_plot,
    width = 8,
    height = 5.2,
    dpi = 300
  )

  summary_tbl <- data.frame(
    draw_id = as.integer(draw_id),
    draw_p = as.numeric(p_draw),
    conditional_mode_p = as.numeric(mode_row$p[[1]]),
    n_observations = as.integer(nrow(valid_obs)),
    max_rel_log_drop_away_from_mode = as.numeric(min(p_grid$rel_log_post))
  )

  write_csv(summary_tbl, "detection_curvature_summary.csv")

  print(summary_tbl)
}
