stan_data <- make_data(
  deployments,
  observations |>
    filter(scientificName %in% c("Fox", "Malleefowl")) |>
    droplevels(),
  survey_length = 14,
  deploymentID = site,
  detection_site_predictors = ter |>
    mutate(
      veg = sample(1:3, 30, replace = T) |>
        factor(labels = c("forest", "grass", "swamp")),
      water = sample(1:2, 30, replace = T) |>
        factor(labels = c("ocean", "river")),
      cloud = sample(1:4, 30, replace = T) |>
        factor(ordered = TRUE),
      rain = sample(1:3, 30, replace = T) |>
        factor(ordered = TRUE)
    ),
  survey_predictors = silo
)
fit <- fit_model(
  stan_data,
  periodic = TRUE,
  overdispersion = "olre",
  iter_warmup = 200,
  iter_sampling = 200,
  chains = 4,
  init = 0.1,
  threads = 2
)

plot_partitions(fit2)

plot_coefficients(fit, level = "mean")

?fit_model

plot_coefficients(fit, component = "survey")
plot_coefficients(
  fit,
  submodel = "detection",
  component = "survey",
  type = "ordinal"
)

plot_sites(fit)

attr(stan_data, "surveys")

plot_site_detection(fits[[2]], species = c("Fox"))

plot_survey_detection(
  fits[[2]],
  species = c("Fox", "Malleefowl"),
  surveys = attr(stan_data, "surveys")[30:100]
) +
  facet_wrap(~species, ncol = 1)

set_priors(psi_bar = c(3, 3))

fit_model(stan_data, iter_warmup = 500, adapt_delta = 0.95)

map(fits[1:2], ~ .$loo()) |>
  loo::loo_compare()


species
doesn
't show occupancy'


deployments
observations


?fit_model

setup_occARU()
