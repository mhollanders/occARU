stan_data <- make_data(
  deployments,
  observations,
  failures = failures,
  survey_length = 14,
  detection_site_predictors = site_predictors,
  survey_predictors = survey_predictors,
  scale_predictors = T,
  survey_summary = list(sqrt_rain = sum)
)
load_all()
library(devtools)

attr(stan_data, "levels")

glimpse(stan_data)

par(mfrow = c(3, 1))
stan_data$X3 |> glimpse()

fit2 <- fit_model(
  stan_data,
  # stan_file = "inst/stan/occARU2.stan",
  periodic_gp = TRUE,
  overdispersion = "nb",
  iter_warmup = 200,
  iter_sampling = 200,
  chains = 4,
  threads = 2,
  pathfinder_args = list(
    max_lbfgs_iters = 200,
    refresh = 100,
    psis_resample = F
  ),
  show_exceptions = F,
  max_treedepth = 10
)
?pathfinder

tibble(X = stan_data$X3[, 1, ] |> c(), D = stan_data$Delta |> c()) |>
  filter(D == 0)

plot_coefficients(fit, component = "site", restricted = F)
install.packages("")

plot_coefficients(fit2, component = "survey", restricted = TRUE)
plot_coefficients(fit3, component = "survey", restricted = TRUE)

bayesplot::pp_check(
  apply(data$y, c(1, 3), sum) |> c(),
  yrep = fit$draws("Qrep", format = "draws_matrix")[1:50,],
  group = rep(attr(data, "species"), each = data$I),
  fun = "ppc_rootogram_grouped"
)
data <- stan_data

devtools::install_github(
  "stan-dev/bayesplot",
  dependencies = TRUE,
  build_vignettes = FALSE
)

stan_data$X3 |>
  reshape2::melt() |>
  as_tibble() |>
  ggplot(aes(ymd(Var3), value)) +
  facet_wrap(~Var2) +
  geom_point()

plot_coefficients(fit, component = "survey", restricted = F)

?fit_model

plot_coefficients(fit, component = "survey", restricted = F)
plot_surveys(fit7, restricted = F, include_predictors = F)


attr(fit, "species") <- str_c("Species ", c(1:9))

attr(fit, "occARU_data") |> attr("species") <- attr(stan_data, "species")
attr(fit, "occARU_data")$XY <- matrix(runif(60, 0, 10), 30, 2, dimnames = list(attr(stan_data, "sites"), c("X", "Y")))

plot_sites(fit, sites = attr(stan_data, "sites")[21:30],
           species = c("Fox", "Malleefowl", "Rabbit", "Tammar Wallaby"), alpha = 0.8) +
  scale_x_continuous(breaks = seq(1, 9, 2)) +
  scale_y_continuous(breaks = seq(1, 9, 2))  +
  facet_wrap(~ factor(species, labels = str_c("Species ", 1:4)), ncol = 2)
ggsave(here("man/figures/sites.png"), width = 8, height = 6)

plot_surveys(fit, species = c("Fox", "Malleefowl"), linewidth = 0.5, show.legend = F) +
  facet_wrap(~ factor(species, labels = str_c("Species ", 1:2)),
             ncol = 1) +
  scale_y_continuous(breaks = seq(0.2, 0.8, 0.2), expand = c(0, 0))
ggsave(here("man/figures/surveys.png"), width = 8, height = 6)

fit |>
  spread_rvars(kappa[s, j]) |>
  ggplot(aes(j, ydist = kappa)) +
  stat_lineribbon() + facet_wrap(~ s, ncol = 2)
fit$draws("iota_bar") |> mcmc_trace()

deployments <- sites |>
  select(
    deploymentID = site,
    deploymentStart = start,
    deploymentEnd = end,
    latitude,
    longitude
  )

observations <- dets |>
  select(
    deploymentID = site,
    eventStart = timestamp,
    scientificName = species,
    count = no_individuals
  )

failures <- failures |>
  select(deploymentID = site, failureStart = early_death, failureEnd = download)

site_predictors <- ter |>
  rename(deploymentID = site)

survey_predictors <- left_join(
  silo |>
    mutate(sqrt_rain = sqrt(rain)) |>
    select(deploymentID = Sitename, date = date_met, maxt, sqrt_rain, vp),
  vi |>
    rename(deploymentID = site) |>
    pivot_wider(names_from = index, values_from = value) |>
    complete(deploymentID, date = seq.Date(min(date), max(date))) |>
    fill(EVI, NDVI, .direction = "down") |>
    select(-NDVI)
) |>
  left_join(deployments |> select(contains("deployment"))) |>
  filter(between(date, deploymentStart, deploymentEnd), .by = deploymentID) |>
  select(-c(deploymentStart, deploymentEnd))

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

Good day,

I'm excited to announce the release of my first R package, occARU, which fits (multispecies) occupancy models for automated recording unit (ARU) data, such as camera traps and acoustic monitors, in Stan. After appropriately thinning species detections and aggregating to surveys of arbitrary lengths, occARU fits occupancy models with count observation models to go beyond binary detection probabilities. The package is somewhat opinionated, with most of the focus on the detection rates instead of occupancy.

A full description of the model can be found in the model vignette, but some of the main features are highlighted below:

* Hierarchical multispecies Gaussian processes (GPs) using matrix normals for site (spatial) and survey (temporal) random effects. If length scales are shared across species, only a single Cholesky decomposition is required for each GP.
* Orthogonal projection of the random effects ensures recovery of fixed effects.
* Fixed effects can include continuous, categorical (with zero-sum vectors), and ordinal (with simplex decomposition) predictors. The detection submodel accommodates both site and site-by-survey varying predictors.
* Global-local shrinkage priors are used for the occupancy and detection submodels. Inspired by R2D priors are set on the explained variance which are simplex decomposed using either Dirichlet or zero-sum logistic-normal distributions.
* Monte Carlo integration (thanks @avehtari) is used to produce a second `log_lik` object to be used for PSIS-LOO-CV. Log likelihoods are stored for each site and species combination, and random effects at this "observation-level" is known to be problematic for loo. The Monte Carlo approach isn't perfect, but it does reduce a lot of the high Pareto k observations.

I'd like to give special thanks to the Stan developers, with the underlying Stan programs showcasing some recent additions to the Stan langauge:

* `sum_to_zero_vector` and `sum_to_zero_matrix` for all random effects for identifiability;
* Ragged arrays of zero-sum vectors and simplexes using the `sum_to_zero_jacobian` and `simplex_jacobian` functions;
* Pathfinder to set initial values;
* Within-chain parallelisation, chunking over sites to increase speed,
* `bayesplot::ppc_rootogram_grouped()` will be useful for posterior predictive checking across species.
