pacman::p_load(here, tidyverse, janitor, readxl, sf)

# sheets in data
loc <- here("scratch/south-coast/Cocrackerup_Peniup_RedMoort_cameradata.xlsx")
sheets <- excel_sheets(loc)

# site locations
locs_raw <- map(
  list.files(here("scratch/south-coast"), pattern = "\\.shp", full.names = T),
  ~ read_sf(.) |> clean_names()
) |>
  bind_rows()

locs <- locs_raw |>
  bind_cols(st_coordinates(locs_raw)) |>
  select(site = name, longitude = X, latitude = Y)

# deployment stuff
sites_raw <- read_excel(loc, sheet = sheets[[2]]) |>
  clean_names() |>
  rename_with(
    ~ paste0("early_death_", seq_along(.x)),
    .cols = starts_with("early_death")
  ) |>
  rename_with(
    ~ paste0("download_", seq_along(.x)),
    .cols = starts_with("download")
  ) |>
  mutate(across(where(is.POSIXct), ymd), site = waypoint)

sites <- sites_raw |>
  mutate(
    start = ymd(date_set_up),
    end = ymd(sites_raw[[ncol(sites_raw)]]),
    .after = waypoint
  ) |>
  left_join(locs, by = "site") |>
  mutate(site = factor(site, levels = str_sort(unique(site), numeric = T)))

failures <- sites |>
  pivot_longer(
    cols = matches("^(early_death_|download_)"),
    names_to = c(".value", "period"),
    names_pattern = "(early_death|download)_(\\d+)"
  ) |>
  drop_na(early_death)

# detections
dets_raw <- read_excel(loc, sheet = sheets[[1]]) |>
  clean_names()

dets <- dets_raw |>
  mutate(
    date = ymd(date),
    time = hms::as_hms(time),
    timestamp = ymd_hms(str_c(date, time)),
    site = factor(location, levels = levels(sites$site)),
    species = case_when(
      species == "Black-gloved wallaby" ~ "Brush wallaby",
      .default = species
    )
  ) |>
  filter(
    !str_detect(spp_category, "Bird|Reptile|unidentified|no ID|Unidentified"),
    !str_detect(species, "no ID|eagle|Honey")
  ) |>
  mutate(species = factor(str_to_title(species))) |>
  arrange(species, timestamp) |>
  select(site, timestamp, species, no_individuals) |>
  drop_na()

# setup; this was a pain
Sys.setenv(
  "RETICULATE_PYTHON" = "/Users/matt.hollanders/.virtualenvs/rgee/bin/python"
)
Sys.setenv(
  "EARTHENGINE_GCLOUD" = "/Users/matt.hollanders/google-cloud-sdk/bin/"
)
use_virtualenv("rgee", required = TRUE)
library(rgee)
ee <- import("ee")
ee$Authenticate()
ee$Initialize(project = "heroic-rain-367400")
ee_dir <- "~/.config/earthengine"
if (!dir.exists(ee_dir)) {
  dir.create(ee_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    data.frame(user = "matthijshollanders", drive_cre = NA, gcs_cre = NA),
    file = path.expand(str_c(ee_dir, "/rgee_sessioninfo.txt")),
    row.names = F
  )
}

deployments
observations
failures

observations_fox <- observations |>
  filter(scientificName == "Fox") |>
  droplevels()

stan_data <- make_data(deployments,
                       ,
                       deploymentID = site, days = 14, thin = 30)

fits[[2]] |>
  spread_rvars(kappa_O[s, ss]) |>
  ggplot(aes(xdist = kappa_O, y = factor(s))) +
  facet_wrap(~ ss) +
  geom_vline(xintercept = 0) +
  stat_pointinterval()

plot_correlations(fit, "site_effects")

observations$scientificName |> levels()

fit_model(stan_data, )

fit <- fit_model(stan_data)

plot_site_detection(fits[[2]])
plot_survey_detection(fits[[2]], species = c("Fox", "Malleefowl"))
fits[[2]]$loo()

map(fits[1:2], ~.$loo("log_lik2")) |>
  loo::loo_compare()

fit_model()

# prepare sites as EE object
deployments_to_fc <- function(
  deployments,
  site = site,
  longitude = "longitude",
  latitude = "latitude",
  crs = 4326
) {
  deployments |>
    st_as_sf(coords = c(longitude, latitude), crs = crs) |>
    select({{ site }}) |>
    sf_as_ee()
}

# terrain function
extract_terrain <- function(deployments, site = site) {
  dem <- ee$Image("USGS/SRTMGL1_003")
  info <- dem$rename("elevation")$addBands(ee$Terrain$slope(dem)$rename(
    "slope"
  ))$addBands(ee$Terrain$aspect(dem)$rename("aspect"))$reduceRegions(
    collection = deployments_to_fc(deployments),
    reducer = ee$Reducer$mean(),
    scale = 30
  )$getInfo()
  site_lvl <- pull(deployments, {{ site }}) |> levels()
  info$features |>
    map(~ as_tibble(.$properties)) |>
    list_rbind() |>
    relocate({{ site }}, .before = 1) |>
    mutate({{ site }} := factor({{ site }}, levels = site_lvl))
}
ter <- extract_terrain(sites)
# ter_scaled <- ter |>
#   mutate(aspect_radian = aspect * pi / 180,
#          northness = cos(aspect_radian),
#          eastness = sin(aspect_radian)) |>
#   select(site, elevation, slope, northness, eastness) |>
#   mutate(across(-site, ~scale(.)[, 1]))

encode_predictors(ter, "site", dep_id_chr = "site")

encode_predictors(
  ter |>
    mutate(
      veg = factor(
        sample(c("forest", "grassland"), 30, replace = T),
        levels = c("forest", "grassland")
      ),
      cool = factor(sample(c("uncool", "cool"), 30, replace = T), ordered = T)
    ),
  "site",
  dep_id_chr = "site",
  scale = T,
  days = 7
)

load_all()
test <- encode_predictors(
  silo,
  "survey",
  "site",
  "date",
  scale = F,
  survey_summary = list(maxt = mean, sqrt_rain = sum),
  reference_date = min(silo$date),
  days = 7
)
glimpse(test)

library(cropgrowdays)
silo_raw <- get_multi_silodata(
  latitude = sites$latitude,
  longitude = sites$longitude,
  Sitename = sites$site,
  email = "matthijs.hollanders@gmail.com",
  START = min(sites$start) |> str_remove_all("-"),
  FINISH = max(sites$end) |> str_remove_all("-")
)

# select some and convert to long format
silo <- silo_raw |>
  mutate(sqrt_rain = sqrt(rain)) |>
  select(site = Sitename, date = date_met, maxt, sqrt_rain, vp)

silo |> arrange(maxt)

deployments <- sites |>
  select(
    site = site,
    deploymentStart = start,
    deploymentEnd = end,
    longitude = longitude,
    latitude = latitude
  )

observations <- dets |>
  select(
    site = site,
    eventStart = timestamp,
    scientificName = species,
    count = no_individuals
  )

failures <- failures |>
  select(site = site, failureStart = early_death, failureEnd = download) |>
  mutate(site = as.character(site))
load_all()
dh <- make_data(
  deployments,
  observations |> filter(scientificName %in% c("Fox", "Malleefowl")) |> droplevels(),
  failures,
  deploymentID = site,
  thin = 30,
  days = 14,
  # occupancy_site_predictors = ter,
  # detection_site_predictors = ter,
  # survey_predictors = silo,
  # survey_summary = list(sqrt_rain = "sum")
)
dh$X_ord3 |> max()
glimpse(dh)

covs <- tibble(
  site = unique(deployments$site),
  temp = rnorm(30),
  rain = rexp(30),
  fire = sample(c("low", "medium", "high"), 30, replace = T) |>
    factor(levels = c("low", "medium", "high"), ordered = T),
  vegetation = sample(c("forest", "grassland", "wetland"), 30, replace = T) |>
    factor()
)
encode_site_predictors(covs, "deploymentID")

encode_site_predictors(
  covs |> mutate(deploymentID = 1:20),
  dep_id_chr = "deploymentID",
  scale = T
)

plot_coefficients(fit, component = "survey", unconditional = TRUE)

rowSums(dh$Delta)
glimpse(dh)
hist(dh$XY)
dh


load_all()
plot_kappa(fit, unconditional = TRUE)
fig <- plot_site_detection(
  fit,
  back_transform = T,
  map = T,
  intercepts = T,
  alpha = 0.6
)
attr(fig, "plot_data")

order

library(devtools)
document()
load_all()
plot_site_coefficients(fit)


test <- expand.grid(
  component = c("occupancy", "detection"),
  predictors = c("continuous", "categorical", "ordinal"),
  level = c("species", "mean"),
  facet_by = c("species", "predictor"),
  stringsAsFactors = F
)
for (i in 24:nrow(test)) {
  print(i)
  p <- plot_site_coefficients(
    fit,
    component = test$component[i],
    predictors = test$predictors[i],
    level = test$level[i],
    facet_by = test$facet_by[i]
  )
  print(attr(p, "plot_data"))
}

test <- plot_site_coefficients(fit)
attr(test, "plot_data")

fit$draws(c("psi_theta", "mu_theta")) |> mcmc_trace()


slice(-(1:5)) |>
  purrr::pwalk(\(component, predictors, level, facet_by) {
    plot_site_coefficients(
      fit,
      component = component,
      predictors = predictors,
      level = level,
      facet_by = facet_by
    )
  })

configs <- expand_grid(spatial = c("gp", "mvn", "none"),
                       temporal = c("gp", "mvn", "none"),
                       overdispersion = c("none", "olre", "nb"),
                       periodic_gp = c(TRUE, FALSE),
                       species_length_scales = c(TRUE, FALSE),
                       variance_decomposition = c("dirichlet", "logistic-normal"),
                       D = c(0, 10)
                       ) |>
  filter(!(periodic_gp == TRUE & temporal != "gp")) |>
  filter(!(D == 10 & (spatial == "none" & overdispersion != "olre"))) |>
  filter(!(species_length_scales == TRUE & (spatial != "gp" & temporal != "gp")))


configs <- expand_grid(spatial = c("gp", "mvn", "none"),
                       temporal = c("gp", "mvn", "none"),
                       overdispersion = c("none", "olre", "nb"),
                       species_length_scales = c(TRUE, FALSE),
) |>
  filter(!(species_length_scales == TRUE & (spatial != "gp" & temporal != "gp")))

fits <- map(1:nrow(configs), ~{
  print(.)
  fit_model(
    dh,
    iter_warmup = 200,
    iter_sampling = 200,
    spatial = configs$spatial[.],
    temporal = configs$temporal[.],
    overdispersion = configs$overdispersion[.],
    species_length_scales = configs$species_length_scales[.],
    latent = T,
    threads = 4,
    chains = 2,
    init = 0.1,
    grainsize = floor(dh$I / 4),
    max_treedepth = 10,
    adapt_delta = 0.8,
    show_exceptions = F)
})

fits2 <- pmap(configs[42,],
             \(spatial, temporal, overdispersion, species_length_scales) {
               fit_model(
                 dh,
                 iter_warmup = 500,
                 iter_sampling = 500,
                 spatial = spatial,
                 temporal = temporal,
                 overdispersion = overdispersion,
                 species_length_scales = species_length_scales,
                 latent = TRUE,
                 threads = 2,
                 chains = 4,
                 init = 0.1,
                 grainsize = floor(dh$I / 4),
                 max_treedepth = 10,
                 adapt_delta = 0.8,
                 show_exceptions = T
               )
             })



dh <- make_data(
  deployments,
  observations,
  failures,
  deploymentID = site,
  thin = 30,
  days = 14,
  # occupancy_site_predictors = ter,
  # detection_site_predictors = ter,
  # survey_predictors = silo,
  # survey_summary = list(sqrt_rain = "sum")
) |>
  glimpse()

fit_model(
  dh,
  iter_warmup = 200,
  iter_sampling = 200,
  spatial = "gp",
  temporal = "gp",
  periodic_gp = F,
  species_length_scales = T,
  overdispersion = "none",
  latent = F,
  threads = 4,
  chains = 2,
  init = 0.1,
  grainsize = floor(dh$I / 4),
  max_treedepth = 10,
  adapt_delta = 0.8
)

configs$spatial[1]

fit <- fit_model(
  dh,
  iter_warmup = 1000,
  iter_sampling = 1000,
  periodic_gp = T,
  species_length_scales = F,
  latent = F,
  overdispersion = "nb",
  threads = 4,
  chains = 2,
  grainsize = floor(dh$I / 4),
  max_treedepth = 10,
  adapt_delta = 0.8
)

plot_coefficients(fit)

int_mode(c(0, 0, 0, 0, 2))

fits <- map(
  c("none", "olre", "nb"),
  ~ fit_model(
    dh,
    iter_warmup = 200,
    iter_sampling = 200,
    periodic_gp = T,
    species_length_scales = F,
    latent = T,
    overdispersion = .,
    threads = 2,
    grainsize = floor(dh$I / 4),
    max_treedepth = 6
  )
)
loos <- map(fits, ~ .$loo("log_lik"))
loos2 <- map(fits, ~ .$loo("log_lik2"))

map(
  fits,
  ~ spread_rvars(., z[i, s]) |>
    summarise(Z = rvar_sum(z), .by = s)
) |>
  list_rbind(names_to = "mod") |>
  arrange(s) |>
  print(n = 30)


set_priors()

bayesplot::pp_check(
  apply(data$y, c(1, 3), sum) |> c(),
  yrep = fit$draws("Qrep", format = "draws_matrix"),
  group = rep(attr(data, "species"), each = data$I),
  fun = ppc_rootogram_grouped
)


failures <- tribble(
  ~deploymentID , ~failureStart     , ~failureEnd       ,
  "cam01"       , ymd("2026-01-01") , ymd("2026-01-06") ,
  "cam01"       , ymd("2026-02-03") , ymd("2026-02-07") ,
  "cam02"       , ymd("2026-01-04") , ymd("2026-01-07")
)
failures

plot_kappa(fit)

load_all()
p <- plot_site_detection(fit)
p
attr(p, "plot_data")


p <- plot_survey_detection(
  fit,
  point_interval = "median_qi",
  species = c("Fox", "Malleefowl"),
  surveys = attr(dh, "surveys")[20:60],
  linewidth = 0.5,
  show.legend = F,
  intercepts = T,
  back_transform = T,
) +
  facet_wrap(~species, ncol = 1)
ggsave("man/figures/kappa.png", width = 12, height = 7, dpi = 600)

attr(p, "plot_data")

?scale_fill_b

plot_iota(
  fit,
  sites = attr(dh, "sites")[1:20],
  intercepts = F,
  back_transform = F,
  map = F
)


attr(dh, "sites")[1:20]
attr(dh, "sites")
dh
glimpse(dh)
plot_iota(fit) +
  coord_fixed(xlim = c(655, 666), ylim = c(6213, 6218))
ggsave("man/figures/iota.png", width = 12, height = 6, dpi = 600)

plot_kappa(fit)

?set_priors
fit$draws("psi_bar") |> mcmc_trace()
fit$loo("log_lik2")

priors = set_priors()
print(priors)
priors

fit$stan_data
attr(dh, "stan_data")

library(cmdstanr)
library(here)
cmdstan_model(here("inst/stan/occARU.stan"))

devtools::document()
