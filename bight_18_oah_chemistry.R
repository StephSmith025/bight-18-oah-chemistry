library(dplyr)
library(tidyr)
library(units)

source("R/helpers.R") # seacarb helper functions

# idea is to use aragonite derived from bottle data to predict aragonite along entire CTD
# profiles, using methods as in Alin et al (2012)
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JC007511

# reference values from Alin et al (2012)
# the linear model centers variables to reduce bias in the interaction term
OFFSET_TEMP <- 10.28
OFFSET_OXY <- 138.46
OFFSET_SAL <- 33.889
OFFSET_SIGMA <- 26.01

# read in bottle data (raw data processed in R/chemistry.R) ####
# associate bottle samples with appropriate CTD samples
joined_bottle_ctd_data <- readRDS("data/bight_chemistry_bottle.rds") |>
  inner_join(
    readRDS("data/bight_chemistry_ctd.rds") |>
      select(Station_ID, Season, Date, Depth, Salinity, Temperature, Oxygen, Density),
    by = join_by(Station_ID, Date, Depth)
  ) |>
  filter(!is.na(Oxygen)) |>
  unnest(c(TAQACode, pHQACode))

# calculate insitu pH and OmegaAragonite
# center Temperature, Oxygen, Density, and Salinity (subtract reference values)
# calculate carbonate chemistry with seacarb
bight_model_data <- joined_bottle_ctd_data |>
  mutate(
    across(where(~ inherits(.x, "units")), as.numeric),
    pres = Depth / 10, # estimate of pressure in bar
    bottle_pH_lab = pH,
    pH = purrr::pmap_dbl(
      .l = list(pH, pHQACode, TA, Temperature, LabTemp, pres, Salinity),
      .f = calc_ph_insitu # seacarb helper, convert pH to insitu pH
    ),
    carb_bottle = purrr::pmap(
      .l = list(pH, TA, Salinity, Temperature, pres),
      .f = calc_carb_bottle # seacarb helper
    ),
    Temperature = Temperature - OFFSET_TEMP,
    Oxygen = Oxygen - OFFSET_OXY,
    Density = Density - OFFSET_SIGMA,
    Salinity = Salinity - OFFSET_SAL
  ) |>
  unnest(carb_bottle) |>
  mutate(
    DIC = set_units(DIC, "mol/kg") |> set_units("umol/kg") |> as.numeric(),
    TA = set_units(TA, "mol/kg") |> set_units("umol/kg") |> as.numeric()
  ) |>
  select(
    Station_ID, Season, Date, Depth, # metadata
    Temperature, Oxygen, Salinity, Density, # CTD
    bottle_Omega = OmegaAragonite, bottle_pH_insitu = pH, bottle_pH_lab, # bottle
    bottle_DIC = DIC, bottle_TA = TA # bottle
  )

# eliminate bad outlier for TA (way outside usual range for TA)
bight_model_data[bight_model_data$bottle_TA > 2300, "bottle_TA"] <- NA_real_

# fit models ####
models <- bight_model_data |>
  summarize(
    Omega = list(lm(bottle_Omega ~ Temperature * Oxygen)),
    pH = list(lm(bottle_pH_insitu ~ Temperature * Oxygen)),
    DIC = list(lm(bottle_DIC ~ Density * Oxygen)),
    TA = list(lm(bottle_TA ~ Temperature * Salinity))
  ) |>
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "model")

model_stats <- models |>
  mutate(
    formula = purrr::map_chr(model, function(x) {
      form <- formula(x$terms)
      paste(form[2], form[1], form[3])
    }),
    tidied = purrr::map(model, broom::glance)
  ) |>
  unnest(tidied) |>
  mutate(
    rmse = sigma * sqrt(df.residual / nobs) # convert residual standard error to rmse
  ) |>
  select(parameter, formula, n = nobs, adj.r.squared, rmse)

readr::write_csv(model_stats, file = "tables/alin_models_fitness_2025-12-11.csv")

model_coefs <- models |>
  mutate(coefficients = purrr::map(model, broom::tidy)) |>
  unnest(coefficients) |>
  select(-model)

readr::write_csv(model_coefs, file = "tables/alin_models_coefficients_2025-12-11.csv")


# estimate ctd profiles ####
all_ctd_predictions <- readRDS("data/bight_chemistry_ctd.rds") |>
  mutate(
    across(where(~ inherits(.x, "units")), as.numeric),
    Month = paste0(lubridate::year(Date), "-", sprintf("%02d", lubridate::month(Date)))
  ) |>
  rowwise() |>
  filter(!all(is.na(c(Temperature, Salinity, Density, pH, Oxygen)))) |> # a few rows of empty data
  ungroup() |>
  filter(!is.na(Oxygen)) |> # need oxygen for model predictions
  nest(data = -c(Station_ID, Season, Month, Date)) |>
  cross_join(models) |> # get one row per model per profile
  mutate(
    preds = purrr::map2(model, data, function(x, y) {
      centered_data <- y |>
        mutate(
          Temperature = Temperature - OFFSET_TEMP,
          Oxygen = Oxygen - OFFSET_OXY,
          Salinity = Salinity - OFFSET_SAL,
          Density = Density - OFFSET_SIGMA
        )
      predict(x, newdata = centered_data, interval = "confidence") |>
        as_tibble()
    })
  )

# predicted profiles data 2019-2020 ####
bight_alin_models_CTD_profiles <- all_ctd_predictions |>
  filter(between(Date, as.Date("2019-03-01"), as.Date("2020-03-01"))) |> # look at one full year
  select(-c(Season, Month, model)) |>
  unnest(c(data, preds)) |>
  select(-c(upr, lwr)) |>
  pivot_wider(
    names_from = "parameter",
    names_glue = "ctd_{parameter}_fit",
    values_from = "fit"
  ) |>
  left_join( # rejoin bottle derived values to final CTD data
    bight_model_data |>
      select(
        Station_ID, Date, Depth,
        bottle_Omega, bottle_TA, bottle_DIC, bottle_pH_insitu, bottle_pH_lab
      ),
    by = join_by(Station_ID, Date, Depth)
  ) |>
  select(
    Station_ID,
    Date,
    Depth_m = Depth,
    pH_ctd_orig = pH,
    Temperature_C = Temperature,
    Salinity,
    Oxygen_umol_kg = Oxygen,
    Density_kg_m3 = Density,
    predicted_OmegaAragonite = ctd_Omega_fit,
    predicted_pH = ctd_pH_fit,
    predicted_DIC_umol_kg = ctd_DIC_fit,
    predicted_TA_umol_kg = ctd_TA_fit,
    bottle_OmegaAragonite = bottle_Omega,
    bottle_pH_insitu,
    bottle_pH_lab,
    bottle_DIC_umol_kg = bottle_DIC,
    bottle_TA_umol_kg = bottle_TA
  )

bight_alin_models_CTD_profiles |>
  mutate(across(where(is.numeric), ~ round(.x, 8))) |>
  readr::write_csv(file = "tables/predicted_CTD_profiles_2025-12-11.csv", na = "")

# average predicted parameter values 15-100 m depth
bight_alin_models_CTD_profiles |>
  filter(Depth_m >= 15, Depth_m <= 100) |>
  tidyr::nest(profile = -c(Station_ID, Date)) |>
  mutate(
    avg = purrr::map(profile, function(x) {
      summarize(x, across(starts_with("predicted"), function(x) mean(x, na.rm = TRUE)))
    })
  ) |>
  tidyr::unnest_wider(avg)
