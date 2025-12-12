library(dplyr)
library(units)

# ingesting and cleaning raw bottle and CTD data

bight_18_bottle <- readxl::read_xlsx("data-raw/Bongo_Bottle_20251211.xlsx") |>
  mutate(
    across(where(is.numeric), function(x) if_else(x == -88 | x == -99, NA_real_, x)),
    Date = lubridate::as_date(Date),
    Station_ID = paste0(Agency, "_", StationID_D),
    `Lab Temp` = tidyr::replace_na(`Lab Temp`, 25)
  ) |>
  filter(pHQACode == "1" | TAQACode == "1", !is.na(`TA (umol/kg)`), !is.na(pH)) |>
  group_by(Station_ID, Date, Depth) |>
  summarize(
    TA = set_units(mean(`TA (umol/kg)`, na.rm = TRUE), "umol/kg") |> set_units("mol/kg"),
    pH = set_units(mean(pH, na.rm = TRUE), 1),
    TAQACode = list(unique(TAQACode)),
    pHQACode = list(unique(pHQACode)),
    `LabSalinity (ppt)` = set_units(mean(`LabSalinity (ppt)`, na.rm = TRUE), 1),
    `Lab Temp` = set_units(mean(`Lab Temp`, na.rm = TRUE), "C")
  ) |>
  mutate(
    Depth = set_units(Depth, "m")
  ) |>
  select(Station_ID, Date, Depth, TA, TAQACode, pH, pHQACode, LabSal = `LabSalinity (ppt)`, LabTemp = `Lab Temp`) |>
  ungroup()

saveRDS(bight_18_bottle, "data/bight_chemistry_bottle.rds")

bight_18_abc_4_missing <- readxl::read_xlsx("data-raw/B18_OA_CTD_ABCL_MissingData_12.2021_v2.xlsx") |>
  filter(lubridate::as_date(sampledate) == lubridate::as_date("2020-02-06"), station == "4-4") |>
  mutate(
    season = "Winter 2020",
    quarter = 4,
    agency = "ABC",
    sampledate,
    sampletime,
    station = 4,
    depth_m = depth,
    temperature_c = temperature,
    conductivity_s_m = `Conductivity (S/m)`,
    salinity_psu = salinity,
    density_kg_m3 = density,
    ph,
    oxygen_mg_L = `Oxygen (mg/L)`,
    .keep = "none"
  )

bight_18_ctd <- readxl::read_xlsx("data-raw/B18_CTD_20250129.xlsx") |>
  bind_rows(bight_18_abc_4_missing) |>
  mutate(
    across(where(is.numeric), function(x) if_else(x %in% c(-88, -1), NA_real_, x))
  ) |>
  mutate(
    Station_ID = paste0(agency, "_", station),
    Season = factor(season, levels = c(
      "Spring 2019", "Summer 2019", "Fall 2019", "Winter 2020", "Summer 2021",
      "Summer 2022", "Spring 2023", "Summer 2023"
    )),
    Date = lubridate::as_date(sampledate),
    Depth = set_units(depth_m, "m"),
    Temperature = set_units(temperature_c, "C"),
    Salinity = set_units(salinity_psu, 1),
    Density = set_units(density_kg_m3, "kg/m3"), # anomaly? probably
    pH = set_units(ph, 1),
    Oxygen = set_units(oxygen_mg_L / ((1000 + density_kg_m3) * 31.9988), "mol/kg") |> set_units("umol/kg"),
    .keep = "none"
  )

# use additional core stations for correction algorithm
all_stations_ctd <- readxl::read_xlsx(
  "data-raw/AllCoreStations_CTD_20251211.xlsx",
  col_types = c(
    "text", "date", "text", "text", "numeric", "numeric", "numeric",
    "numeric", "numeric", "numeric", "numeric", "text"
  ),
  na = c("", "NA")
) |>
  mutate(
    across(where(is.numeric), function(x) if_else(x %in% c(-88, -1), NA_real_, x))
  ) |>
  mutate(
    Station_ID = paste0(agency, "_", station),
    Season = "", # ignoring Season for the core stations
    Date = lubridate::as_date(sampledate),
    Depth = set_units(depth_m, "m"),
    Temperature = set_units(temperature_c, "C"),
    Salinity = set_units(salinity_psu, 1),
    Density = set_units(density_kg_m3, "kg/m3"),
    pH = set_units(ph, 1),
    Oxygen = set_units(oxygen_mg_L / ((1000 + density_kg_m3) * 31.9988), "mol/kg") |> set_units("umol/kg"),
    .keep = "none"
  )

bight_18_ctd |>
  bind_rows(all_stations_ctd) |>
  saveRDS("data/bight_chemistry_ctd.rds")
