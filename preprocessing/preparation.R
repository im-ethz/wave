## ---- message=F------------------------------------------------------------------------------------------
# 0am-5am measurements: HR and RF features are computed directly on the full measurements from 0am to 5am
features_full_period <- read_csv("../data/survival_morning.csv")

# 5min measurements: HRV features are computed in 5min intervals and then averaged over the time from 0am to 5am
features_intervals <- read_csv("../data/wave_features_5min.csv") %>%
  relocate(datetime, ID) %>%
  arrange(ID, datetime)


## --------------------------------------------------------------------------------------------------------
# negative SARS-CoV-2 result
exclude_negative <- c("wave-011")

# technical problems during the recording or nonadherence to the prescribed measurement regime
exclude_measurement <- c("wave-003", "wave-004", "wave-005", "wave-025")

# hospital discharge on the same day of hospitalization
exclude_short_stay <- c("wave-018")

exclude <- unique(c(exclude_negative, exclude_measurement, exclude_short_stay))

# patients eventually admitted to ICU
icu_patients <- c("wave-017", "wave-022", "wave-026", "wave-034", "wave-040", "wave-048", "wave-052")

# wave-020 switched hospitals, wave-037 dropped out
dropouts <- c("wave-020", "wave-037")


## ---- message=F------------------------------------------------------------------------------------------
meta_data <- read_csv("../data/subjects-metadata.csv") %>%
  mutate(start_day = floor_date(start, unit = "day"), .before = start) %>%
  mutate(end_day = floor_date(end, unit = "day"), .before = end) %>%
  mutate(ID = tolower(new_subject_id), .after = new_subject_id)

# data quality assurance
meta_data %<>%
  mutate(start_day = case_when( # in the case of wave-028 and wave-038, wearable measurements started after midnight, so the automatically inferred admission date must be corrected
    new_subject_id == "WAVE-028" ~ as.POSIXct("2021-01-03", tz = "UTC"),
    new_subject_id == "WAVE-038" ~ as.POSIXct("2021-02-16", tz = "UTC"),
    T ~ start_day
  ))

meta_data %<>% mutate(
  total_days = as.integer(difftime(end_day, start_day, units = "days")),
  .after = new_subject_id
) %>%
  mutate(total_days = case_when( # wave-17 was only followed for 1 day
    new_subject_id == "WAVE-017" ~ 1L,
    T ~ total_days
  ))


## ---- message=F------------------------------------------------------------------------------------------
demographics <- read_csv("../data/patient_demographics.csv") %>%
  transmute(ID = str_to_lower(id), sex = factor(sex), age = as.integer(age), age60 = age >= 60)


## --------------------------------------------------------------------------------------------------------
covariates <- names(features_intervals %>% select(-c(datetime, ID)) %>% select(!contains(c("sign_changes", "num_ibis"))))

# 5min measurements
surv_data_intervals <- features_intervals %>%
  inner_join(meta_data %>% transmute(ID = tolower(new_subject_id), total_days, start_day, end_day), by = "ID") %>%
  group_by(ID) %>%
  mutate(
    icu = ID %in% icu_patients,
    rel_day = (as.integer(difftime(floor_date(datetime, unit = "day"), start_day, units = "days"))),
    rem_sync_days = (as.integer(difftime(end_day, floor_date(datetime, unit = "day"), units = "days"))),
    time_of_day = as_hms(floor_date(datetime, unit = "5minutes")),
    .after = datetime
  ) %>%
  filter(time_of_day < as_hms("5:00:00")) %>%
  group_by(ID, icu, rel_day, rem_sync_days) %>%
  summarize(
    datetime = min(datetime, na.rm = T),
    hrv_n_windows_total = n(),
    hrv_n_windows_notna = sum(!is.na(hrv_mean_nni)),
    across(contains(c("sign_changes", "num_ibis", "n_obs", "n_above_mean", "n_below_mean")), sum, na.rm = T),
    across(all_of(covariates), mean, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(ID = as.factor(ID), rem_sync_days = as.integer(rem_sync_days)) %>%
  arrange(ID, desc(rem_sync_days)) %>%
  ungroup() %>%
  mutate(event = case_when(
    rem_sync_days > 0 ~ "stay",
    icu ~ "icu",
    T ~ "dismiss"
  ), .after = ID) %>%
  inner_join(demographics, by = "ID")

# 0am-5am measurements
surv_data_full_period <- features_full_period %>%
  mutate(
    ID = as.factor(ID),
    event = case_when(
      event == "stays" ~ "stay",
      event == "dismissed" ~ "dismiss",
      T ~ event
    )
  ) %>%
  inner_join(demographics, by = "ID")

# make dropouts right-censored
surv_data_intervals %<>%
  mutate(event = ifelse((ID %in% dropouts), "stay", event))
surv_data_full_period %<>%
  mutate(event = ifelse((ID %in% dropouts), "stay", event))

# make event a factor
surv_data_intervals %<>% mutate(event = factor(event, ordered = T, levels = c("dismiss", "stay", "icu")))
surv_data_full_period %<>% mutate(event = factor(event, ordered = T, levels = c("dismiss", "stay", "icu")))


## --------------------------------------------------------------------------------------------------------
surv_data_full_period %<>% filter(!(ID %in% exclude))
surv_data_intervals %<>% filter(!(ID %in% exclude))


## --------------------------------------------------------------------------------------------------------
quality <- function(df) {
  df %>%
    filter(!(ID == "wave-021" & rel_day == 6)) %>% # night 6 has not enough data
    filter(!(ID == "wave-028" & rel_day == 1)) %>% # night 1 has not enough data
    filter(!(ID == "wave-049" & rel_day == 4)) # night 4 has not enough data
}

surv_data_intervals %<>% quality()
surv_data_full_period %<>% quality()


## --------------------------------------------------------------------------------------------------------
surv_data_intervals %<>% relocate(ID, rel_day, rem_sync_days, event, icu, age, age60, sex) %>%
  arrange(ID, rel_day)

surv_data_full_period %<>% relocate(ID, datetime, rel_day, event, age, age60, sex) %>%
  arrange(ID, rel_day)


## --------------------------------------------------------------------------------------------------------
covariates <- c(names(features_intervals %>% select(-c(datetime, ID)) %>% select(!contains(c("sign_changes", "num_ibis", "n_obs")))))
surv_data_intervals %<>%
  mutate(age = scale(age), across(all_of(covariates), scale))

covariates <- names(surv_data_full_period %>% select(-(ID:sex)) %>% select(!contains(c("sign_changes", "num_ibis", "n_obs"))))
surv_data_full_period %<>%
  mutate(age = scale(age), across(all_of(covariates), scale))


## --------------------------------------------------------------------------------------------------------
sanity_check <- surv_data_intervals %>%
  select(ID, rel_day, event) %>%
  full_join(surv_data_full_period %>% select(ID, rel_day, event, datetime), by = c("ID", "rel_day")) %>%
  filter(event.x != event.y | is.na(event.y) | is.na(event.y))

if (nrow(sanity_check) > 0) {
  warning("Not all observations were matched.")
  print(sanity_check)
}


## --------------------------------------------------------------------------------------------------------
surv_data <- surv_data_full_period %>%
  select(-starts_with("hrv_")) %>%
  full_join(surv_data_intervals %>% select(ID:rel_day, starts_with("hrv_")), by = c("ID", "rel_day"))

# check
nrow(surv_data) == nrow(surv_data_full_period) & nrow(surv_data) == nrow(surv_data_intervals)


## --------------------------------------------------------------------------------------------------------
coverage_factor <- 0.5

observation_coverage <- surv_data %>%
  filter(!(ID %in% exclude), !is.na(HEART_RATE_mean), !is.na(RESPIRATION_mean), !is.na(hrv_mean_nni)) %>%
  transmute(ID, rel_day, event, HEART_RATE_n_obs, RESPIRATION_n_obs, hrv_n_windows_notna) %>%
  mutate(HEART_RATE_n_obs = HEART_RATE_n_obs / 18000, RESPIRATION_n_obs = RESPIRATION_n_obs / 18000, hrv_n_windows_notna = hrv_n_windows_notna / 60)
