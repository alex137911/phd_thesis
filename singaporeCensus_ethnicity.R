library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(MASS)
library(sandwich)
library(lmtest)
library(ggplot2)

# Load data
raw <- read_csv("~/Library/CloudStorage/OneDrive-Personal/Documents/McGill Human Genetics/Thesis/Paper 1/SingaporeResidentsBySingleYearOfAgeEthnicGroupAndSexAtEndJuneAnnual.csv")

# Year columns
year_cols <- names(raw)[str_detect(names(raw), "^\\d{4}$")]

# Helper: identify age rows
is_age_row <- function(x) {
  str_detect(x, "^\\d+ Year$|^\\d+ Years & Over$")
}

long <- raw %>%
  mutate(
    row_id = row_number(),
    DataSeries = str_squish(DataSeries),
    
    # age rows are things like "0 Year", "45 Year", "85 Years & Over"
    age_row = is_age_row(DataSeries),
    
    # parent label is only on the non-age rows
    parent_series = if_else(!age_row, DataSeries, NA_character_)
  ) %>%
  fill(parent_series, .direction = "down") %>%
  mutate(
    # optional block id for each repeated section
    block_id = cumsum(!age_row),
    
    # keep the row-level label too
    age_group = case_when(
      age_row ~ DataSeries,
      TRUE ~ "All ages"
    ),
    
    # only assign numeric age for true single-year rows like "0 Year", "1 Year", ...
    age = case_when(
      str_detect(age_group, "^\\d+ Year$") ~ parse_number(age_group),
      TRUE ~ NA_real_
    ),
    
    # optional metadata parsed from the parent label
    sex = case_when(
      str_detect(parent_series, "Male") ~ "Male",
      str_detect(parent_series, "Female") ~ "Female",
      TRUE ~ "Total"
    ),
    ethnicity = case_when(
      str_detect(parent_series, "Chinese") ~ "Chinese",
      str_detect(parent_series, "Malay") ~ "Malay",
      str_detect(parent_series, "Indian") ~ "Indian",
      str_detect(parent_series, "Other Ethnic Groups") ~ "Other Ethnic Groups",
      str_detect(parent_series, "Residents") ~ "All ethnicities",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "year",
    values_to = "count",
    values_transform = list(count = as.character)
  ) %>%
  mutate(
    year = as.integer(year),
    count = parse_number(count)
  ) %>%
  dplyr::select(
    row_id, block_id, parent_series, sex, ethnicity,
    age_group, age, year, count
  ) %>%
  arrange(block_id, row_id, year)

# Check the "parent_series" matches the levels of "raw$DataSeries"
long %>%
  filter(age_group == "0 Year") %>%
  distinct(parent_series)

# Filter to keep only non-sex stratified rows
analysis_df <- long %>%
  filter(age_group != "All ages", sex == "Total") %>%
  dplyr::select(year, age_group, age, ethnicity, count)

# Single-year analysis_df (Filter to keep only single-year ages - i.e., no age bins)
singleAnalysis_df <- analysis_df %>%
  dplyr::filter(
    str_detect(age_group, "^\\d+ Year$")   # keep only single-year ages
  )

# Breakdown Chinese ethnicity by age decile -----------------------------------------
ethnicity_breakdown <- singleAnalysis_df %>%
  filter(
    year %in% c(2005, 2025),
    ethnicity %in% c("Chinese", "Malay", "Indian", "Other Ethnic Groups")
  ) %>%
  group_by(year, ethnicity) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total_4groups = sum(count),
    percent = 100 * count / total_4groups
  ) %>%
  ungroup()

chinese_age_deciles <- singleAnalysis_df %>%
  filter(
    ethnicity == "Chinese",
    year %in% c(2005, 2025)
  ) %>%
  mutate(
    age_decile = paste0(
      floor(age / 10) * 10, "-",
      floor(age / 10) * 10 + 9
    )
  ) %>%
  group_by(year, age_decile) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total_chinese = sum(count),
    percent = 100 * count / total_chinese
  ) %>%
  ungroup()

# PLOTTING 
ethnicity_comp <- singleAnalysis_df %>%
  filter(
    year %in% c(2005, 2025),
    ethnicity %in% c("Chinese", "Malay", "Indian", "Other Ethnic Groups")
  ) %>%
  group_by(year, ethnicity) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    total = sum(count),
    percent = 100 * count / total
  ) %>%
  ungroup()

ggplot(ethnicity_comp, aes(x = factor(year), y = percent, fill = ethnicity)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = NULL,
    y = "Percent of population",
    fill = "Ethnicity",
    title = "Overall ethnic composition is relatively stable"
  ) +
  theme_minimal(base_size = 12)

age_levels <- c("0-9", "10-19", "20-29", "30-39", "40-49",
                "50-59", "60-69", "70-79", "80-89")

chinese_age_deciles_plot <- chinese_age_deciles %>%
  mutate(age_decile = factor(age_decile, levels = age_levels))

chinese_wide <- chinese_age_deciles %>%
  mutate(age_decile = factor(age_decile, levels = age_levels)) %>%
  dplyr::select(year, age_decile, percent) %>%
  pivot_wider(names_from = year, values_from = percent)

ggplot(chinese_wide, aes(y = age_decile)) +
  geom_segment(
    aes(x = `2005`, xend = `2025`, yend = age_decile),
    linewidth = 1, color = "grey70"
  ) +
  geom_point(aes(x = `2005`, color = "2005"), size = 3) +
  geom_point(aes(x = `2025`, color = "2025"), size = 3) +
  scale_color_manual(
    name = "Year",
    values = c("2005" = "#1f78b4", "2025" = "#e31a1c")
  ) +
  labs(
    x = "Percent of Chinese population",
    y = "Age decile",
    title = "Shift in the age distribution of people identifying as Chinese, 2005 to 2025",
    subtitle = "Values farther right in 2025 indicate an older age structure"
  ) +
  theme_minimal(base_size = 12)


plot<-ggplot(chinese_wide, aes(y = age_decile)) +
  geom_segment(
    aes(x = `2005`, xend = `2025`, yend = age_decile),
    linewidth = 0.5, color = "grey20"
  ) +
  geom_point(aes(x = `2005`, color = "2005"), size = 3) +
  geom_point(aes(x = `2025`, color = "2025"), size = 3) +
  scale_color_manual(
    values = c("2005" = "#0DAC99", "2025" = "#8D68B4")
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  coord_cartesian(xlim = c(0, 20)) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "Chinese Dumbell Plot.png",
  plot     = plot,
  width    = 11.88,
  height   = 10.3,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)

# -----------------------------
# Calculate percentage Chinese responses
pct_chinese_df <- analysis_df %>%
  dplyr::filter(
    ethnicity %in% c("All ethnicities", "Chinese"),
    str_detect(age_group, "^\\d+ Year$")   # keep only single-year ages
  ) %>%
  mutate(
    count = parse_number(as.character(count))
  ) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .by = c(year, age_group, age, ethnicity)
  ) %>%
  pivot_wider(
    names_from = ethnicity,
    values_from = count
  ) %>%
  mutate(
    pct_chinese = 100 * .data[["Chinese"]] / .data[["All ethnicities"]],
    prop_chinese = .data[["Chinese"]] / .data[["All ethnicities"]],
    cohort = year - age
  ) %>%
  arrange(year, age)

# Check which ages are not available
age_availability <- pct_chinese_df %>%
  summarise(
    n_years_observed = sum(!is.na(pct_chinese)),
    .by = age
  ) %>%
  arrange(age)

# Restrict to ages observed across all years
n_years_total <- n_distinct(pct_chinese_df$year)

complete_ages <- pct_chinese_df %>%
  summarise(
    n_years_observed = sum(!is.na(pct_chinese)),
    .by = age
  ) %>%
  filter(n_years_observed == n_years_total) %>%
  pull(age)

pct_chinese_complete <- pct_chinese_df %>%
  filter(age %in% complete_ages)

# # Plot heatmap
# ggplot(pct_chinese_complete, aes(x = year, y = age, fill = pct_chinese)) +
#   geom_tile() +
#   scale_y_reverse() +
#   labs(
#     x = "Year",
#     y = "Age",
#     fill = "% Chinese",
#     title = "Share of Singapore residents classified as Chinese"
#   ) +
#   theme_minimal()

#-------------------------------------------------------------------------------
# Make sure cohort exists
pct_chinese_complete <- pct_chinese_complete %>%
  mutate(cohort = year - age)

# Age profile
age_profile <- pct_chinese_complete %>%
  group_by(age) %>%
  summarise(
    pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
    population = sum(`All ethnicities`, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# Period profile
period_profile <- pct_chinese_complete %>%
  group_by(year) %>%
  summarise(
    pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
    population = sum(`All ethnicities`, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# Cohort profile
cohort_profile <- pct_chinese_complete %>%
  group_by(cohort) %>%
  summarise(
    pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
    population = sum(`All ethnicities`, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
   ) #%>%
  # # Optional: remove edge cohorts with very few age-period cells
  # filter(n_cells >= 5)

# Plot age
ggplot(age_profile, aes(x = age, y = pct_chinese)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    title = "Age profile of residents classified as Chinese",
    x = "Age",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

# Plot period
ggplot(period_profile, aes(x = year, y = pct_chinese)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    title = "Period profile of residents classified as Chinese",
    x = "Year",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

# Plot cohort
ggplot(cohort_profile, aes(x = cohort, y = pct_chinese)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    title = "Cohort profile of residents classified as Chinese",
    x = "Birth cohort (year of birth ≈ year - age)",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

# Plot cohort with smoothing (LOESS/local regression)
ggplot(cohort_profile, aes(x = cohort, y = pct_chinese)) +
  geom_line(linewidth = 0.5, color = "grey70") +
  geom_point(size = 1.2, color = "grey50") +
  geom_smooth(
    method = "loess",
    # Increase span to smoothen further
    span = 0.4,
    se = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  labs(
    title = "Cohort profile of residents classified as Chinese",
    x = "Birth cohort (year of birth ≈ year - age)",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

# Smoothing (LOESS) weighted by population (some cohorts represent much 
# larger populations than others)
ggplot(cohort_profile, aes(x = cohort, y = pct_chinese, weight = population)) +
  geom_line(linewidth = 0.5, color = "grey70") +
  geom_point(size = 1.2, color = "grey50") +
  geom_smooth(
    method = "loess",
    span = 0.4,
    se = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  labs(
    title = "Cohort profile of residents classified as Chinese",
    x = "Birth cohort (year of birth ≈ year - age)",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)


# Raw age-specific values for each cohort in the background
cohort_raw <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort
  ) %>%
  filter(!is.na(pct_chinese))

# Weighted cohort summary
cohort_profile <- cohort_raw %>%
  group_by(cohort) %>%
  summarise(
    pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
    population = sum(`All ethnicities`, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

plot<-ggplot() +
  geom_point(
    data = cohort_raw,
    aes(x = cohort, y = pct_chinese),
    color = "#E08DBA",
    alpha = 0.15,
    size = 1
  ) +
  geom_line(
    data = cohort_profile,
    aes(x = cohort, y = pct_chinese),
    linewidth = 1.0,
    lineend = "round",
    color = "#6C1D47"
  ) +
  scale_x_continuous(
    limits = c(1905, 2025),
    breaks = seq(1905, 2025, by = 10),
    expand = expansion(mult = c(0.041, 0.02))
  ) +
  scale_y_continuous(
    limits = c(60,91),
    breaks = seq(60, 90, by = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  # geom_point(
  #   data = cohort_profile,
  #   aes(x = cohort, y = pct_chinese),
  #   size = 1.2,
  #   color = "#6C1D47"
  # ) +
  # geom_smooth(
  #   data = cohort_profile,
  #   aes(x = cohort, y = pct_chinese, weight = population),
  #   method = "loess",
  #   span = 0.4,
  #   se = FALSE,
  #   linewidth = 1,
  #   color = "black"
  # ) +
  # labs(
  #   title = "Cohort profile of residents classified as Chinese",
  #   x = "Birth cohort (year of birth ≈ year - age)",
  #   y = "% Chinese"
  # ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  filename = "Singapore Chinese Responses by Cohort.png",
  plot     = plot,
  width    = 26.16,
  height   = 14.28,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)

# Annotation plot for 1977 cohort -------------
cohort1977_raw     <- subset(cohort_raw, cohort_raw$cohort == 1977)
cohort1977_profile <- subset(cohort_profile, cohort_profile$cohort == 1977)

ggplot() +
  geom_point(
    data = cohort1977_raw,
    aes(x = cohort, y = pct_chinese),
    color = "#E08DBA",
    alpha = 0.3,
    size = 1
  ) +
  geom_point(
    data = cohort1977_profile,
    aes(x = cohort, y = pct_chinese),
    color = "#6C1D47",
    size = 1.5
  ) +
  scale_y_continuous(limits = c(72,80),
                     breaks = seq(72, 80, by = 2),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = 1977, labels = "1977") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = "Singapore Chinese Annotated Weighted Mean.png",
  plot     = plot,
  width    = 3.6,
  height   = 5.56,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)

# Annotation plot for paired 1977 and 1978 cohort --------
cohort_pair_raw <- cohort_raw %>%
  filter(cohort %in% c(1977, 1978))

cohort_pair_profile <- cohort_profile %>%
  filter(cohort %in% c(1977, 1978))

link_df <- cohort_pair_raw %>%
  dplyr::select(cohort, age, pct_chinese) %>%
  pivot_wider(
    names_from = cohort,
    values_from = pct_chinese,
    names_prefix = "cohort_"
  ) %>%
  filter(!is.na(cohort_1977), !is.na(cohort_1978))

# Adjust so lines don't go under/to the centre of the points
# line_pad <- 0.005
line_pad <- 0.04

plot<-ggplot() +
  geom_segment(
    data = link_df,
    aes(
      x = 1977 + line_pad, xend = 1978 - line_pad,
      y = cohort_1977, yend = cohort_1978
    ),
    color = "grey80",
    alpha = 0.3,
    linewidth = 0.35,
    lineend = "round"
  ) +
  geom_point(
    data = cohort_pair_raw,
    aes(x = cohort, y = pct_chinese),
    color = "#E08DBA",
    alpha = 0.3,
    size = 1
  ) +
  geom_point(
    data = cohort_pair_profile,
    aes(x = cohort, y = pct_chinese),
    color = "#6C1D47",
    size = 1.5
  ) +
  scale_y_continuous(
    limits = c(71, 81),
    breaks = seq(71, 81, by = 2),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_continuous(
    breaks = c(1977, 1978),
    labels = c("1977", "1978"),
    expand = expansion(mult = c(1, 1))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = "Singapore Chinese Annotated Weighted Mean Paired.png",
  plot     = plot,
  width    = 3.6,
  height   = 5.56,
  units    = "cm",
  dpi      = 600,
  bg       = "white"
)

# % Chinese by age across cohorts ----------------------------------------------
required_ages <- 20:60

# # Keep only cohorts with complete age coverage
# cohort_complete_raw <- cohort_raw %>%
#   filter(age %in% required_ages) %>%
#   group_by(cohort) %>%
#   filter(all(required_ages %in% age)) %>%
#   ungroup()
# 
# # Optional: cohort-level profile restricted to ages 20-60
# cohort_complete_profile <- cohort_complete_raw %>%
#   group_by(cohort) %>%
#   summarise(
#     pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
#     population = sum(`All ethnicities`, na.rm = TRUE),
#     n_cells = n(),
#     .groups = "drop"
#   )

cohort_partial_raw <- cohort_raw %>%
  filter(age >= 20, age <= 60) %>%
  group_by(cohort) %>%
  filter(n_distinct(age) >= 25) %>%
  ungroup()

cohort_complete_profile <- cohort_partial_raw %>%
  group_by(cohort) %>%
  summarise(
    pct_chinese = weighted.mean(pct_chinese, w = `All ethnicities`, na.rm = TRUE),
    population = sum(`All ethnicities`, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# Connect same age between adjacent cohorts
link_df <- cohort_partial_raw %>%
  dplyr::select(cohort, age, pct_chinese) %>%
  dplyr::rename(y = pct_chinese) %>%
  dplyr::inner_join(
    cohort_partial_raw %>%
      dplyr::select(cohort, age, pct_chinese) %>%
      dplyr::transmute(
        cohort = cohort - 1,   # shift so cohort t links to cohort t+1
        age = age,
        yend = pct_chinese
      ),
    by = c("cohort", "age")
  ) %>%
  mutate(
    x = cohort + line_pad,
    xend = cohort + 1 - line_pad
  )

# Plot
line_pad <- 0.04

ggplot() +
  geom_segment(
    data = link_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend
    ),
    color = "grey80",
    alpha = 0.25,
    linewidth = 0.3,
    lineend = "round"
  ) +
  geom_point(
    data = cohort_partial_raw,
    aes(x = cohort, y = pct_chinese),
    color = "#E08DBA",
    alpha = 0.25,
    size = 1
  ) +
  geom_point(
    data = cohort_complete_profile,
    aes(x = cohort, y = pct_chinese),
    color = "#6C1D47",
    size = 1.4
  ) +
  scale_y_continuous(
    limits = c(60, 90),
    breaks = seq(60, 90, by = 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    breaks = seq(
      min(cohort_complete_profile$cohort),
      max(cohort_complete_profile$cohort),
      by = 5
    ),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )



# Testing for statistically significant decline in % Chinese across birth cohort ----
# Chinese = numerator, All ethnicites = denominator
# "Non-Chinese" = 'failures' in the binomial model

# Currently filtering for cohorts that contain ages between 20-60
# Perhaps better to filter by >= 20 in cohort_profile$n_cells? i.e., cohorts which include
# greater than or equal to 20 ages

# Cohort Trend Model -----------------------------------------------------------
# Chinese = numerator, All ethnicites = denominator
# "Non-Chinese" = 'failures' in the binomial model

# Cohort measured in 10-year units
trend_df <- pct_chinese_complete %>%
  filter(
    age >= 20,
    age <= 60,
    !is.na(Chinese),
    !is.na(`All ethnicities`)
  ) %>%
  mutate(
    cohort = year - age,
    non_chinese = `All ethnicities` - Chinese,
    # Centre cohort by subtracting the average cohort year from every cohort value
    cohort10 = (cohort - mean(cohort, na.rm = TRUE)) / 10
  )

# Simplest model: tests whether later cohorts have lower % Chinese (without adjusting for age)
m1 <- glm(
  cbind(Chinese, non_chinese) ~ cohort10,
  family = quasibinomial(),
  data = trend_df
)

summary(m1)

# Odds ratio: 0.907 (<1 = 10-year-later cohort has lower odds of being classified as Chinese)
exp(coef(m1)["cohort10"])


# Use a quasibinomial model: When a response variable is binary but has a higher variance 
# than would be predicted by a binomial distribution, the quasibinomial model is utilized. 
# This could happen if the response variable has excessive dispersion or additional variation 
# that the model is not taking into account. 

# Model adjusting for age (i.e., comparing counts at the same age, do later birth cohorts have
# lower %Chinese?)

# glm() estimates coefficients
m2 <- glm(
  cbind(Chinese, non_chinese) ~ cohort10 + factor(age),
  family = quasibinomial(),
  data = trend_df
)

# cohort10      -0.107223   0.002212 -48.479   <2e-16 *** (Estimate, Std. Error, t-value, Pr(>|t|))
# The quasibinomial dispersion is 126.15 (very large) = there is substantial 
# overdispersion relative to an ordinary binomial mode (justifying quasibinomimal model)
summary(m2)

# vcovCL adjusts SEs, CIs, and p-values to account for within-cohort dependence
# cluster = ~ cohort treats rows with the same 'cohort' as belonging to the same cluster
coeftest(m2, vcov = vcovCL(m2, cluster = ~ cohort))

# Odds ratio: 0.898
exp(coef(m2)["cohort10"])


# Fixed-age cohort trend (i.e., model fitted at one age only)
age30_df <- trend_df %>%
  filter(age == 30)

m_age30 <- glm(
  cbind(Chinese, non_chinese) ~ cohort10,
  family = quasibinomial(),
  data = age30_df
)

# cohort10    -0.085045   0.009863  -8.623 5.24e-11 *** (Estimate, Std. Error, t-value, Pr(>|t|))
summary(m_age30)

# Odds ratio: 0.918
exp(coef(m_age30)["cohort10"])

# Age-standardized fitted model plot -------------------------------------------
# NOT WORKING
# Set the same age range as in m2
ages_std <- 20:60

# Standardized weights based on the number of people observed at each age
std_weights <- trend_df %>%
  filter(age %in% ages_std) %>%
  group_by(age) %>%
  summarise(w = sum(`All ethnicities`, na.rm = TRUE), .groups = "drop") %>%
  mutate(w = w / sum(w))

# Keep cohorts that have all ages between 20-60
complete_cohorts <- trend_df %>%
  filter(age %in% ages_std) %>%
  group_by(cohort) %>%
  summarise(n_ages = n_distinct(age), .groups = "drop") %>%
  filter(n_ages == length(ages_std)) %>%
  pull(cohort)

# Calculate age-standardized cohort profile
obs_cohort_std <- trend_df %>%
  filter(cohort %in% complete_cohorts, age %in% ages_std) %>%
  mutate(p = Chinese / `All ethnicities`) %>%
  left_join(std_weights, by = "age") %>%
  group_by(cohort) %>%
  summarise(
    pct_chinese_std = 100 * sum(w * p, na.rm = TRUE),
    .groups = "drop"
  )

# Model-based age-standardized fitted values from m2
pred_grid <- expand_grid(
  cohort = complete_cohorts,
  age = ages_std
) %>%
  left_join(std_weights, by = "age") %>%
  mutate(
    cohort10 = (cohort - mean(trend_df$cohort, na.rm = TRUE)) / 10
  )

# Fitted values from m2
X <- model.matrix(delete.response(terms(m2)), data = pred_grid)

beta_hat <- coef(m2)
V <- vcov(m2)

eta_hat <- as.vector(X %*% beta_hat)

pred_grid <- pred_grid %>%
  mutate(
    fit = plogis(eta_hat)
  )

fit_cohort_std <- pred_grid %>%
  group_by(cohort) %>%
  summarise(
    fit_pct = 100 * sum(w * fit),
    .groups = "drop"
  )

# Calculate 95% CI by simulating from coefficient covariance matrix
set.seed(123)

B <- 2000
beta_sim <- MASS::mvrnorm(B, mu = beta_hat, Sigma = V)

eta_sim <- X %*% t(beta_sim)       # rows = age/cohort cells, cols = simulations
p_sim <- plogis(eta_sim)

cohort_index <- split(seq_len(nrow(pred_grid)), pred_grid$cohort)

ci_df <- lapply(names(cohort_index), function(cc) {
  idx <- cohort_index[[cc]]
  w_i <- pred_grid$w[idx]
  
  std_sim <- colSums(p_sim[idx, , drop = FALSE] * w_i)
  
  tibble(
    cohort = as.numeric(cc),
    lower_pct = 100 * quantile(std_sim, 0.025),
    upper_pct = 100 * quantile(std_sim, 0.975)
  )
}) %>%
  bind_rows()

# Combine and plot
plot_df <- fit_cohort_std %>%
  left_join(ci_df, by = "cohort") %>%
  left_join(obs_cohort_std, by = "cohort")

ggplot() +
  geom_ribbon(
    data = plot_df,
    aes(x = cohort, ymin = lower_pct, ymax = upper_pct),
    fill = "grey75",
    alpha = 0.35
  ) +
  geom_line(
    data = plot_df,
    aes(x = cohort, y = fit_pct),
    linewidth = 1,
    color = "black"
  ) +
  geom_line(
    data = plot_df,
    aes(x = cohort, y = pct_chinese_std),
    linewidth = 0.8,
    color = "grey35"
  ) +
  geom_point(
    data = plot_df,
    aes(x = cohort, y = pct_chinese_std),
    size = 1.2,
    color = "grey35"
  ) +
  labs(
    title = "Age-standardized cohort profile of residents classified as Chinese",
    subtitle = "Observed standardized profile with age-adjusted fitted trend and 95% CI",
    x = "Birth cohort",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

#-------------------------------------------------------------------------------
# Pct-Chinese x Age (across all cohorts)
cohort_age_df <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort
  ) %>%
  arrange(cohort, age)

# Check the number of ages per cohort
cohort_availability <- cohort_age_df %>%
  summarise(
    n_years_observed = sum(!is.na(age)),
    .by = cohort
  ) %>%
  arrange(cohort)

ggplot(cohort_age_df, aes(x = age, y = pct_chinese, group = cohort, colour = cohort)) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  scale_color_viridis_c(guide = "none") +
  labs(
    title = "% Chinese by age across birth cohorts",
    x = "Age",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

ggplot(cohort_age_df, aes(x = year, y = pct_chinese, group = cohort, colour = cohort)) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  scale_color_viridis_c(guide = "none") +
  labs(
    title = "% Chinese by age across birth cohorts",
    x = "Age",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

# Line chart with small multiple -----------------------------------------------
# https://r-graph-gallery.com/web-line-chart-small-multiple-all-group-greyed-out.html
# % Chinese by Age (across all cohorts)

pct_chinese_age30 <- cohort_age_df %>%
  filter(age == 30) %>%
  distinct(cohort, .keep_all = TRUE) %>%
  select(cohort, year, age, pct_chinese) %>%
  arrange(cohort)

pct_chinese_age40 <- cohort_age_df %>%
  filter(age == 40) %>%
  distinct(cohort, .keep_all = TRUE) %>%
  select(cohort, year, age, pct_chinese) %>%
  arrange(cohort)

# Calculate the population-weighted mean pct_chinese
# Weighted by the age-structure (i.e., population per age group)
overall_mean_wt <- weighted.mean(
  pct_chinese_complete$pct_chinese,
  w = pct_chinese_complete$`All ethnicities`,
  na.rm = TRUE
)

# 1950 Cohort
cohort_age_df <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    highlight = if_else(cohort == 1950, "1950 cohort", "Other cohorts")
  ) %>%
  arrange(cohort, age)

# Add a point at age 30
point_30 <- cohort_age_df %>%
  filter(highlight == "1950 cohort", age == 30) %>%
  select(year, age, cohort, pct_chinese)

ggplot(cohort_age_df, aes(x = age, y = pct_chinese, group = cohort)) +
  geom_line(
    data = subset(cohort_age_df, highlight == "Other cohorts"),
    color = "grey85",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_line(
    data = subset(cohort_age_df, highlight == "1950 cohort"),
    color = "black",
    linewidth = 1
  ) +
  # Draw point at age 30
  geom_point(
    data = point_30,
    color = "black",
    size = 2
  ) +
  # Draw population-weighted mean (at 75.45%)
  geom_hline(
    yintercept = overall_mean_wt,
    color = "black",
    linewidth = 0.25
  ) +
  theme_minimal(base_size = 13)


# 1965 Cohort
cohort_age_df <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    highlight = if_else(cohort == 1965, "1965 cohort", "Other cohorts")
  ) %>%
  arrange(cohort, age)

# Add a point at age 30
point_30 <- cohort_age_df %>%
  filter(highlight == "1965 cohort", age == 30) %>%
  select(year, age, cohort, pct_chinese)

ggplot(cohort_age_df, aes(x = age, y = pct_chinese, group = cohort)) +
  geom_line(
    data = subset(cohort_age_df, highlight == "Other cohorts"),
    color = "grey85",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_line(
    data = subset(cohort_age_df, highlight == "1965 cohort"),
    color = "black",
    linewidth = 1
  ) +
  # Draw point at age 30
  geom_point(
    data = point_30,
    color = "black",
    size = 2
  ) +
  # Draw population-weighted mean (at 75.45%)
  geom_hline(
    yintercept = overall_mean_wt,
    color = "black",
    linewidth = 0.25
  ) +
  theme_minimal(base_size = 13)


# 1980 Cohort
cohort_age_df <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    highlight = if_else(cohort == 1980, "1980 cohort", "Other cohorts")
  ) %>%
  arrange(cohort, age)

# Add a point at age 30
point_30 <- cohort_age_df %>%
  filter(highlight == "1980 cohort", age == 30) %>%
  select(year, age, cohort, pct_chinese)

ggplot(cohort_age_df, aes(x = age, y = pct_chinese, group = cohort)) +
  geom_line(
    data = subset(cohort_age_df, highlight == "Other cohorts"),
    color = "grey85",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_line(
    data = subset(cohort_age_df, highlight == "1980 cohort"),
    color = "black",
    linewidth = 1
  ) +
  # Draw point at age 30
  geom_point(
    data = point_30,
    color = "black",
    size = 2
  ) +
  # Draw population-weighted mean (at 75.45%)
  geom_hline(
    yintercept = overall_mean_wt,
    color = "black",
    linewidth = 0.25
  ) +
  theme_minimal(base_size = 13)









#----------------

cohort_age_5yr <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    cohort_5yr = floor(cohort / 5) * 5
  ) %>%
  group_by(cohort_5yr, age) %>%
  summarise(
    pct_chinese = mean(pct_chinese, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(cohort_age_5yr, aes(x = age, y = pct_chinese, group = cohort_5yr, colour = cohort_5yr)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_color_viridis_c(guide = "none") +
  labs(
    title = "% Chinese by age across 5-year birth cohorts",
    x = "Age",
    y = "% Chinese",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 13)



# 1980 Cohort
cohort_age_df <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    highlight = if_else(cohort == 1980, "1980 cohort", "Other cohorts")
  ) %>%
  arrange(cohort, age)

# Add a point at age 30
point_30 <- cohort_age_df %>%
  filter(highlight == "1980 cohort", age == 30) %>%
  select(year, age, cohort, pct_chinese)

ggplot(cohort_age_df, aes(x = age, y = pct_chinese, group = cohort)) +
  geom_line(
    data = subset(cohort_age_df, highlight == "Other cohorts"),
    color = "grey85",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_line(
    data = subset(cohort_age_df, highlight == "1980 cohort"),
    color = "black",
    linewidth = 1
  ) +
  # Draw point at age 30
  geom_point(
    data = point_30,
    color = "black",
    size = 2
  ) +
  # Draw population-weighted mean (at 75.45%)
  geom_hline(
    yintercept = overall_mean_wt,
    color = "black",
    linewidth = 0.25
  ) +
  theme_minimal(base_size = 13)




#----------------

cohort_age_5yr <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort,
    cohort_5yr = floor(cohort / 5) * 5
  ) %>%
  group_by(cohort_5yr, age) %>%
  summarise(
    pct_chinese = mean(pct_chinese, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(cohort_age_5yr, aes(x = age, y = pct_chinese, group = cohort_5yr, colour = cohort_5yr)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_color_viridis_c(guide = "none") +
  labs(
    title = "% Chinese by age across 5-year birth cohorts",
    x = "Age",
    y = "% Chinese",
    colour = "Cohort"
  ) +
  theme_minimal(base_size = 13)

# selected_cohorts <- c(1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000)
# 
# cohort_age_selected <- pct_chinese_complete %>%
#   mutate(
#     cohort = if (!"cohort" %in% names(.)) year - age else cohort
#   ) %>%
#   filter(cohort %in% selected_cohorts)
# 
# ggplot(cohort_age_selected, aes(x = age, y = pct_chinese, group = cohort, colour = factor(cohort))) +
#   geom_line(linewidth = 0.9) +
#   labs(
#     title = "% Chinese by age for selected birth cohorts",
#     x = "Age",
#     y = "% Chinese",
#     colour = "Birth cohort"
#   ) +
#   theme_minimal(base_size = 13)

#-------------------------------------------------------------------------------
# Version with only cohorts that have observations for every age in the retained age range
required_ages <- 20:50

cohort_age_df_complete <- pct_chinese_complete %>%
  mutate(
    cohort = if (!"cohort" %in% names(.)) year - age else cohort
  ) %>%
  filter(
    !is.na(pct_chinese),
    age %in% required_ages
  ) %>%
  group_by(cohort) %>%
  filter(all(required_ages %in% age)) %>%
  ungroup() %>%
  arrange(cohort, age)

ggplot(cohort_age_df_complete, aes(x = age, y = pct_chinese, group = cohort, colour = cohort)) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  scale_color_viridis_c(guide = "none") +
  labs(
    title = "% Chinese by age across birth cohorts (ages 20–50 only)",
    x = "Age",
    y = "% Chinese"
  ) +
  theme_minimal(base_size = 13)

cohort_age_df_complete %>%
  distinct(cohort) %>%
  arrange(cohort)
  
  
  
  theme_minimal(base_size = 13)