#=============================================================================
# Rapid shift of a ground-nesting shorebird from natural grassland to
# agricultural habitats in Central Europe
#
# Collared Pratincole, Hortobagy National Park, 1950-2005.

#   MODULE 1  Population dynamics
#             1a  number of colonies       ~ year + year^2     (negative binomial)
#             1b  number of breeding pairs ~ year + year^2     (negative binomial)
#
#   MODULE 2  Livestock association, controlling for year
#             2a  number of colonies       ~ livestock + year  (negative binomial)
#             2b  number of breeding pairs ~ livestock + year  (negative binomial)
#
#   MODULE 3  Grazing context (natural habitats only)
#             3a  colonies at high-grazing sites, per year     (binomial)
#             3b  breeding pairs, high vs low grazing x year   (negative binomial)
#
#   MODULE 4  Habitat use (all habitats)
#             4a  colonies in artificial habitat, per year     (beta-binomial)
#             4b  breeding pairs, natural vs artificial x year (negative binomial)
#
# UNIT OF ANALYSIS
#   A colony is identified by its Colonycode, not by a data row. Several rows
#   can describe alternative sites of the same colony in the same year and
#   repeat that colony's size on each row, so counting rows would double-count
#   both colonies and breeding pairs.
#
# REPORTING CONVENTION
#   Every estimate is a posterior MEDIAN with a 95% credible interval (CrI),
#   plus pd (probability of direction: how consistently the effect sits on one
#   side of zero). We never use brms' default posterior means.
#=============================================================================


#=============================================================================
# SETUP
#=============================================================================

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(brms)
library(posterior)
library(ggplot2)
library(scales)
library(patchwork)
library(knitr)
library(kableExtra)

# Other packages also define select/filter/mutate; make sure we get dplyr's.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate


# ---------------------------------------------------------------------------
# summarise a fitted model's coefficients
# ---------------------------------------------------------------------------
report_fixef <- function(model, term_labels = NULL) {
  draws <- as_draws_df(model)
  coef_names <- grep("^b_", variables(model), value = TRUE)   # fixed effects only
  
  summaries <- lapply(coef_names, function(coef_name) {
    x <- draws[[coef_name]]
    data.frame(
      Term     = sub("^b_", "", coef_name),
      Estimate = median(x),                                   # posterior median
      CI_low   = unname(quantile(x, 0.025)),
      CI_high  = unname(quantile(x, 0.975)),
      P_gt0    = mean(x > 0),                                 # P(effect > 0)
      pd       = max(mean(x > 0), mean(x < 0)),               # probability of direction
      row.names = NULL
    )
  })
  
  result <- do.call(rbind, summaries)
  if (!is.null(term_labels)) result$Term <- term_labels
  result
}


# ---------------------------------------------------------------------------
# turn a summary into formatted table text
# ---------------------------------------------------------------------------
format_for_table <- function(summary_df) {
  data.frame(
    Term        = summary_df$Term,
    Estimate    = sprintf("%.2f", summary_df$Estimate),
    `95% CrI`   = sprintf("[%.2f, %.2f]", summary_df$CI_low, summary_df$CI_high),
    `P(>0)`     = sprintf("%.3f", summary_df$P_gt0),
    pd          = sprintf("%.3f", summary_df$pd),
    check.names = FALSE
  )
}


# ---------------------------------------------------------------------------
# Shared plot theme
# ---------------------------------------------------------------------------
theme_panel <- theme_classic(base_size = 13) +
  theme(
    plot.title  = element_text(size = 14, face = "plain", hjust = 0, margin = margin(b = 6)),
    axis.line   = element_line(colour = "black", linewidth = 0.6),
    axis.ticks  = element_line(colour = "black", linewidth = 0.6),
    axis.title  = element_text(size = 12),
    axis.text   = element_text(size = 11, colour = "black"),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )


# ---------------------------------------------------------------------------
# Priors and MCMC settings (shared by every model)
# ---------------------------------------------------------------------------
# Weakly informative priors: wide enough not to drive the result, tight enough
# to keep sampling stable.
priors_negbinom <- c(prior(normal(0, 5),        class = "Intercept"),
                     prior(normal(0, 3),        class = "b"),
                     prior(gamma(0.01, 0.01),   class = "shape"))

priors_binomial <- c(prior(normal(0, 5),        class = "Intercept"),
                     prior(normal(0, 3),        class = "b"))

# 4 chains x 4000 iterations, first 2000 discarded as warm-up.
MCMC <- list(chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 123,
             control = list(adapt_delta = 0.95))

# Thin wrapper so every model uses the same settings.
fit_model <- function(...) {
  brm(..., chains = MCMC$chains, iter = MCMC$iter, warmup = MCMC$warmup,
      cores = MCMC$cores, seed = MCMC$seed, control = MCMC$control, refresh = 0)
}


#=============================================================================
# DATA PREPARATION
#=============================================================================

# ---------------------------------------------------------------------------
# STEP 1. Read the colony records
# ---------------------------------------------------------------------------
raw_data <- read_excel("Saduakas et al._Data1.xlsx",sheet = "Sheet1")
names(raw_data) <- trimws(names(raw_data))          # some headers have stray spaces
raw_data$Annual_colony_size_mean <- as.numeric(raw_data$Annual_colony_size_mean)
raw_data$grazing_intensity       <- as.numeric(raw_data$grazing_intensity)


# ---------------------------------------------------------------------------
# STEP 2. One shared year scaling, used everywhere
# ---------------------------------------------------------------------------
# Standardising year (mean 0, SD 1) helps sampling and makes the year
# coefficient comparable across all modules. Defining it once guarantees every
# module uses the same scale.
year_mean <- mean(raw_data$Nesting_activity_year)
year_sd   <- sd(raw_data$Nesting_activity_year)
scale_year <- function(year) (year - year_mean) / year_sd


# ---------------------------------------------------------------------------
# STEP 3. Collapse to one row per colony per year
# ---------------------------------------------------------------------------
# This is the key step. Rows sharing a Colonycode within a year are alternative
# sites of one colony, each carrying that colony's size. Taking `unique()`
# before summing counts the colony's size once instead of repeating it.
# Habitat and grazing are constant within a colony-year, so `first()` is safe.
colony_year <- raw_data %>%
  mutate(colony_id = factor(paste0("C", Colonycode))) %>%
  group_by(colony_id, Nesting_activity_year) %>%
  summarise(
    colony_size       = sum(unique(Annual_colony_size_mean[!is.na(Annual_colony_size_mean)])),
    Habitat_type      = dplyr::first(Habitat_type),
    grazing_intensity = dplyr::first(grazing_intensity),
    Sitename          = dplyr::first(Sitename),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------
# STEP 4. Yearly totals (used by Modules 1 and 2)
# ---------------------------------------------------------------------------
annual <- colony_year %>%
  group_by(Nesting_activity_year) %>%
  summarise(
    n_colonies  = n_distinct(colony_id),
    total_pairs = round(sum(colony_size, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Year = Nesting_activity_year) %>%
  mutate(
    year_scaled  = scale_year(Year),
    year_scaled2 = year_scaled^2          # quadratic term: allows rise-then-fall
  )


# ---------------------------------------------------------------------------
# STEP 5. Colony-level subsets (used by Modules 3 and 4)
# ---------------------------------------------------------------------------
# Grazing intensity is recorded 0-3; we collapse it to Low (0-1) vs High (2-3).
# Artificial-habitat colonies were all ungrazed, so including them would make
# habitat and grazing indistinguishable. Module 3 therefore uses natural only.
natural_colony_years <- colony_year %>%
  filter(Habitat_type == "Natural") %>%
  mutate(
    grazing_cat = case_when(
      grazing_intensity %in% c(0, 1) ~ "Low",
      grazing_intensity %in% c(2, 3) ~ "High",
      TRUE                           ~ NA_character_
    ),
    grazing_cat = factor(grazing_cat, levels = c("Low", "High")),
    year_scaled = scale_year(Nesting_activity_year)
  ) %>%
  filter(!is.na(grazing_cat)) %>%
  mutate(colony_id = droplevels(colony_id))

# Module 4 uses every colony, natural and artificial.
all_colony_years <- colony_year %>%
  filter(!is.na(Habitat_type)) %>%
  mutate(
    Habitat_type = factor(Habitat_type, levels = c("Natural", "Artificial")),
    colony_id    = droplevels(colony_id),
    year_scaled  = scale_year(Nesting_activity_year)
  )


#=============================================================================
# MODULE 1 - POPULATION DYNAMICS
#=============================================================================
#
# QUESTION
#   How did colony number and breeding-pair number change from 1950 to 2005?
#
#   Both counts vary much more than a Poisson allows (variance/mean ~2.9 for
#   colonies, far higher for pairs). Negative binomial adds a parameter that
#   absorbs this extra variability.
#=============================================================================

model_colony_trend <- fit_model(
  n_colonies ~ year_scaled + year_scaled2,
  family = negbinomial(), data = annual, prior = priors_negbinom
)

model_pair_trend <- fit_model(
  total_pairs ~ year_scaled + year_scaled2,
  family = negbinomial(), data = annual, prior = priors_negbinom
)

trend_terms <- c("Intercept", "Year (scaled)", "Year (scaled)\u00b2")
res_colony_trend <- report_fixef(model_colony_trend, trend_terms)
res_pair_trend   <- report_fixef(model_pair_trend,   trend_terms)

cat("\n== Module 1: population dynamics ==\n")
print(res_colony_trend)
print(res_pair_trend)


# ---------------------------------------------------------------------------
# Figure 1: fitted trajectories with the observed yearly totals
# ---------------------------------------------------------------------------
# Build a smooth grid of years, predict on it, and convert back to real years
# for the x-axis.
year_grid <- data.frame(
  year_scaled = seq(min(annual$year_scaled), max(annual$year_scaled), length.out = 300)
)
year_grid$year_scaled2 <- year_grid$year_scaled^2
year_grid$Year <- year_mean + year_sd * year_grid$year_scaled

# predict on the grid and return a tidy curve with its CrI ribbon.
predict_trend <- function(model) {
  prediction <- fitted(model, newdata = year_grid)
  cbind(year_grid,
        predicted = prediction[, "Estimate"],
        conf.low  = prediction[, "Q2.5"],
        conf.high = prediction[, "Q97.5"])
}

curve_colony_trend <- predict_trend(model_colony_trend)
curve_pair_trend   <- predict_trend(model_pair_trend)

fig_1a <- ggplot() +
  geom_ribbon(data = curve_colony_trend, aes(Year, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_colony_trend, aes(Year, predicted), linewidth = 1.15) +
  geom_point(data = annual, aes(Year, n_colonies), size = 2, alpha = 0.8) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(5), expand = expansion(mult = c(0, 0.06))) +
  labs(title = "(a)", x = "Year", y = "Number of colonies") +
  theme_panel

fig_1b <- ggplot() +
  geom_ribbon(data = curve_pair_trend, aes(Year, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_pair_trend, aes(Year, predicted), linewidth = 1.15) +
  geom_point(data = annual, aes(Year, total_pairs), size = 2, alpha = 0.8) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(5), expand = expansion(mult = c(0, 0.06))) +
  labs(title = "(b)", x = "Year", y = "Number of breeding pairs") +
  theme_panel

fig_population <- fig_1a | fig_1b
fig_population

# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
n_trend_rows <- nrow(res_colony_trend)

tab1 <- bind_rows(format_for_table(res_colony_trend),
                  format_for_table(res_pair_trend)) %>%
  kbl(caption = "Temporal change in the breeding population of Collared Pratincoles, 1950-2005",
      align = "lcccc") %>%
  pack_rows("Model 1: Number of colonies (negative binomial)",
            1, n_trend_rows) %>%
  pack_rows("Model 2: Number of breeding pairs (negative binomial)",
            n_trend_rows + 1, n_trend_rows + nrow(res_pair_trend)) %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  footnote(general = paste(
    "Posterior medians with 95% credible intervals (CrI).",
    "Year is standardised; a negative quadratic term indicates a rise-then-decline trajectory.",
    "P(>0): posterior probability the effect is positive. pd: probability of direction.",
    "All Rhat = 1.00; effective sample sizes > 2000."))


#=============================================================================
# MODULE 2 - LIVESTOCK ASSOCIATION (controlling for year)
#=============================================================================
#
# QUESTION
#   Both the pratincole population and regional livestock declined over time.
#   Does livestock still explain bird numbers once we account for that shared
#   downward trend? Every model below therefore includes year as a control,
#   and we read off the livestock coefficient.
#
# WHAT WE FIT
#   For each livestock type (sheep, cattle, horse, pig, and the combined
#   livestock-unit index):
#       number of colonies       ~ livestock + year
#       number of breeding pairs ~ livestock + year
#   That is 5 livestock types x 2 responses = 10 models.
#
# HOW TO READ THE RESULT
#   Both predictors are standardised, so a coefficient is "change in log-count
#   per one SD more livestock, at a fixed year". Positive = more livestock goes
#   with more birds.
#=============================================================================

# ---------------------------------------------------------------------------
# STEP 1. Load the regional livestock counts
# ---------------------------------------------------------------------------
# One row per year. Some early years have no pig record; treat those as zero.
livestock <- read_csv("Saduakas et al._Data2.csv") %>%
  filter(Year <= 2005) %>%
  mutate(Pig = ifelse(is.na(Pig), 0, Pig))

livestock_types <- c("Sheep", "Cattle", "Horse", "Pig", "Livestock_Units")
livestock_label <- c(Sheep = "Sheep", Cattle = "Cattle", Horse = "Horse",
                     Pig = "Pig", Livestock_Units = "Total (livestock units)")


# ---------------------------------------------------------------------------
# STEP 2. Helper: pull a tidy summary of the livestock effect from one model
# ---------------------------------------------------------------------------
summarise_livestock_effect <- function(model, response_name, livestock_name) {
  effect <- as_draws_df(model)[["b_livestock_z"]]   # draws for the livestock slope
  
  data.frame(
    Response  = response_name,
    Livestock = livestock_name,
    Estimate  = median(effect),
    CI_low    = unname(quantile(effect, 0.025)),
    CI_high   = unname(quantile(effect, 0.975)),
    pd        = max(mean(effect > 0), mean(effect < 0))
  )
}


# ---------------------------------------------------------------------------
# STEP 3. Fit the 10 models, one livestock type at a time
# ---------------------------------------------------------------------------
colony_models <- list()   # colonies ~ livestock + year
pair_models   <- list()   # breeding pairs ~ livestock + year
effect_rows   <- list()   # tidy summaries for the table

for (type in livestock_types) {
  
  # 3a. Join yearly bird totals to this livestock type, keeping only years
  #     present in both datasets (1977-2005).
  model_data <- livestock %>%
    select(Year, livestock_count = all_of(type)) %>%
    inner_join(select(annual, Year, n_colonies, total_pairs), by = "Year") %>%
    arrange(Year)
  
  # 3b. Standardise both predictors so coefficients are comparable.
  model_data$livestock_z <- as.numeric(scale(model_data$livestock_count))
  model_data$year_z      <- scale_year(model_data$Year)
  
  # 3c. Fit both responses.
  m_colonies <- fit_model(n_colonies  ~ livestock_z + year_z,
                          family = negbinomial(), data = model_data, prior = priors_negbinom)
  m_pairs    <- fit_model(total_pairs ~ livestock_z + year_z,
                          family = negbinomial(), data = model_data, prior = priors_negbinom)
  
  # 3d. Store the models and their summaries.
  colony_models[[type]] <- m_colonies
  pair_models[[type]]   <- m_pairs
  effect_rows[[type]]   <- bind_rows(
    summarise_livestock_effect(m_colonies, "Number of colonies",       livestock_label[[type]]),
    summarise_livestock_effect(m_pairs,    "Number of breeding pairs", livestock_label[[type]])
  )
}

livestock_effects <- bind_rows(effect_rows)
cat("\n== Module 2: livestock association ==\n")
print(livestock_effects)


# ---------------------------------------------------------------------------
# STEP 4. Table 2 - colonies and pairs side by side
# ---------------------------------------------------------------------------
tab2 <- livestock_effects %>%
  mutate(cell = sprintf("%.2f [%.2f, %.2f] (%.2f)", Estimate, CI_low, CI_high, pd)) %>%
  select(Livestock, Response, cell) %>%
  pivot_wider(names_from = Response, values_from = cell) %>%
  kbl(caption = "Year-controlled associations between regional livestock numbers and breeding abundance, 1977-2005",
      col.names = c("Livestock", "Number of colonies", "Number of breeding pairs"),
      align = "lcc") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  footnote(general = paste(
    "Each cell: posterior median beta [95% CrI] (probability of direction, pd).",
    "Coefficients are standardised (per SD of the livestock predictor), on the log scale,",
    "from models including year to control for the shared temporal trend.",
    "Colonies and breeding pairs: negative binomial. Overlap 1977-2005 (21 years).",
    "Livestock numbers are regional totals (Borza et al. 2017), a coarse index of local grazing.",
    "All Rhat = 1.00; effective sample sizes > 2000."))


# ---------------------------------------------------------------------------
# STEP 5. Figure 2, panel (a): livestock time series
# ---------------------------------------------------------------------------
livestock_palette <- c("Sheep" = "grey40", "Cattle" = "#1f77b4", "Horse" = "#d62728",
                       "Pig" = "#2ca02c", "Total (livestock units)" = "black")
livestock_order   <- c("Sheep", "Cattle", "Horse", "Pig", "Total (livestock units)")

# Sheep and total units are large; cattle/horse/pig are much smaller. Put the
# big series on the left axis and rescale the small ones onto a right axis so
# all five fit in one plot.
right_axis_max <- max(c(livestock$Cattle, livestock$Horse, livestock$Pig))
left_axis_max  <- max(c(livestock$Sheep, livestock$Livestock_Units))
scale_factor   <- left_axis_max / right_axis_max

fig_2a <- ggplot(livestock, aes(x = Year)) +
  # Left-axis series, drawn at their true values:
  geom_line(aes(y = Sheep, colour = "Sheep"), linewidth = 1) +
  geom_point(aes(y = Sheep, colour = "Sheep"), size = 1.8) +
  geom_line(aes(y = Livestock_Units, colour = "Total (livestock units)"), linewidth = 2) +
  # Right-axis series, multiplied by scale_factor to share the plot:
  geom_line(aes(y = Cattle * scale_factor, colour = "Cattle"), linewidth = 0.8) +
  geom_point(aes(y = Cattle * scale_factor, colour = "Cattle"), size = 1.5, shape = 15) +
  geom_line(aes(y = Horse * scale_factor, colour = "Horse"), linewidth = 0.8) +
  geom_point(aes(y = Horse * scale_factor, colour = "Horse"), size = 1.5, shape = 17) +
  geom_line(aes(y = Pig * scale_factor, colour = "Pig"), linewidth = 0.8) +
  geom_point(aes(y = Pig * scale_factor, colour = "Pig"), size = 1.5, shape = 18) +
  scale_x_continuous(breaks = pretty_breaks(7)) +
  scale_y_continuous(
    name = "Number of sheep or livestock units", labels = comma, limits = c(0, NA),
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of animals", labels = comma)) +
  scale_colour_manual(name = NULL, values = livestock_palette, breaks = livestock_order) +
  labs(title = "(a)", x = "Year") +
  theme_panel + theme(legend.position = "none")

# A throwaway plot used only to harvest a clean legend for the 2x2 layout.
legend_source <- ggplot(
  data.frame(Year = rep(range(livestock$Year), 5), val = 1,
             grp = factor(rep(livestock_order, each = 2), levels = livestock_order)),
  aes(Year, val, colour = grp)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(name = NULL, values = livestock_palette, breaks = livestock_order) +
  theme_panel +
  theme(legend.position = "right", legend.text = element_text(size = 10),
        legend.key.width = unit(1.1, "cm"))

extract_legend <- function(plot) {
  grobs <- ggplotGrob(plot)$grobs
  grobs[[which(sapply(grobs, function(g) g$name) == "guide-box")[1]]]
}
legend_cell <- wrap_elements(full = extract_legend(legend_source))


# ---------------------------------------------------------------------------
# STEP 6. Figure 2, panels (b) and (c): the livestock effect alone
# ---------------------------------------------------------------------------
# These panels use the combined livestock-unit index. Rebuild exactly the data
# those two models were fitted on, so the x-axis spans the observed range.
lu_data <- livestock %>%
  select(Year, Livestock_Units) %>%
  inner_join(select(annual, Year, n_colonies, total_pairs), by = "Year") %>%
  arrange(Year)
lu_data$livestock_z <- as.numeric(scale(lu_data$Livestock_Units))
lu_data$year_z      <- scale_year(lu_data$Year)

# 6a. Fitted curve: predict across the livestock range with year fixed at its
#     mean (year_z = 0), isolating the livestock effect.
fitted_curve <- function(model) {
  lu_mean <- mean(lu_data$Livestock_Units)
  lu_sd   <- sd(lu_data$Livestock_Units)
  
  grid <- data.frame(Livestock_Units = seq(min(lu_data$Livestock_Units),
                                           max(lu_data$Livestock_Units),
                                           length.out = 200))
  grid$livestock_z <- (grid$Livestock_Units - lu_mean) / lu_sd
  grid$year_z      <- 0                       # hold year at its mean
  
  prediction <- fitted(model, newdata = grid)
  data.frame(Livestock_Units = grid$Livestock_Units,
             predicted = prediction[, "Estimate"],
             conf.low  = prediction[, "Q2.5"],
             conf.high = prediction[, "Q97.5"])
}

# 6b. Partial residuals: each observed year's count, adjusted to the mean year
#     so the dots are comparable to the year-held-constant curve.
#       adjusted = (prediction at mean year) x (observed / prediction at true year)
#     Note these adjusted values can exceed any real yearly count.
partial_residuals <- function(model, response_col) {
  pred_at_true_year <- fitted(model, newdata = lu_data)[, "Estimate"]
  
  data_at_mean_year <- lu_data
  data_at_mean_year$year_z <- 0
  pred_at_mean_year <- fitted(model, newdata = data_at_mean_year)[, "Estimate"]
  
  data.frame(
    Livestock_Units = lu_data$Livestock_Units,
    resid = pred_at_mean_year * (lu_data[[response_col]] / pred_at_true_year)
  )
}

curve_colonies  <- fitted_curve(colony_models[["Livestock_Units"]])
curve_pairs     <- fitted_curve(pair_models[["Livestock_Units"]])
points_colonies <- partial_residuals(colony_models[["Livestock_Units"]], "n_colonies")
points_pairs    <- partial_residuals(pair_models[["Livestock_Units"]],   "total_pairs")

fig_2b <- ggplot() +
  geom_ribbon(data = curve_colonies, aes(Livestock_Units, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_colonies, aes(Livestock_Units, predicted), linewidth = 1.2) +
  geom_point(data = points_colonies, aes(Livestock_Units, resid), size = 2, alpha = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(5), labels = comma, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(5), expand = expansion(mult = c(0, 0.06))) +
  labs(title = "(b)", x = NULL, y = "Colonies (adj.)") +
  theme_panel + theme(plot.margin = margin(6, 18, 6, 6))

fig_2c <- ggplot() +
  geom_ribbon(data = curve_pairs, aes(Livestock_Units, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_pairs, aes(Livestock_Units, predicted), linewidth = 1.2) +
  geom_point(data = points_pairs, aes(Livestock_Units, resid), size = 2, alpha = 0.7) +
  scale_x_continuous(breaks = pretty_breaks(5), labels = comma, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, NA), breaks = pretty_breaks(5), expand = expansion(mult = c(0, 0.06))) +
  labs(title = "(c)", x = NULL, y = "Breeding pairs (adj.)") +
  theme_panel + theme(plot.margin = margin(6, 14, 6, 18))


# ---------------------------------------------------------------------------
# STEP 7. Assemble Figure 2
# ---------------------------------------------------------------------------
shared_xlab <- wrap_elements(full = grid::textGrob(
  "Total livestock units (year held at mean)", gp = grid::gpar(fontsize = 12)))

bottom_row <- (fig_2b | plot_spacer() | fig_2c) + plot_layout(widths = c(1, 0.06, 1))
bottom_row <- bottom_row / shared_xlab + plot_layout(heights = c(1, 0.06))

fig_livestock <- (fig_2a | legend_cell) / bottom_row + plot_layout(heights = c(1, 1.1))

fig_livestock 

#=============================================================================
# MODULE 3 - GRAZING CONTEXT (natural habitats only)
#=============================================================================
#
# QUESTION
#   Did colonies increasingly sit at heavily grazed sites, and where were the
#   breeding pairs concentrated?
#
#   Grazing category barely changes within a colony (it varies across years in
#   only 1 of 34 natural colonies) or within a site (4 of 15). A colony or site
#   random intercept would therefore be completely separated: every group would
#   be all-High or all-Low, so the year effect would be driven by the prior
#   rather than the data. Aggregating to the year avoids this, and matches how
#   model 3b already treats the pairs.
#=============================================================================

# 3a. Of the N natural colonies present in a year, k sat at high-grazing sites.
grazing_annual <- natural_colony_years %>%
  group_by(Nesting_activity_year) %>%
  summarise(
    N = n_distinct(colony_id),
    k = n_distinct(colony_id[grazing_cat == "High"]),
    .groups = "drop"
  ) %>%
  mutate(year_scaled = scale_year(Nesting_activity_year))

model_grazing_colonies <- fit_model(
  k | trials(N) ~ year_scaled,
  family = binomial(), data = grazing_annual, prior = priors_binomial
)

# 3b. Total breeding pairs per year in each grazing category. `complete()` fills
#     absent category-years with 0 so the series has no gaps.
grazing_pairs_annual <- natural_colony_years %>%
  filter(!is.na(colony_size)) %>%
  group_by(Nesting_activity_year, grazing_cat) %>%
  summarise(total_pairs = round(sum(colony_size)), .groups = "drop") %>%
  complete(Nesting_activity_year, grazing_cat, fill = list(total_pairs = 0)) %>%
  mutate(year_scaled = scale_year(Nesting_activity_year))

# The interaction asks whether the two categories changed at different rates.
model_grazing_pairs <- fit_model(
  total_pairs ~ grazing_cat * year_scaled,
  family = negbinomial(), data = grazing_pairs_annual, prior = priors_negbinom
)

res_grazing_colonies <- report_fixef(model_grazing_colonies,
                                     c("Intercept", "Year (scaled)"))
res_grazing_pairs    <- report_fixef(model_grazing_pairs,
                                     c("Intercept (Low)", "Grazing: High",
                                       "Year (scaled)", "Grazing High x Year"))

cat("\n== Module 3: grazing ==\n")
print(res_grazing_colonies)
print(res_grazing_pairs)


#=============================================================================
# MODULE 4 - HABITAT USE (all habitats)
#=============================================================================
#
# QUESTION
#   Did colonies shift from natural grassland to artificial (agricultural)
#   habitat, and did the breeding pairs follow?
#
# WHY BETA-BINOMIAL FOR 4a
#   Habitat type never changes within a colony (0 of 35), so the same
#   separation problem as Module 3 applies and we again aggregate to the year.
#   The yearly proportions also scatter more than a plain binomial allows
#   (Pearson dispersion ~4), so beta-binomial adds the needed flexibility.
#=============================================================================

# 4a. Of the N colonies present in a year, k were in artificial habitat.
habitat_annual <- all_colony_years %>%
  group_by(Nesting_activity_year) %>%
  summarise(
    N = n_distinct(colony_id),
    k = n_distinct(colony_id[Habitat_type == "Artificial"]),
    .groups = "drop"
  ) %>%
  mutate(year_scaled = scale_year(Nesting_activity_year))

model_habitat_colonies <- fit_model(
  k | trials(N) ~ year_scaled,
  family = beta_binomial(), data = habitat_annual, prior = priors_binomial
)

# 4b. Total breeding pairs per year in each habitat type.
habitat_pairs_annual <- all_colony_years %>%
  filter(!is.na(colony_size)) %>%
  group_by(Nesting_activity_year, Habitat_type) %>%
  summarise(total_pairs = round(sum(colony_size)), .groups = "drop") %>%
  complete(Nesting_activity_year, Habitat_type, fill = list(total_pairs = 0)) %>%
  mutate(year_scaled = scale_year(Nesting_activity_year))

model_habitat_pairs <- fit_model(
  total_pairs ~ Habitat_type * year_scaled,
  family = negbinomial(), data = habitat_pairs_annual, prior = priors_negbinom
)

res_habitat_colonies <- report_fixef(model_habitat_colonies,
                                     c("Intercept", "Year (scaled)"))
res_habitat_pairs    <- report_fixef(model_habitat_pairs,
                                     c("Intercept (Natural)", "Habitat: Artificial",
                                       "Year (scaled)", "Artificial x Year"))

cat("\n== Module 4: habitat ==\n")
print(res_habitat_colonies)
print(res_habitat_pairs)



#=============================================================================
# FIGURE 3 - grazing (top row) and habitat (bottom row)
#=============================================================================

# Helper A: predicted proportion over time, for the binomial models (3a, 4a).
# Setting N = 1 makes `fitted()` return a proportion rather than a count.
predict_proportion <- function(model, model_data) {
  grid <- data.frame(
    year_scaled = seq(min(model_data$year_scaled), max(model_data$year_scaled),
                      length.out = 200)
  )
  grid$Year <- year_mean + year_sd * grid$year_scaled
  grid$N <- 1
  
  prediction <- fitted(model, newdata = grid)
  cbind(grid,
        predicted = prediction[, "Estimate"],
        conf.low  = prediction[, "Q2.5"],
        conf.high = prediction[, "Q97.5"])
}

# Helper B: predicted counts over time for each group, for the interaction
# models (3b, 4b). Produces one curve per group level.
predict_by_group <- function(model, model_data, group_column, group_levels) {
  grid <- expand.grid(
    year_scaled = seq(min(model_data$year_scaled), max(model_data$year_scaled),
                      length.out = 200),
    group = factor(group_levels, levels = group_levels)
  )
  names(grid)[2] <- group_column         # rename to the column the model expects
  grid$Year <- year_mean + year_sd * grid$year_scaled
  
  prediction <- fitted(model, newdata = grid)
  cbind(grid,
        predicted = prediction[, "Estimate"],
        conf.low  = prediction[, "Q2.5"],
        conf.high = prediction[, "Q97.5"])
}

curve_grazing_colonies <- predict_proportion(model_grazing_colonies, grazing_annual)
curve_grazing_pairs    <- predict_by_group(model_grazing_pairs, grazing_pairs_annual,
                                           "grazing_cat", c("Low", "High"))
curve_habitat_colonies <- predict_proportion(model_habitat_colonies, habitat_annual)
curve_habitat_pairs    <- predict_by_group(model_habitat_pairs, habitat_pairs_annual,
                                           "Habitat_type", c("Natural", "Artificial"))

grazing_palette <- c("Low" = "#2166ac", "High" = "#b2182b")
habitat_palette <- c("Natural" = "#1b7837", "Artificial" = "#762a83")

# Observed colony records as 0/1 outcomes, for the points in panels (a) and (c).
# The models are fitted on yearly counts, but a fitted probability is the same
# quantity either way, so the individual records can be plotted against it.
grazing_points <- natural_colony_years %>%
  mutate(occurred_at_high_grazing = ifelse(grazing_cat == "High", 1, 0))

habitat_points <- all_colony_years %>%
  mutate(occurred_in_artificial = ifelse(Habitat_type == "Artificial", 1, 0))

# Shared legend styling for the two right-hand panels.
small_legend <- theme(
  legend.position = "top", legend.justification = "right",
  legend.direction = "horizontal",
  legend.title = element_text(size = 8),
  legend.text  = element_text(size = 7),
  legend.key.size = unit(0.35, "cm"),
  legend.margin = margin(0, 0, 0, 0)
)

# Panel (a): predicted probability that a colony was at a high-grazing site.
# Points are the observed colony records (0 = low, 1 = high), jittered so that
# overlapping records stay visible.
fig_3a <- ggplot() +
  geom_ribbon(data = curve_grazing_colonies, aes(Year, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_grazing_colonies, aes(Year, predicted), linewidth = 1.15) +
  geom_jitter(data = grazing_points,
              aes(Nesting_activity_year, occurred_at_high_grazing),
              height = 0.05, width = 0.2, alpha = 0.30, size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Prob. of colony at high-grazing site",
                     limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1),
                     labels = c("0", "", "1"), expand = expansion(mult = c(0, 0))) +
  labs(title = "(a)", x = "Year") +
  theme_panel

# Panel (b): breeding pairs by grazing category. Log scale because the counts
# span three orders of magnitude; log1p keeps the zeros visible.
fig_3b <- ggplot() +
  geom_ribbon(data = curve_grazing_pairs,
              aes(Year, ymin = conf.low, ymax = conf.high, group = grazing_cat),
              fill = "grey70", alpha = 0.30) +
  geom_line(data = curve_grazing_pairs, aes(Year, predicted, colour = grazing_cat), linewidth = 1.15) +
  geom_point(data = grazing_pairs_annual,
             aes(Nesting_activity_year, total_pairs, colour = grazing_cat), size = 1.5, alpha = 0.5) +
  scale_colour_manual(name = "Grazing", values = grazing_palette) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Breeding pairs (log scale)", trans = "log1p",
                     breaks = c(0, 3, 10, 30, 100, 300, 1000),
                     expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "(b)", x = "Year") +
  theme_panel + small_legend

# Panel (c): predicted probability that a colony was in artificial habitat.
# Points are the observed colony records (0 = natural, 1 = artificial).
fig_3c <- ggplot() +
  geom_ribbon(data = curve_habitat_colonies, aes(Year, ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(data = curve_habitat_colonies, aes(Year, predicted), linewidth = 1.15) +
  geom_jitter(data = habitat_points,
              aes(Nesting_activity_year, occurred_in_artificial),
              height = 0.05, width = 0.2, alpha = 0.30, size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Prob. of colony in artificial habitat",
                     limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1),
                     labels = c("0", "", "1"), expand = expansion(mult = c(0, 0))) +
  labs(title = "(c)", x = "Year") +
  theme_panel

# Panel (d): breeding pairs by habitat type.
fig_3d <- ggplot() +
  geom_ribbon(data = curve_habitat_pairs,
              aes(Year, ymin = conf.low, ymax = conf.high, group = Habitat_type),
              fill = "grey70", alpha = 0.30) +
  geom_line(data = curve_habitat_pairs, aes(Year, predicted, colour = Habitat_type), linewidth = 1.15) +
  geom_point(data = habitat_pairs_annual,
             aes(Nesting_activity_year, total_pairs, colour = Habitat_type), size = 1.5, alpha = 0.5) +
  scale_colour_manual(name = "Habitat", values = habitat_palette) +
  scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Breeding pairs (log scale)", trans = "log1p",
                     breaks = c(0, 3, 10, 30, 100, 300, 1000),
                     expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "(d)", x = "Year") +
  theme_panel + small_legend

fig_grazing_habitat <- (fig_3a | fig_3b) / (fig_3c | fig_3d)
fig_grazing_habitat 
#=============================================================================
# TABLE 3 - all four grazing and habitat models
#=============================================================================

# Row counts, used to place the section headers in the stacked table.
n_rows_3a <- nrow(res_grazing_colonies)
n_rows_3b <- nrow(res_grazing_pairs)
n_rows_4a <- nrow(res_habitat_colonies)
n_rows_4b <- nrow(res_habitat_pairs)

tab3 <- bind_rows(
  format_for_table(res_grazing_colonies),
  format_for_table(res_grazing_pairs),
  format_for_table(res_habitat_colonies),
  format_for_table(res_habitat_pairs)
) %>%
  kbl(caption = "Colony occurrence and breeding-pair distribution with respect to grazing intensity and habitat type, 1950-2005",
      align = "lcccc") %>%
  pack_rows("Model 1: Colony occurrence in high-grazing sites",
            1,
            n_rows_3a) %>%
  pack_rows("Model 2: Breeding pairs in high- vs low-grazing sites",
            n_rows_3a + 1,
            n_rows_3a + n_rows_3b) %>%
  pack_rows("Model 3: Colony occurrence in artificial habitats",
            n_rows_3a + n_rows_3b + 1,
            n_rows_3a + n_rows_3b + n_rows_4a) %>%
  pack_rows("Model 4: Breeding pairs in natural vs artificial habitats",
            n_rows_3a + n_rows_3b + n_rows_4a + 1,
            n_rows_3a + n_rows_3b + n_rows_4a + n_rows_4b) %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  footnote(general = paste(
    "Posterior medians with 95% credible intervals (CrI).",
    "P(>0): posterior probability the effect is positive. pd: probability of direction.",
    "Grazing models (1-2): natural habitats only; colony model is binomial on annual colony counts.",
    "Habitat models (3-4): all habitats; colony model is beta-binomial on annual colony counts.",
    "Colony models on the logit scale; breeding-pair models on the log scale.",
    "All Rhat = 1.00; effective sample sizes > 2000."))


#=============================================================================
# OUTPUT
#=============================================================================

fig_population
fig_livestock
fig_grazing_habitat

tab1
tab2
tab3


#=============================================================================
# CONVERGENCE CHECK
#=============================================================================
# Rhat compares variation between chains to variation within them. Values at
# or very near 1.00 mean the chains agree and the posterior can be trusted.
# Anything above about 1.01 means the model needs re-running.
check_convergence <- function(model, model_name) {
  model_summary <- summary(model)
  rhat_values <- model_summary$fixed[, "Rhat"]
  cat(sprintf("[%-24s] max Rhat = %.3f\n", model_name, max(rhat_values, na.rm = TRUE)))
}

cat("\n== Convergence ==\n")
check_convergence(model_colony_trend,     "1a colonies ~ year")
check_convergence(model_pair_trend,       "1b pairs ~ year")
check_convergence(model_grazing_colonies, "3a colonies ~ grazing")
check_convergence(model_grazing_pairs,    "3b pairs ~ grazing")
check_convergence(model_habitat_colonies, "4a colonies ~ habitat")
check_convergence(model_habitat_pairs,    "4b pairs ~ habitat")


#=============================================================================
# QUICK CHECKS - observed peaks, for cross-checking the text
#=============================================================================
cat("\n== Observed peaks ==\n")
print(annual[which.max(annual$n_colonies),  c("Year", "n_colonies")])   # peak colony year
print(annual[order(-annual$n_colonies), c("Year", "n_colonies")][1:5, ]) # top five years
print(annual[which.max(annual$total_pairs), c("Year", "total_pairs")])  # peak pair year


#=============================================================================
# All analyses in R (R Core Team 2025). Spatial visualisation in QGIS 3.38.0.
#=============================================================================
sessionInfo()
