library(lubridate)
library(ggplot2)
library(moments)
library(boot)
library(dunn.test)
library(car)
library(dplyr)
library(nlme)
library(caret)
library(rpart)
library(randomForest)
library(e1071)
library(InformationValue)
# Load the dataset
df <- read.csv("bfm_data.csv")

# Check for missing values
colSums(is.na(df))

# Convert timestamp column from UTC to Asia/Kuala_Lumpur (MYT)
t_utc <- as.POSIXct(df$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
t_myt <- with_tz(t_utc, "Asia/Kuala_Lumpur")
df$timestamp <- format(t_myt, "%Y-%m-%d %H:%M:%OS6+08:00")

# Identify Magnitude, Phase and Address columns
scidx_cols <- grep("^SCIDX", colnames(df), value = TRUE)
mag_cols <- grep("Mag$|Magnitude$", scidx_cols, value = TRUE)
phase_cols <- grep("Phase$|phase$", scidx_cols, value = TRUE)
address_cols <- grep("address$", colnames(df), value = TRUE)

# Calculate row-wise mean and standard deviation for Magnitude
X_mag <- as.matrix(df[, mag_cols])
Mean_Magnitude <- rowMeans(X_mag, na.rm = TRUE)
N_mag <- ncol(X_mag)
Std_Magnitude <- sqrt((rowSums(X_mag^2, na.rm = TRUE) - N_mag * Mean_Magnitude^2) / (N_mag - 1))

# Calculate row-wise mean, standard deviation, and phase coherence for Phase
X_phase <- as.matrix(df[, phase_cols])
Mean_Phase <- rowMeans(X_phase, na.rm = TRUE)
N_phase <- ncol(X_phase)

cos_mean <- rowMeans(cos(X_phase), na.rm = TRUE)
sin_mean <- rowMeans(sin(X_phase), na.rm = TRUE)
Phase_Coherence <- sqrt(cos_mean^2 + sin_mean^2)

# Bind calculated statistics to dataframe
df$Mean_Magnitude <- Mean_Magnitude
df$Std_Magnitude <- Std_Magnitude
df$Mean_Phase <- Mean_Phase
df$Phase_Coherence <- Phase_Coherence
df$timestamp <- as.POSIXct(df$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "Asia/Kuala_Lumpur")
df$timestamp <- ymd_hms(df$timestamp, tz = "Asia/Kuala_Lumpur")

cols_to_drop <- c(scidx_cols, address_cols)
df <- df[, !colnames(df) %in% cols_to_drop]

df_foil <- df[df$environment == "foil", ]
df_nofoil <- df[df$environment == "nofoil", ]
df_open <- df[df$environment == "open", ]
sum(df$environment == "foil") / nrow(df)
sum(df$environment == "nofoil") / nrow(df)
sum(df$environment == "open") / nrow(df)

# Function for Histogram of a Feature for an Environment
plot_histogram <- function(data, env, x_col, xlim = c(0, 30), ylim = c(0, 50), breaks = 30, color = "steelblue", by = 1) {
  hist(
    data[[x_col]],
    breaks = breaks,
    col = color,
    border = "black",
    main = paste("Distribution of", x_col, "(", env, ")"),
    xlab = x_col,
    ylab = "Count",
    xaxt = "n",
    xlim = xlim,
    ylim = ylim
  )
  
  axis(1, at = seq(xlim[1], xlim[2], by = by))
}

# Function to remove outliers using IQR method
remove_outliers_iqr <- function(data, x_col) {
  x <- data[[x_col]]
  
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_fence <- Q1 - 1.5 * IQR_val
  upper_fence <- Q3 + 1.5 * IQR_val
  
  cat("Lower fence:", lower_fence, "\n")
  cat("Number of values below lower fence:", sum(x < lower_fence, na.rm = TRUE), "\n")
  cat("Percentage of values dropped:", round(sum(x < lower_fence, na.rm = TRUE) / length(x) * 100, 2))
  cleaned_data <- data[x >= lower_fence, ]
  return(cleaned_data)
}

# Function to create statistical properties summary table
cols <- c("mean_mag", "std_mag", "mean_pha", "pha_coh")
describe_features <- function(data, cols) {
  stats_table <- data.frame(
    Feature = cols,
    Mean = sapply(cols, function(c) mean(data[[c]], na.rm = TRUE)),
    Median = sapply(cols, function(c) median(data[[c]], na.rm = TRUE)),
    Std_Dev = sapply(cols, function(c) sd(data[[c]], na.rm = TRUE)),
    Skewness = sapply(cols, function(c) skewness(data[[c]], na.rm = TRUE)),
    Kurtosis = sapply(cols, function(c) kurtosis(data[[c]], na.rm = TRUE)),
    row.names = NULL
  )
  return(stats_table)
}

# Directory for all saved figures
fig_dir <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Histogram Plot of Features in Open Environment Before and After Data Cleaning
png(file.path(fig_dir, "hist_open_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open, "Open", "Mean_Magnitude", xlim = c(0,22), ylim = c(0, 35000))
dev.off()
png(file.path(fig_dir, "hist_open_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open, "Open", "Std_Magnitude", ylim = c(0, 35000))
dev.off()
png(file.path(fig_dir, "hist_open_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open, "Open", "Mean_Phase", xlim = c(-5, 4), ylim = c(0, 40000))
dev.off()
png(file.path(fig_dir, "hist_open_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open, "Open", "Phase_Coherence", xlim = c(0.6,1), ylim = c(0, 35000), by = 0.1)
dev.off()

df_open_clean <- remove_outliers_iqr(df_open, "Mean_Magnitude")

png(file.path(fig_dir, "hist_open_clean_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open_clean, "Cleaned Open", "Mean_Magnitude", xlim = c(10, 22), ylim = c(0, 35000))
dev.off()
png(file.path(fig_dir, "hist_open_clean_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open_clean, "Cleaned Open", "Std_Magnitude", ylim = c(0, 35000))
dev.off()
png(file.path(fig_dir, "hist_open_clean_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open_clean, "Cleaned Open", "Mean_Phase",xlim = c(0.37, 0.42), ylim = c(0, 10000), by = 0.01)
dev.off()
png(file.path(fig_dir, "hist_open_clean_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_open_clean, "Cleaned Open", "Phase_Coherence", xlim = c(0.98,1), ylim = c(0, 25000), by = 0.001)
dev.off()

# Histogram Plot of Features in Foil Environment Before and After Data Cleaning
png(file.path(fig_dir, "hist_foil_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil, "Foil", "Mean_Magnitude", ylim = c(0, 6000))
dev.off()
png(file.path(fig_dir, "hist_foil_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil, "Foil", "Std_Magnitude", xlim = c(0, 40), ylim = c(0, 6000))
dev.off()
png(file.path(fig_dir, "hist_foil_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil, "Foil", "Mean_Phase", xlim = c(-15, 15), ylim = c(0, 50000))
dev.off()
png(file.path(fig_dir, "hist_foil_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil, "Foil", "Phase_Coherence", xlim = c(0,1), ylim = c(0, 50000))
dev.off()

df_foil_clean <- remove_outliers_iqr(df_foil, "Mean_Magnitude")

png(file.path(fig_dir, "hist_foil_clean_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil_clean, "Cleaned Foil", "Mean_Magnitude", ylim = c(0, 4000))
dev.off()
png(file.path(fig_dir, "hist_foil_clean_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil_clean, "Cleaned Foil", "Std_Magnitude", xlim = c(0, 40), ylim = c(0, 6000))
dev.off()
png(file.path(fig_dir, "hist_foil_clean_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil_clean, "Cleaned Foil", "Mean_Phase",xlim = c(0.26, 0.52), ylim = c(0, 5000), by = 0.01)
dev.off()
png(file.path(fig_dir, "hist_foil_clean_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_foil_clean, "Cleaned Foil", "Phase_Coherence", xlim = c(0.90,1), ylim = c(0, 5000), by = 0.05)
dev.off()

# Histogram Plot of Features in No Foil Environment Before and After Data Cleaning
png(file.path(fig_dir, "hist_nofoil_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil, "No Foil", "Mean_Magnitude", ylim = c(0, 15000))
dev.off()
png(file.path(fig_dir, "hist_nofoil_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil, "No Foil", "Std_Magnitude", xlim = c(0, 40), ylim = c(0, 12000))
dev.off()
png(file.path(fig_dir, "hist_nofoil_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil, "No Foil", "Mean_Phase", xlim = c(-10, 10), ylim = c(0, 50000))
dev.off()
png(file.path(fig_dir, "hist_nofoil_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil, "No Foil", "Phase_Coherence", xlim = c(0,1), ylim = c(0, 50000))
dev.off()

df_nofoil_clean <- remove_outliers_iqr(df_nofoil, "Mean_Magnitude")

png(file.path(fig_dir, "hist_nofoil_clean_Mean_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil_clean, "Cleaned No Foil", "Mean_Magnitude", ylim = c(0, 10000))
dev.off()
png(file.path(fig_dir, "hist_nofoil_clean_Std_Magnitude.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil_clean, "Cleaned No Foil", "Std_Magnitude", xlim = c(0, 40), ylim = c(0, 12000))
dev.off()
png(file.path(fig_dir, "hist_nofoil_clean_Mean_Phase.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil_clean, "Cleaned No Foil", "Mean_Phase",xlim = c(0.31, 0.50), ylim = c(0, 10000), by = 0.01)
dev.off()
png(file.path(fig_dir, "hist_nofoil_clean_Phase_Coherence.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
plot_histogram(df_nofoil_clean, "Cleaned No Foil", "Phase_Coherence", xlim = c(0.93,1), ylim = c(0, 15000), by = 0.01)
dev.off()

# Desciptive Analysis
nrow(df_open_clean)
nrow(df_foil_clean)
nrow(df_nofoil_clean)

# Categorical Variable Distribution
sum(df_open_clean$activity == "standing") / nrow(df_open_clean)
sum(df_open_clean$activity == "walking") / nrow(df_open_clean)

sum(df_foil_clean$activity == "standing")
sum(df_foil_clean$activity == "walking")
sum(df_foil_clean$activity == "standing") / nrow(df_foil_clean)
sum(df_foil_clean$activity == "walking") / nrow(df_foil_clean)

sum(df_nofoil_clean$activity == "standing")
sum(df_nofoil_clean$activity == "walking")
sum(df_nofoil_clean$activity == "standing") / nrow(df_nofoil_clean)
sum(df_nofoil_clean$activity == "walking") / nrow(df_nofoil_clean)

standing_counts <- c(
  sum(df_open_clean$activity == "standing", na.rm = TRUE),
  sum(df_foil_clean$activity == "standing", na.rm = TRUE),
  sum(df_nofoil_clean$activity == "standing", na.rm = TRUE)
)

walking_counts <- c(
  sum(df_open_clean$activity == "walking", na.rm = TRUE),
  sum(df_foil_clean$activity == "walking", na.rm = TRUE),
  sum(df_nofoil_clean$activity == "walking", na.rm = TRUE)
)

# Group Bar Chart for Activity Count per Environment
bar_data <- rbind(standing_counts, walking_counts)
colnames(bar_data) <- c("Open", "Foil", "No Foil")
rownames(bar_data) <- c("Standing", "Walking")

png(file.path(fig_dir, "barplot_activity_counts.png"), width = 8, height = 6, units = "in", res = 300, pointsize = 16)
par(mar = c(5, 4, 4, 8), xpd = TRUE)   # expand right margin to draw legend

barplot(
  bar_data,
  beside = TRUE,
  col = c("steelblue", "tomato"),
  main = "Activity Counts Across Environments",
  xlab = "Environment",
  ylab = "Count",
  legend.text = rownames(bar_data),
  args.legend = list(
    x = par("usr")[2] + 1.5,
    y = par("usr")[4],
    xpd = TRUE,
    bty = "n"
  )
)
dev.off()

# SESSION-LEVEL AGGREGATION
# Collapse each session to one row (mean of each feature).
df_clean <- rbind(df_open_clean, df_foil_clean, df_nofoil_clean)
environments <- c("open", "foil", "nofoil")

session_level <- df_clean %>%
  group_by(session_id, environment) %>%
  summarise(
    mean_mag = mean(Mean_Magnitude, na.rm = TRUE),
    mean_pha = mean(Mean_Phase, na.rm = TRUE),
    std_mag = mean(Std_Magnitude, na.rm = TRUE),
    pha_coh = mean(Phase_Coherence, na.rm = TRUE),
    .groups = "drop"
  )

# Session counts per environment (sanity check)
table(session_level$environment)

# Descriptive Statistics Summary Table (Open)
stats_table_open <- describe_features(session_level[session_level$environment == "open", ], cols)
stats_table_open

# Descriptive Statistics Summary Table (Foil)
stats_table_foil <- describe_features(session_level[session_level$environment == "foil", ], cols)
stats_table_foil

# Descriptive Statistics Summary Table (No Foil)
stats_table_nofoil <- describe_features(session_level[session_level$environment == "nofoil", ], cols)
stats_table_nofoil

session_features <- c("mean_mag", "mean_pha", "std_mag", "pha_coh")
feature_labels <- c("Mean_Magnitude", "Mean_Phase", "Std_Magnitude", "Phase_Coherence")


# ASSUMPTION CHECKING
# Normality (Shapiro-Wilk Test, per environment)
# H0: data are normally distributed.
# H1: data are not normally distributed

cat("\n=== Shapiro-Wilk Normality Tests (session-level, per environment) ===\n")
for (fidx in seq_along(session_features)) {
  feat <- session_features[fidx]
  cat("\n-- Feature:", feature_labels[fidx], "--\n")
  for (env in environments) {
    x <- session_level[session_level$environment == env, feat, drop = TRUE]
    sw <- shapiro.test(x)
    cat(sprintf("  %-8s : W = %.4f, p = %.4f %s\n",
                env, sw$statistic, sw$p.value,
                ifelse(sw$p.value > 0.05, "(normal)", "(NON-normal)")))
  }
}

# Levene's test (Equal of Variance)
# H0: variances are equal across environments.
# H1: at least one environment pair does not have equal variances
# Required assumption for standard one-way ANOVA.

cat("\n=== Levene's Test for Homogeneity of Variance (session-level) ===\n")
for (fidx in seq_along(session_features)) {
  feat <- session_features[fidx]
  cat("\n-- Feature:", feature_labels[fidx], "--\n")
  print(leveneTest(as.formula(paste(feat, "~ environment")),
                     data = session_level))
}


# HYPOTHESIS TESTING (session-level)
# Kruskal-Wallis is used since not all features are normally distributed
# and they do not have equal variances
# Post-hoc pairwise comparisons are run for all features

cat("\n=== Kruskal-Wallis Hypothesis Tests (session-level) ===\n")
for (fidx in seq_along(session_features)) {
  feat <- session_features[fidx]
  cat("\n---------- Feature:", feature_labels[fidx], "----------\n")
  
  # Kruskal-Wallis
  cat("[Kruskal-Wallis]\n")
  print(kruskal.test(as.formula(paste(feat, "~ environment")),
                     data = session_level))
  
  # Dunn's post-hoc pairwise comparisons (Bonferroni)
  cat("[Dunn's post-hoc pairwise (Bonferroni)]\n")
  dunn.test(session_level[[feat]], session_level$environment,
            method = "bonferroni")
}

# STANDING vs WALKING COMPARISON (per environment, session-level)
session_activity <- df_clean %>%
  group_by(session_id, environment, activity) %>%
  summarise(
    mean_mag = mean(Mean_Magnitude, na.rm = TRUE),
    mean_pha = mean(Mean_Phase, na.rm = TRUE),
    std_mag = mean(Std_Magnitude, na.rm = TRUE),
    pha_coh = mean(Phase_Coherence, na.rm = TRUE),
    .groups = "drop"
  )

# Descriptive statistics: standing vs walking, per environment
# Violin plots + Box Plots: standing vs walking per feature per environment
for (fidx in seq_along(session_features)) {
  feat <- session_features[fidx]
  p <- ggplot(session_activity,
              aes(x = activity, y = .data[[feat]], fill = activity)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 fill = "white", alpha = 0.6) +
    stat_summary(fun = mean, geom = "point",
                 shape = 18, size = 3, color = "black") +
    facet_wrap(~ environment, nrow = 1) +
    scale_fill_manual(values = c(standing = "steelblue", walking = "tomato")) +
    labs(
      title = paste("Standing vs Walking -", feature_labels[fidx],
                    "(session-level, by environment)"),
      x = "Activity", y = feature_labels[fidx]
    ) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, face = "bold"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          strip.text = element_text(size = 20, face = "bold"))
  print(p)
  ggsave(file.path(fig_dir, paste0("violin_activity_", feat, ".png")),
         plot = p, width = 14, height = 8, dpi = 300)
}

# Hypothesis tests (per environment)
cat("\n=== Standing vs Walking: Hypothesis Tests ===\n")
for (env in environments) {
  cat("\n========== Environment:", toupper(env), "==========\n")
  sub <- session_activity[session_activity$environment == env, ]
  
  for (fidx in seq_along(session_features)) {
    feat <- session_features[fidx]
    cat("\n -- Feature:", feature_labels[fidx], "--\n")
    
    x_stand <- sub[sub$activity == "standing", feat, drop = TRUE]
    x_walk <- sub[sub$activity == "walking",  feat, drop = TRUE]
    
    # Wilcoxon rank-sum test
    wx <- wilcox.test(x_stand, x_walk, exact = FALSE)
    cat(sprintf("Wilcoxon rank-sum     : W = %.0f, p = %.4f\n",
                wx$statistic, wx$p.value))
    if (wx$p.value < 0.05) {
      cat("Reject null hypothesis")
    } else {
      cat("Fail to reject null hypothesis")
    }
  }
}


# TIME SERIES PLOTS — 1-MINUTE MIDDLE WINDOW PER SESSION

# Build the windowed dataset across all sessions
window_duration <- 60   # total window in seconds

ts_window <- df_clean %>%
  group_by(session_id, environment, activity) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    session_start = min(timestamp),
    session_end = max(timestamp),
    session_mid = session_start + as.numeric(session_end - session_start,
                                                units = "secs") / 2,
    window_start = session_mid - window_duration / 2,
    window_end = session_mid + window_duration / 2
  ) %>%
  filter(timestamp >= window_start & timestamp <= window_end) %>%
  mutate(
    # Relative time in seconds from the start of this session's window
    rel_time = as.numeric(timestamp - window_start, units = "secs")
  ) %>%
  ungroup()

# Time Series Plot for Each Feature
ts_features <- c("Mean_Magnitude", "Mean_Phase", "Std_Magnitude", "Phase_Coherence")

for (feat in ts_features) {
  p <- ggplot(ts_window,
              aes(x = rel_time,
                  y = .data[[feat]],
                  group = session_id,
                  color = activity)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    facet_wrap(~ environment, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c(standing = "steelblue", walking = "tomato")) +
    scale_x_continuous(
      breaks = seq(0, window_duration, by = 10),
      labels = seq(0, window_duration, by = 10)
    ) +
    labs(
      title = paste("Time Series —", feat,
                     "(1-min middle window, all sessions)"),
      x = "Relative time (seconds)",
      y = feat,
      color = "Activity"
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 4, alpha = 1))) +
    theme_minimal(base_size = 18) +
    theme(
      legend.key.width = unit(2, "cm"),
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  print(p)
  ggsave(file.path(fig_dir, paste0("timeseries_", feat, ".png")),
         plot = p, width = 15, height = 8, dpi = 300)
}

# FEATURE ENGINEERING — SESSION-LEVEL TIME SERIES FEATURES
raw_features <- c("Mean_Magnitude", "Mean_Phase", "Std_Magnitude", "Phase_Coherence")

session_engineered <- df_clean %>%
  group_by(session_id, environment, activity) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  summarise(
    # Variance of rate of change
    var_roc_mean_mag = {
      dt <- as.numeric(diff(timestamp), units = "secs")
      dx <- diff(Mean_Magnitude)
      var(dx / dt, na.rm = TRUE)
    },
    var_roc_mean_pha = {
      dt <- as.numeric(diff(timestamp), units = "secs")
      dx <- diff(Mean_Phase)
      var(dx / dt, na.rm = TRUE)
    },
    var_roc_std_mag  = {
      dt <- as.numeric(diff(timestamp), units = "secs")
      dx <- diff(Std_Magnitude)
      var(dx / dt, na.rm = TRUE)
    },
    var_roc_pha_coh  = {
      dt <- as.numeric(diff(timestamp), units = "secs")
      dx <- diff(Phase_Coherence)
      var(dx / dt, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  mutate(activity_num = ifelse(activity == "walking", 1, 0))

# Descriptive statistics for engineered features
eng_features <- c("var_roc_mean_mag", "var_roc_mean_pha", "var_roc_std_mag", "var_roc_pha_coh")
eng_labels   <- c("Var_Roc(Mean_Mag)", "Var_Roc(Mean_Phase)", "Var_Roc(Std_Mag)", "Var_Roc(Phase_Coh)")

for (env in environments) {
  session_engineered_env <- df_clean %>%
    filter(environment == env) %>%
    group_by(session_id, environment, activity) %>%
    arrange(timestamp, .by_group = TRUE) %>%
    summarise(
      var_roc_mean_mag = {
        dt <- as.numeric(diff(timestamp), units = "secs")
        dx <- diff(Mean_Magnitude)
        var(dx / dt, na.rm = TRUE)
      },
      var_roc_mean_pha = {
        dt <- as.numeric(diff(timestamp), units = "secs")
        dx <- diff(Mean_Phase)
        var(dx / dt, na.rm = TRUE)
      },
      var_roc_std_mag  = {
        dt <- as.numeric(diff(timestamp), units = "secs")
        dx <- diff(Std_Magnitude)
        var(dx / dt, na.rm = TRUE)
      },
      var_roc_pha_coh  = {
        dt <- as.numeric(diff(timestamp), units = "secs")
        dx <- diff(Phase_Coherence)
        var(dx / dt, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  for (fidx in seq_along(eng_features)) {
    feat <- eng_features[fidx]
    p <- ggplot(session_engineered_env,
                aes(x = activity, y = .data[[feat]], fill = activity)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.shape = NA,
                   fill = "white", alpha = 0.6) +
      # Plot mean point
      stat_summary(fun = mean, geom = "point",
                   shape = 18, size = 3, color = "black") +
      facet_wrap(~ environment, nrow = 1) +
      scale_fill_manual(values = c(standing = "steelblue", walking = "tomato")) +
      labs(
        title = paste("Standing vs Walking —", eng_labels[fidx],
                      "(session-level)"),
        x = "Activity", y = eng_labels[fidx]
      ) +
      theme_minimal(base_size = 18) +
      theme(legend.position = "none",
            plot.title = element_text(size = 18, face = "bold"),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 18, face = "bold"))
    print(p)
    ggsave(file.path(fig_dir, paste0("violin_engineered_", env, "_", feat, ".png")),
           plot = p, width = 11, height = 8, dpi = 300)
  }
}

# Assumption checks + hypothesis tests for engineered features
cat("\n=== Standing vs Walking: Tests on Engineered Features ===\n")

for (env in environments) {
  cat("\n========== Environment:", toupper(env), "==========\n")
  sub <- session_engineered[session_engineered$environment == env, ]
  
  for (fidx in seq_along(eng_features)) {
    feat <- eng_features[fidx]
    cat("\n  --", eng_labels[fidx], "--\n")
    
    x_stand <- sub[sub$activity == "standing", feat, drop = TRUE]
    x_walk  <- sub[sub$activity == "walking",  feat, drop = TRUE]
    
    sw_s <- shapiro.test(x_stand)
    sw_w <- shapiro.test(x_walk)
    
    cat(sprintf("Shapiro-Wilk standing : W = %.4f, p = %.4f %s\n",
                sw_s$statistic, sw_s$p.value,
                ifelse(sw_s$p.value > 0.05, "(normal)", "(NON-normal)")))
    cat(sprintf("Shapiro-Wilk walking  : W = %.4f, p = %.4f %s\n",
                sw_w$statistic, sw_w$p.value,
                ifelse(sw_w$p.value > 0.05, "(normal)", "(NON-normal)")))
    
    both_normal <- sw_s$p.value > 0.05 && sw_w$p.value > 0.05
    
    if (both_normal) {
      tt <- t.test(x_stand, x_walk, var.equal = FALSE)
      cat(sprintf("Welch's t-test        : t = %.4f, df = %.2f, p = %.4f\n",
                  tt$statistic, tt$parameter, tt$p.value))
    } else {
      wx <- wilcox.test(x_stand, x_walk, exact = FALSE)
      cat(sprintf("Wilcoxon rank-sum     : W = %.0f, p = %.4f\n",
                  wx$statistic, wx$p.value))
    }
  }
}

# Variance Inflation Factor (VIF)
vif_model <- lm(activity_num ~ var_roc_mean_mag + var_roc_mean_pha + var_roc_std_mag + var_roc_pha_coh, data = session_engineered)

vif_values <- vif(vif_model)
print(vif_values)

# Compute and plot Spearman correlation matrix per environment
for (env in environments) {
  sub <- session_engineered[session_engineered$environment == env, ]
  
  cor_matrix <- cor(
    sub[, c("var_roc_mean_mag", "var_roc_mean_pha", "var_roc_std_mag", "var_roc_pha_coh", "activity_num")],
    method = "spearman"
  )
  
  rownames(cor_matrix) <- colnames(cor_matrix) <-
    c("var_roc(Mean_Mag)", "var_roc(Mean_Phase)", "var_roc(Std_Mag)", "var_roc(Phase_Coh)", "Activity")
  
  cor_long <- as.data.frame(as.table(cor_matrix))
  names(cor_long) <- c("Var1", "Var2", "Correlation")
  
  p <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)),
              size = 6, fontface = "bold") +
    scale_fill_gradient2(
      low = "steelblue",
      mid = "white",
      high = "tomato",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(
      title = paste("Spearman Correlation Heatmap —", toupper(env), "environment"),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      panel.grid = element_blank()
    )
  print(p)
  ggsave(file.path(fig_dir, paste0("spearman_heatmap_", env, ".png")),
         plot = p, width = 11, height = 8, dpi = 300)
}


# Build session-level features (unscaled) for each environment.
# Min-max scaling is applied per LOSO fold below, using the training
# set's own min/max, so it happens AFTER the train/test split.
build_session_data <- function(env) {
  df_clean %>%
    filter(environment == env) %>%
    group_by(session_id, activity, subject) %>%
    arrange(timestamp, .by_group = TRUE) %>%
    summarise(
      var_roc_pha_coh = {
        dt <- as.numeric(diff(timestamp), units = "secs")
        dx <- diff(Phase_Coherence)
        var(dx / dt, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(
      activity_num    = ifelse(activity == "walking", 1, 0),
      activity_factor = factor(activity, levels = c("standing", "walking"))
    )
}

session_data_open   <- build_session_data("open")
session_data_foil   <- build_session_data("foil")
session_data_nofoil <- build_session_data("nofoil")

session_data_by_env <- list(
  "Open"    = session_data_open,
  "Foil"    = session_data_foil,
  "No Foil" = session_data_nofoil
)

subjects <- c("collin", "kenny", "abel", "matthew", "ivan")


# CROSS-ENVIRONMENT LOSO EVALUATION — LOGISTIC REGRESSION, DECISION TREE, RANDOM FOREST, SVM
# For each model and each ordered pair of (train environment, test environment),
# run leave-one-subject-out CV: train on all subjects but one in the train
# environment (min-max scaled on the training fold only), test on the held-out
# subject's data from the test environment (scaled with the training fold's
# min/max), then average Accuracy / Sensitivity / Specificity across subjects.

# Fit a classifier of the given type on a (already scaled) training fold
fit_model <- function(model_type, train_data) {
  switch(model_type,
    "lr"   = glm(activity_num ~ var_roc_pha_coh_scaled, data = train_data, family = "binomial"),
    "tree" = rpart(activity_factor ~ var_roc_pha_coh_scaled, data = train_data, method = "class"),
    "rf"   = randomForest(activity_factor ~ var_roc_pha_coh_scaled, data = train_data, ntree = 100),
    "svm"  = svm(activity_factor ~ var_roc_pha_coh_scaled, data = train_data, kernel = "radial", probability = TRUE)
  )
}

# Get predicted probability of the "walking" class for a fitted classifier
predict_walking_prob <- function(model, model_type, newdata) {
  switch(model_type,
    "lr"   = predict(model, newdata = newdata, type = "response"),
    "tree" = predict(model, newdata = newdata, type = "prob")[, "walking"],
    "rf"   = predict(model, newdata = newdata, type = "prob")[, "walking"],
    "svm"  = attr(predict(model, newdata = newdata, probability = TRUE), "probabilities")[, "walking"]
  )
}

model_types <- c("Logistic Regression" = "lr", "Decision Tree" = "tree",
                  "Random Forest" = "rf", "SVM" = "svm")
env_labels <- c("Open", "Foil", "No Foil")

cross_env_loso_results <- data.frame()

for (model_name in names(model_types)) {
  model_type <- model_types[[model_name]]

  for (train_label in env_labels) {
    for (test_label in env_labels) {
      fold_metrics <- data.frame()

      for (test_subject in subjects) {
        train_fold <- session_data_by_env[[train_label]]
        train_fold <- train_fold[train_fold$subject != test_subject, ]
        test_fold  <- session_data_by_env[[test_label]]
        test_fold  <- test_fold[test_fold$subject == test_subject, ]

        if (nrow(train_fold) == 0 || nrow(test_fold) == 0) next
        if (length(unique(train_fold$activity_factor)) < 2) next

        actual_factor <- factor(test_fold$activity_factor, levels = c("standing", "walking"))
        if (length(unique(actual_factor)) < 2) next

        # Min-max scale using the training fold's own range only
        train_min <- min(train_fold$var_roc_pha_coh, na.rm = TRUE)
        train_max <- max(train_fold$var_roc_pha_coh, na.rm = TRUE)
        train_fold$var_roc_pha_coh_scaled <- (train_fold$var_roc_pha_coh - train_min) / (train_max - train_min)
        test_fold$var_roc_pha_coh_scaled  <- (test_fold$var_roc_pha_coh - train_min) / (train_max - train_min)

        model <- fit_model(model_type, train_fold)
        probs <- predict_walking_prob(model, model_type, test_fold)
        pred  <- factor(ifelse(probs > 0.5, "walking", "standing"), levels = c("standing", "walking"))
        cm    <- caret::confusionMatrix(pred, actual_factor, positive = "walking")

        fold_metrics <- rbind(fold_metrics, data.frame(
          accuracy    = as.numeric(cm$overall["Accuracy"]),
          sensitivity = as.numeric(cm$byClass["Sensitivity"]),
          specificity = as.numeric(cm$byClass["Specificity"])
        ))
      }

      mean_accuracy    <- round(mean(fold_metrics$accuracy, na.rm = TRUE), 4)
      mean_sensitivity <- round(mean(fold_metrics$sensitivity, na.rm = TRUE), 4)
      mean_specificity <- round(mean(fold_metrics$specificity, na.rm = TRUE), 4)

      cat(sprintf("\n%s trained in %s Environment (LOSO), tested with %s environment data (%d/%d folds):\n",
                  model_name, train_label, test_label, nrow(fold_metrics), length(subjects)))
      cat("Mean Accuracy: ", mean_accuracy,
          " Mean Sensitivity: ", mean_sensitivity,
          " Mean Specificity: ", mean_specificity, "\n")

      cross_env_loso_results <- rbind(cross_env_loso_results, data.frame(
        model             = model_name,
        train_environment = train_label,
        test_environment  = test_label,
        accuracy          = mean_accuracy,
        sensitivity       = mean_sensitivity,
        specificity       = mean_specificity,
        n_folds           = nrow(fold_metrics)
      ))
    }
  }
}

print(cross_env_loso_results)

write.csv(cross_env_loso_results, "results/loso_results.csv", row.names = FALSE)
