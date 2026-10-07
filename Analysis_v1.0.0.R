# ANALYSIS: data preparation, outlier removal, descriptive statistics,
# hypothesis testing and feature engineering.
# Please run this file first as it saves the packet-level features to
# data/bfm_packet_features.rds, which Modelling_v1.0.0.R reads.

library(lubridate)
library(ggplot2)
library(moments)
library(dunn.test)
library(car)
library(dplyr)

# Load the dataset
df <- read.csv("bfm_data.csv")

# Total number of rows and columns in dataset
nrow(df)
length(df)

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

# Save packet-level features (before outlier removal) for Modelling_v1.0.0.R
dir.create("data", showWarnings = FALSE)
saveRDS(df, "data/bfm_packet_features.rds")

df_foil <- df[df$environment == "foil", ]
df_nofoil <- df[df$environment == "nofoil", ]
df_open <- df[df$environment == "open", ]

# Number and percentage of data per environment
sum(df$environment == "foil")
round(sum(df$environment == "foil") / nrow(df) * 100, 2)
sum(df$environment == "nofoil")
round(sum(df$environment == "nofoil") / nrow(df) * 100, 2)
sum(df$environment == "open")
round(sum(df$environment == "open") / nrow(df) * 100, 2)

# Number and percentage of data per activity
sum(df$activity == "standing")
round(sum(df$activity == "standing") / nrow(df) * 100, 2)
sum(df$activity == "walking")
round(sum(df$activity == "walking") / nrow(df) * 100, 2)

# Number and percentage of data per subject
sum(df$subject == "abel")
round(sum(df$subject == "abel") / nrow(df) * 100, 2)
sum(df$subject == "collin")
round(sum(df$subject == "collin") / nrow(df) * 100, 2)
sum(df$subject == "ivan")
round(sum(df$subject == "ivan") / nrow(df) * 100, 2)
sum(df$subject == "kenny")
round(sum(df$subject == "kenny") / nrow(df) * 100, 2)
sum(df$subject == "matthew")
round(sum(df$subject == "matthew") / nrow(df) * 100, 2)

# Table of number of sessions per environment, subject and activity (raw data, before outlier removal)
sessions_per_subject <- df %>%
  group_by(environment, subject, activity) %>%
  summarise(n_sessions = n_distinct(session_id), .groups = "drop") %>%
  arrange(factor(environment, levels = c("open", "foil", "nofoil")), subject, activity)

dir.create("results", showWarnings = FALSE)
write.csv(sessions_per_subject, "results/sessions_per_subject.csv", row.names = FALSE)
print(as.data.frame(sessions_per_subject), row.names = FALSE)

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
  
  cat("Lower fence:", lower_fence, "\n")
  cat("Number of values below lower fence:", sum(x < lower_fence, na.rm = TRUE), "\n")
  cat("Percentage of values dropped:", round(sum(x < lower_fence, na.rm = TRUE) / length(x) * 100, 2), "\n")
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
  ylab = "Count"
)

# Legend drawn after the bars so par("usr") holds this plot's coordinates
legend(
  x = par("usr")[2] + 0.3,
  y = par("usr")[4],
  legend = rownames(bar_data),
  fill = c("steelblue", "tomato"),
  xpd = TRUE,
  bty = "n"
)
dev.off()

# SESSION-LEVEL AGGREGATION
# Collapse each session to one row (mean of each feature).
df_clean <- rbind(df_open_clean, df_foil_clean, df_nofoil_clean)
environments <- c("open", "foil", "nofoil")

session_level <- df_clean %>%
  group_by(session_id, environment, subject) %>%
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


# HYPOTHESIS TESTING (session-level, per subject)
# Kruskal-Wallis is used since not all features are normally distributed
# and they do not have equal variances
# Sessions from the same subject are not independent, so the environments
# are compared within each subject separately (10 sessions per environment
# per subject) instead of pooling all subjects into one test.
# H0: feature distribution is the same across environments for this subject.
# H1: at least one environment differs for this subject.
# Holm correction is applied across the subjects of each feature.
# Post-hoc pairwise comparisons (Dunn's, Bonferroni) are run per subject.

session_subjects <- sort(unique(session_level$subject))
kw_subject_results <- data.frame()
dunn_subject_results <- data.frame()

for (fidx in seq_along(session_features)) {
  feat <- session_features[fidx]
  for (subj in session_subjects) {
    sub <- session_level[session_level$subject == subj, ]

    # Kruskal-Wallis
    kw <- kruskal.test(sub[[feat]], factor(sub$environment, levels = environments))
    n <- nrow(sub)
    kw_subject_results <- rbind(kw_subject_results, data.frame(
      feature = feature_labels[fidx],
      subject = subj,
      n_sessions = n,
      chi_squared = round(unname(kw$statistic), 3),
      df = unname(kw$parameter),
      p_value = kw$p.value,
      # Epsilon-squared effect size: H / (n - 1)
      epsilon_sq = round(unname(kw$statistic) / (n - 1), 3)
    ))

    # Dunn's post-hoc pairwise comparisons (Bonferroni, two-sided p-values)
    invisible(capture.output(
      dn <- dunn.test(sub[[feat]], sub$environment, method = "bonferroni",
                      kw = FALSE, table = FALSE, altp = TRUE)
    ))
    dunn_subject_results <- rbind(dunn_subject_results, data.frame(
      feature = feature_labels[fidx],
      subject = subj,
      comparison = dn$comparisons,
      z = round(dn$Z, 3),
      p_adj_bonferroni = round(dn$altP.adjusted, 4)
    ))
  }
}

kw_subject_results <- kw_subject_results %>%
  group_by(feature) %>%
  mutate(p_adj_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup() %>%
  mutate(significant = ifelse(p_adj_holm < 0.05, "Yes", "No"),
         p_value = round(p_value, 4),
         p_adj_holm = round(p_adj_holm, 4))

# Number of subjects in which each feature differs across environments
kw_subject_summary <- kw_subject_results %>%
  group_by(feature) %>%
  summarise(n_subjects_significant = sum(significant == "Yes"),
            n_subjects = n(),
            median_epsilon_sq = median(epsilon_sq),
            .groups = "drop")

cat("\n=== Kruskal-Wallis Hypothesis Tests (session-level, per subject) ===\n")
print(as.data.frame(kw_subject_results), row.names = FALSE)
cat("\n=== Dunn's Post-hoc Pairwise Tests (per subject, Bonferroni) ===\n")
print(as.data.frame(dunn_subject_results), row.names = FALSE)
cat("\n=== Subjects with significant environment effect (Holm-adjusted p < 0.05) ===\n")
print(as.data.frame(kw_subject_summary), row.names = FALSE)

write.csv(kw_subject_results, "results/kruskal_wallis_per_subject.csv", row.names = FALSE)
write.csv(dunn_subject_results, "results/dunn_per_subject.csv", row.names = FALSE)

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

# Function to calculate variance of rate of change of features
var_roc <- function(x, t) {
  dt <- as.numeric(diff(t), units = "secs")
  var(diff(x) / dt, na.rm = TRUE)
}

raw_features <- c("Mean_Magnitude", "Mean_Phase", "Std_Magnitude", "Phase_Coherence")

session_engineered <- df_clean %>%
  group_by(session_id, environment, activity) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  summarise(
    # Variance of rate of change
    var_roc_mean_mag = var_roc(Mean_Magnitude, timestamp),
    var_roc_mean_pha = var_roc(Mean_Phase, timestamp),
    var_roc_std_mag  = var_roc(Std_Magnitude, timestamp),
    var_roc_pha_coh  = var_roc(Phase_Coherence, timestamp),
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
      var_roc_mean_mag = var_roc(Mean_Magnitude, timestamp),
      var_roc_mean_pha = var_roc(Mean_Phase, timestamp),
      var_roc_std_mag  = var_roc(Std_Magnitude, timestamp),
      var_roc_pha_coh  = var_roc(Phase_Coherence, timestamp),
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
