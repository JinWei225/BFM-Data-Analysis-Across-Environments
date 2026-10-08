# MODELLING: leave-one-subject-out (LOSO) classification of standing vs walking.
# Run Analysis_v1.0.0.R first to create data/bfm_packet_features.rds.

library(dplyr)
library(ggplot2)
library(caret)
library(rpart)
library(randomForest)
library(e1071)

# Load packet-level features saved by Analysis_v1.0.0.R
features_file <- "data/bfm_packet_features.rds"
if (!file.exists(features_file))
  stop("Run Analysis_v1.0.0.R first to create ", features_file)
df <- readRDS(features_file)

dir.create("results", showWarnings = FALSE)

# Function to calculate variance of rate of change of features
var_roc <- function(x, t) {
  dt <- as.numeric(diff(t), units = "secs")
  var(diff(x) / dt, na.rm = TRUE)
}

# Set random seed so that the results are reproducible
set.seed(42)

subjects <- sort(unique(df$subject))
envs <- c("open", "foil", "nofoil")
env_labels <- c(open = "Open", foil = "Foil", nofoil = "No Foil")
potential_features <- c("var_roc_mean_mag", "var_roc_mean_pha", 
                        "var_roc_std_mag", "var_roc_pha_coh")
model_types <- c("Logistic Regression" = "lr", "Decision Tree" = "tree",
                 "Random Forest" = "rf", "SVM" = "svm")

FENCE_MODE    <- "train-only"
FIXED_FEATURE <- NULL

lower_fence <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  unname(Q1 - 1.5 * IQR_val)
}

# Function to create session-level engineered features
build_session_features <- function(packets) {
  packets %>%
    group_by(session_id, environment, activity, subject) %>%
    arrange(timestamp, .by_group = TRUE) %>%
    summarise(
      var_roc_mean_mag = var_roc(Mean_Magnitude, timestamp),
      var_roc_mean_pha = var_roc(Mean_Phase, timestamp),
      var_roc_std_mag  = var_roc(Std_Magnitude, timestamp),
      var_roc_pha_coh  = var_roc(Phase_Coherence, timestamp),
      .groups = "drop"
    ) %>%
    mutate(
      activity_num    = ifelse(activity == "walking", 1, 0),
      activity_factor = factor(activity, levels = c("standing", "walking"))
    )
}

# Function to apply IQR-based filtering on each fold
prepare_env_fold <- function(env, test_subject) {
  env_packets <- df[df$environment == env, ]
  train_packets <- env_packets[env_packets$subject != test_subject, ]
  fence <- lower_fence(train_packets$Mean_Magnitude)
  keep <- env_packets$Mean_Magnitude >= fence
  removed <- env_packets$Mean_Magnitude < fence
  pct_removed <- tapply(removed, env_packets$subject, mean) * 100
  cat("Fence:", round(fence, 2), "| % of packets removed per subject:\n")
  print(round(pct_removed, 2))
  sessions <- build_session_features(env_packets[keep, ])
  return(list(sessions = sessions, fence = fence, removed = pct_removed))
}

# Feature(s) with the largest absolute Spearman correlation with activity
# All features tied at the maximum are returned so that they share the fold's vote
select_feature <- function(train_sessions) {
  rho <- sapply(potential_features, function(f)
    cor(train_sessions[[f]], train_sessions$activity_num, method = "spearman"))
  abs_rho <- abs(rho)
  tied <- names(abs_rho)[abs(abs_rho - max(abs_rho)) < 1e-9]
  list(feature = tied, rho = rho)
}

# Fit a classifier of the given type on a (already scaled) training fold
# predictors: names of the scaled predictor columns
fit_model <- function(model_type, train_data, predictors) {
  f_num <- reformulate(predictors, response = "activity_num")
  f_fac <- reformulate(predictors, response = "activity_factor")
  switch(model_type,
    "lr"   = glm(f_num, data = train_data, family = "binomial"),
    "tree" = rpart(f_fac, data = train_data, method = "class"),
    "rf"   = randomForest(f_fac, data = train_data, ntree = 100),
    "svm"  = svm(f_fac, data = train_data, kernel = "radial", probability = TRUE)
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

# Run LOSO for every held-out subject, training environment, model and test
# environment using the given feature(s) as predictors
# Each feature is min-max scaled with the training fold's min and max
run_loso <- function(features) {
  predictors <- paste0(features, "_scaled")
  pred_log <- list()

  for (test_subject in subjects) {
    env_tab <- env_folds[[test_subject]]

    for (train_env in envs) {
      train <- env_tab[[train_env]]$sessions
      train <- train[train$subject != test_subject, ]

      low  <- sapply(features, function(f) min(train[[f]]))
      high <- sapply(features, function(f) max(train[[f]]))
      for (f in features)
        train[[paste0(f, "_scaled")]] <- (train[[f]] - low[[f]]) / (high[[f]] - low[[f]])
      train$activity_factor <- factor(train$activity, levels = c("standing", "walking"))

      for (model_name in names(model_types)) {
        model_type <- model_types[[model_name]]
        model <- fit_model(model_type, train, predictors)

        for (test_env in envs) {
          test <- env_tab[[test_env]]$sessions
          test <- test[test$subject == test_subject, ]
          for (f in features)
            test[[paste0(f, "_scaled")]] <- (test[[f]] - low[[f]]) / (high[[f]] - low[[f]])
          test$activity_factor <- factor(test$activity, levels = c("standing", "walking"))

          probs <- predict_walking_prob(model, model_type, test)

          pred_log[[length(pred_log) + 1]] <- data.frame(
            test_subject = test_subject,
            model = model_name,
            train_env = train_env,
            test_env = test_env,
            selected_feature = paste(features, collapse = " + "),
            session_id = test$session_id,
            actual = test$activity_factor,
            prob = unname(probs),
            pred = factor(ifelse(probs > 0.5, "walking", "standing"),
                          levels = c("standing", "walking"))
          )
        }
      }
    }
  }

  do.call(rbind, pred_log)
}

# Accuracy, sensitivity, specificity, F1, balanced accuracy and MCC for every
# model / training environment / test environment / held-out subject
# The confusion matrix of every fold is also saved as a PNG
compute_fold_metrics <- function(pred_log, cm_dir) {
  if (!is.null(cm_dir)) dir.create(cm_dir, recursive = TRUE, showWarnings = FALSE)
  fold_metrics_all <- data.frame()

  for (model_name in names(model_types)) {
    for (train_env in envs) {
      for (test_env in envs) {
        for (test_subject in subjects) {
          d <- pred_log[pred_log$model == model_name &
                        pred_log$train_env == train_env &
                        pred_log$test_env == test_env &
                        pred_log$test_subject == test_subject, ]
          cm <- caret::confusionMatrix(d$pred, d$actual, positive = "walking")

          # Matthews correlation coefficient (rows of cm$table = predicted, columns = actual)
          TP <- as.numeric(cm$table["walking", "walking"])
          TN <- as.numeric(cm$table["standing", "standing"])
          FP <- as.numeric(cm$table["walking", "standing"])
          FN <- as.numeric(cm$table["standing", "walking"])
          mcc_denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
          mcc <- ifelse(mcc_denom == 0, 0, (TP * TN - FP * FN) / mcc_denom)

          if (!is.null(cm_dir)) {
            cm_df <- as.data.frame(cm$table)
            p_cm <- ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
              geom_tile(color = "white") +
              geom_text(aes(label = Freq), size = 8, fontface = "bold") +
              scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 10)) +
              scale_y_discrete(limits = rev(levels(cm_df$Reference))) +
              labs(
                title = paste0(model_name, ": ", env_labels[[train_env]], " -> ",
                               env_labels[[test_env]]),
                subtitle = paste("Held-out subject:", test_subject),
                x = "Predicted", y = "Actual"
              ) +
              theme_minimal(base_size = 16) +
              theme(legend.position = "none", panel.grid = element_blank())
            ggsave(
              file.path(cm_dir, paste0("cm_", gsub("[^A-Za-z0-9]+", "", model_name), "_",
                                       train_env, "_", test_env, "_", test_subject, ".png")),
              plot = p_cm, width = 5, height = 4.5, dpi = 300
            )
          }

          fold_metrics_all <- rbind(fold_metrics_all, data.frame(
            model = model_name,
            train_environment = env_labels[[train_env]],
            test_environment = env_labels[[test_env]],
            test_subject = test_subject,
            selected_feature = unique(d$selected_feature),
            accuracy    = as.numeric(cm$overall["Accuracy"]),
            sensitivity = as.numeric(cm$byClass["Sensitivity"]),
            specificity = as.numeric(cm$byClass["Specificity"]),
            f1_score = as.numeric(cm$byClass["F1"]),
            balanced_accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
            mcc = mcc
          ))
        }
      }
    }
  }

  fold_metrics_all
}

# STEP 1: Spearman feature selection in every training fold
# (one fold per held-out subject and training environment = 15 folds)
env_folds <- list()
fold_info <- list()

for (test_subject in subjects) {
  env_folds[[test_subject]] <- lapply(setNames(envs, envs), prepare_env_fold, test_subject = test_subject)
  env_tab <- env_folds[[test_subject]]
  fold_info[[test_subject]] <- list(
    fence = sapply(env_tab, function(e) e$fence),
    removed = lapply(env_tab, function(e) e$removed),
    chosen = list(), rho = list()
  )

  for (train_env in envs) {
    train <- env_tab[[train_env]]$sessions
    train <- train[train$subject != test_subject, ]

    selection <- select_feature(train)
    fold_info[[test_subject]]$chosen[[train_env]] <- selection$feature
    fold_info[[test_subject]]$rho[[train_env]] <- selection$rho
  }
}

# Table that shows percentage of data removed by IQR-based filtering
# per activity per subject per environment
removal_table <- do.call(rbind, lapply(subjects, function(s) {
  do.call(rbind, lapply(envs, function(e) {
    p <- fold_info[[s]]$removed[[e]]
    data.frame(held_out = s,
               environment = e,
               fence = round(unname(fold_info[[s]]$fence[e]), 2),
               subject = names(p),
               pct_removed = as.numeric(p))
  }))
}))

removal_table
write.csv(removal_table, "results/removal_table.csv", row.names = FALSE)

# Table that shows feature selection result per fold
# When there are more than one feature have the highest rho, the first feature
# count from left to right is chosen
feature_selection_table <- do.call(rbind, lapply(subjects, function(s) {
  chosen <- fold_info[[s]]$chosen          # named list: one entry per training environment
  rho    <- fold_info[[s]]$rho
  rho_mat <- t(sapply(names(chosen), function(e) round(rho[[e]], 2)))

  data.frame(test_subject     = s,
             train_env        = names(chosen),
             selected_feature = sapply(chosen, paste, collapse = " / "),
             n_tied           = sapply(chosen, length),
             rho_mat,
             row.names = NULL)
}))

feature_selection_table
write.csv(feature_selection_table, "results/feature_selection_per_fold.csv", row.names = FALSE)

# STEP 2: Optimal predictor = feature with the most votes across training folds
# Each fold casts one vote; if k features tie for the largest |rho|, each gets 1/k
# Ties in total votes are broken by the larger mean |Spearman rho| across all folds
fold_votes <- do.call(rbind, lapply(subjects, function(s) {
  do.call(rbind, lapply(envs, function(e) {
    tied <- fold_info[[s]]$chosen[[e]]
    data.frame(feature = tied, vote = 1 / length(tied))
  }))
}))

# Table that shows the votes and mean absolute rho for each feature
feature_votes <- data.frame(
  feature = potential_features,
  votes = sapply(potential_features, function(f)
    round(sum(fold_votes$vote[fold_votes$feature == f]), 2)),
  n_folds_top = sapply(potential_features, function(f)
    sum(fold_votes$feature == f)),
  mean_abs_rho = sapply(potential_features, function(f)
    round(mean(sapply(subjects, function(s)
      sapply(envs, function(e) abs(fold_info[[s]]$rho[[e]][f])))), 3)),
  row.names = NULL
)

feature_votes <- feature_votes[order(-feature_votes$votes, -feature_votes$mean_abs_rho), ]
optimal_feature <- feature_votes$feature[1]
cat("Optimal predictor used in every fold:", optimal_feature, "\n")
write.csv(feature_votes, "results/feature_selection_votes.csv", row.names = FALSE)

# STEP 3: Train and evaluate every fold with the optimal predictor
pred_log <- run_loso(optimal_feature)
fold_metrics_all <- compute_fold_metrics(pred_log, cm_dir = "results/confusion_matrices")
write.csv(fold_metrics_all, "results/loso_fold_results.csv", row.names = FALSE)

# Mean, SD and 95% CI (t-distribution, df = folds - 1) of each metric across LOSO folds
metric_cols <- c("accuracy", "balanced_accuracy", "sensitivity", "specificity", "mcc")

ci_bound <- function(x, side) {
  x <- x[!is.na(x)]
  mean(x) + side * qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
}

fold_ci_table <- fold_metrics_all %>%
  group_by(model, train_environment, test_environment) %>%
  summarise(
    n_folds = n(),
    across(all_of(metric_cols),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                ci_lower = ~ci_bound(.x, -1),
                ci_upper = ~ci_bound(.x, 1)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

# Keep the CI within the valid range of each metric (MCC: -1 to 1, others: 0 to 1)
for (m in metric_cols) {
  lower_limit <- if (m == "mcc") -1 else 0
  fold_ci_table[[paste0(m, "_ci_lower")]] <- pmax(fold_ci_table[[paste0(m, "_ci_lower")]], lower_limit)
  fold_ci_table[[paste0(m, "_ci_upper")]] <- pmin(fold_ci_table[[paste0(m, "_ci_upper")]], 1)
}

fold_ci_table <- fold_ci_table %>% mutate(across(where(is.numeric), ~round(.x, 4)))
fold_ci_table
write.csv(fold_ci_table, "results/loso_fold_ci.csv", row.names = FALSE)

# Mean LOSO performance per model / training environment / test environment
loso_summary <- fold_metrics_all %>%
  group_by(model, train_environment, test_environment) %>%
  summarise(across(c(accuracy, sensitivity, specificity, f1_score, balanced_accuracy, mcc),
                   ~ round(mean(.x, na.rm = TRUE), 4)),
            .groups = "drop")

print(as.data.frame(loso_summary), row.names = FALSE)
write.csv(loso_summary, "results/loso_summary.csv", row.names = FALSE)


# ABLATION STUDY: individual features and combinations of features
# All 15 non-empty subsets of the four candidate features are evaluated with the
# same LOSO folds, outlier fences, scaling and models as the main evaluation.
# No feature selection is done here; each subset is used as given.
set.seed(42)

feature_subsets <- unlist(lapply(seq_along(potential_features), function(k)
  combn(potential_features, k, simplify = FALSE)), recursive = FALSE)
feature_set_names <- sapply(feature_subsets, paste, collapse = " + ")

ablation_fold_results <- do.call(rbind, lapply(feature_subsets, function(fs) {
  cat("Ablation:", paste(fs, collapse = " + "), "\n")
  data.frame(n_features = length(fs), compute_fold_metrics(run_loso(fs)))
})) %>%
  rename(feature_set = selected_feature) %>%
  relocate(n_features, feature_set)

write.csv(ablation_fold_results, "results/ablation_fold_results.csv", row.names = FALSE)

# Mean performance per feature set / model / training environment / test environment
ablation_summary <- ablation_fold_results %>%
  group_by(n_features, feature_set, model, train_environment, test_environment) %>%
  summarise(across(c(accuracy, sensitivity, specificity, f1_score, balanced_accuracy, mcc),
                   ~ round(mean(.x, na.rm = TRUE), 4)),
            .groups = "drop") %>%
  arrange(n_features, factor(feature_set, levels = feature_set_names), model)

write.csv(ablation_summary, "results/ablation_summary.csv", row.names = FALSE)

# Overview: mean accuracy within the training environment (train = test) and
# across environments (train != test) for each feature set and model
ablation_overview <- ablation_fold_results %>%
  mutate(setting = ifelse(train_environment == test_environment,
                          "within_env", "cross_env")) %>%
  group_by(n_features, feature_set, model, setting) %>%
  summarise(accuracy = round(mean(accuracy, na.rm = TRUE), 4), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = setting, names_prefix = "accuracy_",
                     values_from = accuracy) %>%
  relocate(accuracy_within_env, .before = accuracy_cross_env) %>%
  arrange(n_features, factor(feature_set, levels = feature_set_names), model)

print(as.data.frame(ablation_overview), row.names = FALSE)
write.csv(ablation_overview, "results/ablation_overview.csv", row.names = FALSE)
