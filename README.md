# Analysis Code for Wi-Fi Beamforming Feedback Matrix (BFM) Activity Recognition

## 1. Contents

| File | Description |
| --- | --- |
| `Analysis_v1.0.0.R` | Statistical analysis (preprocessing → outlier removal → descriptive analysis → hypothesis testing → feature engineering). Saves the packet-level features to `data/bfm_packet_features.rds`. |
| `Modelling_v1.0.0.R` | Data preprocessing, normalization, model training and evaluation using Leave-one-subject-out (LOSO) approach. Reads `data/bfm_packet_features.rds`, so `Analysis_v1.0.0.R` must be run first. |
| `results` | Directory that stores the figures generated, statistical analysis and final classification results: mean Accuracy / Sensitivity / Specificity / F1-Score / Balanced Accuracy for every (model, train environment, test environment) combination under leave-one-subject-out cross-validation. |
| `data` | Directory that provides link and description for the public dataset used in the study and contains packet-level features from `Analysis_v1.0.0.R`. |

---

## 2. Software requirements

R with the following packages:

```r
# Analysis_v1.0.0.R
install.packages(c("lubridate", "ggplot2", "moments", "dunn.test", "car", "dplyr"))

# Modelling_v1.0.0.R
install.packages(c("dplyr", "ggplot2", "caret", "rpart", "randomForest", "e1071"))
```

---

## 3. Pipeline overview

### 3.1 Preprocessing (`Analysis_v1.0.0.R`)

1. Missing-value check across all columns.
2. Timestamps converted from UTC to `Asia/Kuala_Lumpur` (MYT, UTC+08:00).
3. Subcarrier columns identified by pattern (`^SCIDX`, `Mag$`, `Phase$`).
4. Four per-packet (row-level) features are derived from the subcarrier vectors:
   - `Mean_Magnitude` — row-wise mean of the subcarrier magnitudes
   - `Std_Magnitude` — row-wise standard deviation of the subcarrier magnitudes
   - `Mean_Phase` — row-wise mean of the subcarrier phases
   - `Phase_Coherence` — circular resultant length,
     `sqrt(mean(cos φ)² + mean(sin φ)²)`, bounded in [0, 1]
5. The raw subcarrier and MAC-address columns are dropped and the packet-level
   features are saved to `data/bfm_packet_features.rds` for `Modelling_v1.0.0.R`;
   the data are split by environment into `df_open`, `df_foil`, `df_nofoil`.

### 3.2 Outlier removal (`Analysis_v1.0.0.R`)

`remove_outliers_iqr()` applies the standard IQR rule to `Mean_Magnitude` within
each environment and retains rows at or above the lower fence
(`Q1 − 1.5 × IQR`). The rows dropped also helps to remove noise in the other three
features. The code also reports the lower fence, the number of rows removed
and the percentage dropped. Histograms of all four features
are plotted before and after cleaning for each environment, which shows why
only the lower tail is removed.

### 3.3 Descriptive analysis (`Analysis_v1.0.0.R`)

- Per-environment row counts, activity proportions and per-subject counts /
  proportions in percentage.
- Grouped bar chart of standing vs. walking counts across environments.
- **Session-level aggregation**: each session is collapsed to a single row by
  averaging the four features (`mean_mag`, `std_mag`, `mean_pha`, `pha_coh`).
  All subsequent statistics and models operate at this session level, so that the
  unit of analysis is a recording session rather than an individual packet.
- `describe_features()` produces the mean, median, standard deviation, skewness
  and kurtosis tables reported per environment.

### 3.4 Assumption checks and hypothesis testing (`Analysis_v1.0.0.R`)

- **Shapiro–Wilk** normality test per feature per environment.
- **Levene's test** for homogeneity of variance across environments.
- Because normality and equal variance are not jointly satisfied,
  **Kruskal–Wallis** is used to test for differences across the three
  environments. Sessions from the same subject are not independent, so the
  test is run **within each subject** (10 sessions per environment), with Holm
  correction across subjects and epsilon-squared as the effect size. It is
  followed by **Dunn's post-hoc** pairwise comparisons per subject with
  Bonferroni correction (`kruskal_wallis_per_subject.csv`, `dunn_per_subject.csv`).
- Standing vs. walking is compared within each environment using
  **Wilcoxon rank-sum** tests, accompanied by violin + box plots.

### 3.5 Time-series inspection (`Analysis_v1.0.0.R`)

For each session, a 60-second window centred on the session midpoint is
extracted and all four features are plotted against relative time, faceted by
environment and coloured by activity.

### 3.6 Feature engineering (`Analysis_v1.0.0.R`)

Four session-level dynamic features are computed:

- **Variance of the rate of change**: `var_roc_*`, i.e. `var(Δx / Δt)` where the
  differences are taken over time-ordered packets within a session

These are compared between standing and walking per environment (Shapiro–Wilk,
then Welch's t-test or Wilcoxon rank-sum depending on normality), screened for
multicollinearity with the **variance inflation factor**, and visualised as
Spearman correlation heatmaps against the activity label.

### 3.7 Leave-one-subject-out (LOSO) classification (`Modelling_v1.0.0.R`)

Previous multicollinearity checking shows that multicollinearity issue exists,
so only one single predictor will be chosen.The single predictor is chosen by 
voting from all folds using Spearman correlation in two steps:

1. In each of the 15 training folds (5 held-out subjects × 3 training
   environments), the feature with the largest |Spearman ρ| with activity is
   selected (`feature_selection_per_fold.csv`). When k features tie, each
   receives 1/k of that fold's vote.
2. The feature with the most votes across all folds is used as the predictor in
   every fold (`feature_selection_votes.csv`). `var_roc_mean_mag` is actually 
   chosen after voting.

Four classifiers are evaluated without hyperparameter tuning, each using the 
min–max scaled `var_roc_mean_mag` as the sole predictor and `activity` 
(standing / walking) as the target:

| Model | Implementation |
| --- | --- |
| Logistic Regression | `glm(..., family = "binomial")`, threshold 0.5 |
| Decision Tree | `rpart(..., method = "class")` |
| Random Forest | `randomForest(..., ntree = 100)` |
| SVM | `e1071::svm(..., kernel = "radial", probability = TRUE)` |

A single loop covers every ordered pair of training and testing environments,
**including the matched pairs where the two are the same**. In each fold the
model is trained on four subjects' sessions from the training environment and
tested on the held-out subject's sessions from the testing environment. No
test information leaks into training:

- The outlier lower fence for each environment is computed from the training
  subjects' packets only and then applied to all subjects in that environment.
  The fence and the percentage of packets removed per subject in every fold are
  saved to `removal_table.csv`.
- Min–max scaling is fitted on the training fold and applied to the test fold
  with the training fold's min/max.

### 3.8 Ablation study (`Modelling_v1.0.0.R`)

Every non-empty subset of the four candidate features (4 single features,
6 pairs, 4 triples and all four together = 15 feature sets) is evaluated with
the same LOSO folds, outlier fences, scaling and models as in 3.7. No feature
selection is done; each subset is used as given. The ablation results are
provided in this repository only and are not reported in the paper.

---

## 4. Results files

`loso_fold_results.csv` has one row per (model, train environment, test
environment, held-out subject): 4 models × 3 training environments × 3 testing
environments × 5 subjects = 180 rows. `loso_summary.csv` averages these over the
five folds, giving 36 rows.

| Column | Meaning |
| --- | --- |
| `model` | Logistic Regression, Decision Tree, Random Forest, or SVM |
| `train_environment` | Environment supplying the training sessions (`Open`, `Foil`, `No Foil`) |
| `test_environment` | Environment supplying the held-out subject's test sessions |
| `test_subject` | Held-out subject (`loso_fold_results.csv` only) |
| `selected_feature` | Predictor used in the fold (`loso_fold_results.csv` only) |
| `accuracy` | Accuracy |
| `sensitivity` | Sensitivity (walking correctly identified) |
| `specificity` | Specificity (standing correctly identified) |
| `f1_score` | F1 score for the walking class (`NA` when no session is predicted as walking) |
| `balanced_accuracy` | Mean of sensitivity and specificity |
[ `mcc` | Matthews correlation coefficient |

Rows where `train_environment == test_environment` are the within-environment
LOSO results; the remaining rows are the cross-environment transfer results.

Other files written by `Modelling_v1.0.0.R`:

| File | Contents |
| --- | --- |
| `removal_table.csv` | Lower fence and percentage of packets removed per subject, for every held-out subject and environment |
| `feature_selection_per_fold.csv` | Spearman ρ of each feature and the selected (or tied) feature(s) in each training fold |
| `feature_selection_votes.csv` | Split votes per feature across the 15 training folds |
| `ablation_fold_results.csv` | Same as `loso_fold_results.csv` for every feature set, with `n_features` and `feature_set` columns |
| `ablation_summary.csv` | Same as `loso_summary.csv` for every feature set |
| `ablation_overview.csv` | Mean accuracy within the training environment (`accuracy_within_env`) and across environments (`accuracy_cross_env`) per feature set and model |
