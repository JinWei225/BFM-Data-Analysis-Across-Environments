# Analysis Code for Wi-Fi Beamforming Feedback Matrix (BFM) Activity Recognition

## 1. Contents

| File | Description |
| --- | --- |
| `Source Code_v1.0.0.R` | Full analysis pipeline (data cleaning → EDA → hypothesis testing → feature engineering → LOSO classification). |
| `results` | Directory that stores the figures generated and final classification results: mean Accuracy / Sensitivity / Specificity for every (model, train environment, test environment) combination under leave-one-subject-out cross-validation. |
| `data` | Directory that provides link and description for the public dataset used in the study. |

---

## 2. Software requirements

R with the following packages:

```r
install.packages(c("lubridate", "ggplot2", "moments", "boot", "dunn.test",
                   "car", "dplyr", "nlme", "caret", "rpart", "randomForest",
                   "e1071", "InformationValue"))
```

Note: `InformationValue` has been archived on CRAN and may need to be installed
from the archive or from source.

---

## 3. Pipeline overview

### 3.1 Preprocessing (lines 14–62)

1. Missing-value check across all columns.
2. Timestamps converted from UTC to `Asia/Kuala_Lumpur` (MYT, UTC+08:00).
3. Subcarrier columns identified by pattern (`^SCIDX`, `Mag$`, `Phase$`).
4. Four per-packet (row-level) features are derived from the subcarrier vectors:
   - `Mean_Magnitude` — row-wise mean of the subcarrier magnitudes
   - `Std_Magnitude` — row-wise standard deviation of the subcarrier magnitudes
   - `Mean_Phase` — row-wise mean of the subcarrier phases
   - `Phase_Coherence` — circular resultant length,
     `sqrt(mean(cos φ)² + mean(sin φ)²)`, bounded in [0, 1]
5. The raw subcarrier and MAC-address columns are dropped; the data are split by
   environment into `df_open`, `df_foil`, `df_nofoil`.

### 3.2 Outlier removal (lines 82–151)

`remove_outliers_iqr()` applies the standard IQR rule to `Mean_Magnitude` within
each environment and retains rows at or above the lower fence
(`Q1 − 1.5 × IQR`); the function reports the fence, the number of rows removed
and the percentage dropped. Only the lower tail is trimmed, because the artefact
of interest is drop-outs in reported magnitude. Histograms of all four features
are plotted before and after cleaning for each environment.

### 3.3 Descriptive analysis (lines 153–270)

- Per-environment row counts, activity proportions and per-subject counts /
  proportions.
- Grouped bar chart of standing vs. walking counts across environments.
- **Session-level aggregation**: each session is collapsed to a single row by
  averaging the four features (`mean_mag`, `std_mag`, `mean_pha`, `pha_coh`).
  All subsequent statistics and models operate at this session level, so that the
  unit of analysis is a recording session rather than an individual packet.
- `describe_features()` produces the mean, median, standard deviation, skewness
  and kurtosis tables reported per environment.

### 3.4 Assumption checks and hypothesis testing (lines 276–387)

- **Shapiro–Wilk** normality test per feature per environment.
- **Levene's test** for homogeneity of variance across environments.
- Because normality and equal variance are not jointly satisfied,
  **Kruskal–Wallis** is used to test for differences across the three
  environments, followed by **Dunn's post-hoc** pairwise comparisons with
  Bonferroni correction.
- Standing vs. walking is compared within each environment using
  **Wilcoxon rank-sum** tests, accompanied by violin + box plots.

### 3.5 Time-series inspection (lines 390–443)

For each session, a 60-second window centred on the session midpoint is
extracted and all four features are plotted against relative time, faceted by
environment and coloured by activity.

### 3.6 Feature engineering (lines 445–687)

Four session-level dynamic features are computed:

- **Variance of the rate of change**: `var_roc_*`, i.e. `var(Δx / Δt)` where the
  differences are taken over time-ordered packets within a session

These are compared between standing and walking per environment (Shapiro–Wilk,
then Welch's t-test or Wilcoxon rank-sum depending on normality), screened for
multicollinearity with the **variance inflation factor**, and visualised as
Spearman correlation heatmaps against the activity label.

`var_roc_pha_coh` — the variance of the rate of change of phase coherence — is
selected as the single predictor carried into the classification stage, and is
min–max scaled.

### 3.7 Leave-one-subject-out (LOSO) classification (lines 690–1029)

Four classifiers are evaluated, each using `var_roc_pha_coh_scaled` as the sole
predictor and `activity` (standing / walking) as the target:

| Model | Implementation |
| --- | --- |
| Logistic Regression | `glm(..., family = "binomial")`, threshold 0.5 |
| Decision Tree | `rpart(..., method = "class")` |
| Random Forest | `randomForest(..., ntree = 100)` |
| SVM | `e1071::svm(..., kernel = "radial", probability = TRUE)` |

Two evaluation protocols are run:

1. **Within-environment LOSO** (lines 690–897) — for each environment, train on
   four subjects and test on the held-out subject; metrics are averaged over the
   five folds and plotted as grouped bar charts of mean accuracy, sensitivity and
   specificity.
2. **Cross-environment LOSO** (lines 899–1030) — for every ordered pair of
   training and testing environments, **including the matched pairs where the
   two are the same**, train on four subjects' sessions from the training
   environment and test on the held-out subject's sessions from the testing
   environment. Min–max scaling is fitted on the training fold only and applied
   to the test fold with the training fold's min/max, so no test information
   leaks into the scaling step. Folds are skipped when a fold would contain only
   one class. This loop alone produces all 36 rows of `loso_results.csv`.

---

## 4. Results file — `loso_results.csv`

One row per (model, train environment, test environment) combination, 36 rows in
total: 4 models × 3 training environments × 3 testing environments.

| Column | Meaning |
| --- | --- |
| `model` | Logistic Regression, Decision Tree, Random Forest, or SVM |
| `train_environment` | Environment supplying the training sessions (`Open`, `Foil`, `No Foil`) |
| `test_environment` | Environment supplying the held-out subject's test sessions |
| `accuracy` | Mean accuracy across the LOSO folds |
| `sensitivity` | Mean sensitivity (walking correctly identified) |
| `specificity` | Mean specificity (standing correctly identified) |
| `n_folds` | Number of folds contributing to the averages (5 = all subjects used) |

Rows where `train_environment == test_environment` are the within-environment
LOSO results; the remaining rows are the cross-environment transfer results.
Both are produced by the same loop, so every row in the file —
diagonal included — uses the leakage-free protocol in which min–max scaling is
fitted on the training fold alone and the held-out subject never appears in
training. The file is the contents of `cross_env_loso_results` at the end of the
script.
