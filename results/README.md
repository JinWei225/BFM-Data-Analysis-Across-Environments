# Results

This folder holds every table, figure and confusion matrix produced by
`Analysis_v1.0.0.R` and `Modelling_v1.0.0.R`. It also holds the supplementary
tables that the paper refers to GitHub for: the model hyperparameters, IQR
filtering per LOSO fold, Dunn's post-hoc tests, the uncertainty of the LOSO
results and the ablation study.

All numbers below are taken directly from the CSV files in this folder. Values are
rounded to two decimal places (three for z-statistics) for readability; the CSV
files keep the full precision.

## Contents

1. [Files in this folder](#1-files-in-this-folder)
2. [Model hyperparameters](#2-model-hyperparameters)
3. [IQR-based filtering per LOSO fold](#3-iqr-based-filtering-per-loso-fold)
4. [Dunn's post-hoc comparisons per subject](#4-dunns-post-hoc-comparisons-per-subject)
5. [LOSO results with standard deviation and 95% confidence interval](#5-loso-results-with-standard-deviation-and-95-confidence-interval)
6. [Ablation study](#6-ablation-study)
7. [Confusion matrices](#7-confusion-matrices)

---

## 1. Files in this folder

| File / folder | Written by | Contents |
| --- | --- | --- |
| `sessions_per_subject.csv` | `Analysis_v1.0.0.R` | Number of sessions per environment, subject and activity (5 standing + 5 walking per subject per environment) |
| `kruskal_wallis_per_subject.csv` | `Analysis_v1.0.0.R` | Kruskal–Wallis test of each feature across the three environments within each subject: H, p, Holm-adjusted p, ε² (Table II of the paper) |
| `dunn_per_subject.csv` | `Analysis_v1.0.0.R` | Dunn's post-hoc pairwise comparisons per subject (Section 4 below) |
| `figures/` | `Analysis_v1.0.0.R` | Histograms before/after IQR cleaning, activity bar chart, violin plots, time-series plots and Spearman heatmaps |
| `removal_table.csv` | `Modelling_v1.0.0.R` | IQR lower fence and % of packets removed per subject in every LOSO fold (Section 3 below) |
| `feature_selection_per_fold.csv` | `Modelling_v1.0.0.R` | Spearman ρ of each feature with activity and the selected (or tied) feature(s) in each of the 15 training folds |
| `feature_selection_votes.csv` | `Modelling_v1.0.0.R` | Split votes, number of folds ranked top and mean \|ρ\| per feature (Table III of the paper) |
| `loso_fold_results.csv` | `Modelling_v1.0.0.R` | Accuracy, sensitivity, specificity, F1, balanced accuracy and MCC for each of the 180 folds (4 models × 3 training × 3 testing environments × 5 held-out subjects) |
| `loso_summary.csv` | `Modelling_v1.0.0.R` | Mean of each metric over the 5 folds, 36 rows (Table IV of the paper) |
| `loso_fold_ci.csv` | `Modelling_v1.0.0.R` | Mean, SD and 95% CI of each metric over the 5 folds (Section 5 below) |
| `ablation_fold_results.csv` | `Modelling_v1.0.0.R` | Same as `loso_fold_results.csv` for each of the 15 feature sets |
| `ablation_summary.csv` | `Modelling_v1.0.0.R` | Same as `loso_summary.csv` for each of the 15 feature sets |
| `ablation_overview.csv` | `Modelling_v1.0.0.R` | Mean accuracy within and across environments per feature set and model (Section 6 below) |
| `confusion_matrices/` | `Modelling_v1.0.0.R` | 180 confusion-matrix plots, one per fold (Section 7 below) |

---

## 2. Model hyperparameters

No hyperparameter tuning (grid search, random search or nested cross-validation)
was done. Every model is fitted once per LOSO fold on the four training subjects'
sessions. Hyperparameters are left at the R package defaults except for the two
that are marked with **(changed)**. `set.seed(42)` is called once before the main
LOSO loop and once before the ablation study, so the random forest and the
SVM probability model can be reproduced.

The predictor is min–max scaled to [0, 1] using the training fold's minimum and
maximum, and the same values are then applied to the held-out subject. A session
is classified as walking when the predicted probability of walking is greater than 0.5.

| Model | R function (package version used) | Hyperparameter | Value |
| --- | --- | --- | --- |
| Logistic Regression | `glm(activity_num ~ x, family = "binomial")` (stats, R 4.3) | Link function | logit |
| | | Regularisation | none |
| | | Fitting method / max. iterations / convergence tolerance | IRLS (`glm.fit`) / 25 / 1e-8 |
| | | Decision threshold | 0.5 |
| Decision Tree | `rpart(activity_factor ~ x, method = "class")` (rpart 4.1.23) | Split criterion | Gini index |
| | | `minsplit` (min. sessions in a node to try a split) | 20 |
| | | `minbucket` (min. sessions in a leaf) | 7 (= round(20 / 3)) |
| | | `cp` (complexity parameter) | 0.01 |
| | | `maxdepth` | 30 |
| | | `xval` (internal CV folds, only used for the cp table) | 10 — no pruning is applied |
| | | Class priors | Training-set class proportions (balanced, 0.5 / 0.5) |
| Random Forest | `randomForest(activity_factor ~ x, ntree = 100)` (randomForest 4.7-1.1) | `ntree` | **100 (changed; default 500)** — a single predictor does not need 500 trees for a binary task |
| | | `mtry` | floor(√p) = 1 for the single predictor (1–2 in the ablation study) |
| | | `nodesize` (min. leaf size) | 1 |
| | | `maxnodes` | no limit |
| | | Bootstrap sampling | with replacement, sample size = number of training sessions |
| SVM | `svm(activity_factor ~ x, kernel = "radial", probability = TRUE)` (e1071 1.7-14) | Type | C-classification |
| | | Kernel | radial basis function (RBF) |
| | | `cost` (C) | 1 |
| | | `gamma` | 1 / p = 1 for the single predictor (1/2 to 1/4 in the ablation study) |
| | | `scale` | TRUE — e1071 also standardises the inputs internally using the training data |
| | | `tolerance` / `shrinking` / `cachesize` | 0.001 / TRUE / 40 MB |
| | | `class.weights` | none |
| | | `probability` | **TRUE (changed; default FALSE)** — needed so that the SVM outputs a probability that can be compared with the 0.5 threshold (Platt scaling, fitted on the training fold) |

The package versions in the table were used to regenerate the results on
13 Oct 2026. They reproduce the committed CSV files exactly.

---

## 3. IQR-based filtering per LOSO fold

Noise appears as a cluster of packets with a very low `Mean_Magnitude`. Only the
lower tail is removed: packets with `Mean_Magnitude < Q1 − 1.5 × IQR` are dropped.
The filtering is done at two separate stages, both starting from the original
uncleaned data:

**(a) Analysis stage (`Analysis_v1.0.0.R`).** For descriptive statistics and hypothesis
testing only, one fence is computed per environment using all subjects:

| Environment | Lower fence | Packets removed | % removed |
| --- | ---: | ---: | ---: |
| Open | 10.72 | 1,172 of 47,236 | 2.48 |
| Foil | 4.65 | 3,306 of 48,438 | 6.83 |
| No Foil | 7.53 | 1,415 of 47,540 | 2.98 |

**(b) Modelling stage (`Modelling_v1.0.0.R`).** To avoid information leakage, in
every LOSO fold the fence for each environment is computed from the **four training
subjects only**. That fence is then applied to all five subjects, including the
held-out one. The tables below show the fence and the percentage of each subject's
packets removed. Each row is one LOSO fold, and the bold cell is the held-out
(test) subject.

**Open**

| Held-out subject | Lower fence | abel | collin | ivan | kenny | matthew |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| abel | 10.62 | **2.27** | 2.34 | 2.94 | 2.45 | 1.60 |
| collin | 10.66 | 2.27 | **2.34** | 3.05 | 2.47 | 1.62 |
| ivan | 10.79 | 2.28 | 2.34 | **5.59** | 2.50 | 3.03 |
| kenny | 10.66 | 2.27 | 2.34 | 3.09 | **2.47** | 1.62 |
| matthew | 10.80 | 2.28 | 2.34 | 6.00 | 2.51 | **3.40** |

**Foil**

| Held-out subject | Lower fence | abel | collin | ivan | kenny | matthew |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| abel | 4.69 | **5.99** | 3.61 | 3.31 | 3.20 | 16.32 |
| collin | 6.37 | 5.99 | **3.61** | 3.31 | 3.20 | 16.32 |
| ivan | 5.00 | 5.99 | 3.61 | **3.31** | 3.20 | 16.32 |
| kenny | 3.58 | 5.94 | 3.61 | 3.31 | **3.20** | 16.32 |
| matthew | 3.97 | 5.98 | 3.61 | 3.31 | 3.20 | **16.32** |

**No Foil**

| Held-out subject | Lower fence | abel | collin | ivan | kenny | matthew |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| abel | 6.75 | **1.37** | 1.41 | 3.21 | 4.72 | 4.07 |
| collin | 7.21 | 1.37 | **1.41** | 3.21 | 4.72 | 4.07 |
| ivan | 9.10 | 1.37 | 1.41 | **3.21** | 4.72 | 4.07 |
| kenny | 7.55 | 1.37 | 1.41 | 3.21 | **4.72** | 4.07 |
| matthew | 6.33 | 1.37 | 1.41 | 3.21 | 4.72 | **4.07** |

---

## 4. Dunn's post-hoc comparisons per subject

The Kruskal–Wallis test (Table II of the paper) is run within each subject, on the
30 session-level observations of that subject (10 per environment). Dunn's test then
compares each pair of environments. P-values are two-sided and Bonferroni-adjusted
over the three comparisons within each subject and feature.

The z-statistic for "A − B" is based on the mean rank of A minus the mean rank of B.
**A positive z means the feature tends to be larger in environment A.** Bold cells
have an adjusted p < 0.05. Cells show z (adjusted p).

| Feature | Subject | Foil − No Foil | Foil − Open | No Foil − Open |
| --- | --- | ---: | ---: | ---: |
| Mean Magnitude | abel | **4.724 (< 0.001)** | **2.896 (0.011)** | -1.829 (0.202) |
|  | collin | **2.667 (0.023)** | -0.686 (1.000) | **-3.353 (0.002)** |
|  | ivan | **-2.642 (0.025)** | 1.727 (0.252) | **4.369 (< 0.001)** |
|  | kenny | **2.413 (0.048)** | **4.978 (< 0.001)** | **2.565 (0.031)** |
|  | matthew | **2.540 (0.033)** | **5.080 (< 0.001)** | **2.540 (0.033)** |
| Mean Phase | abel | -0.152 (1.000) | **3.505 (0.001)** | **3.658 (< 0.001)** |
|  | collin | -2.159 (0.092) | **-3.327 (0.003)** | -1.168 (0.728) |
|  | ivan | **4.166 (< 0.001)** | 1.016 (0.929) | **-3.150 (0.005)** |
|  | kenny | **3.302 (0.003)** | -0.940 (1.000) | **-4.242 (< 0.001)** |
|  | matthew | -1.727 (0.252) | 0.432 (1.000) | 2.159 (0.092) |
| Std Magnitude | abel | **4.775 (< 0.001)** | **2.845 (0.013)** | -1.930 (0.161) |
|  | collin | **3.480 (0.002)** | 1.245 (0.640) | -2.235 (0.076) |
|  | ivan | -2.032 (0.127) | 2.261 (0.071) | **4.293 (< 0.001)** |
|  | kenny | **2.540 (0.033)** | **5.080 (< 0.001)** | **2.540 (0.033)** |
|  | matthew | **2.540 (0.033)** | **5.080 (< 0.001)** | **2.540 (0.033)** |
| Phase Coherence | abel | **-4.928 (< 0.001)** | **-2.692 (0.021)** | 2.235 (0.076) |
|  | collin | 0.508 (1.000) | -1.346 (0.535) | -1.854 (0.191) |
|  | ivan | **2.616 (0.027)** | -2.388 (0.051) | **-5.004 (< 0.001)** |
|  | kenny | **2.540 (0.033)** | **-2.540 (0.033)** | **-5.080 (< 0.001)** |
|  | matthew | **-2.540 (0.033)** | **-5.080 (< 0.001)** | **-2.540 (0.033)** |

---

## 5. LOSO results with standard deviation and 95% confidence interval

Table IV of the paper reports only the mean over the five LOSO folds. Here each
cell is **mean ± SD [95% CI]** over the five folds, where each fold is one held-out
subject. The CI uses the t-distribution with 4 degrees of freedom and is clipped
to the valid range of the metric. With only five folds the intervals are wide and
should be read as a rough indication of uncertainty.

* Balanced accuracy is identical to accuracy in every fold, because each held-out
  subject has 5 standing and 5 walking sessions per environment, so it is not repeated
  here. It is in `loso_fold_ci.csv`.
* F1-score is undefined in the 50 folds where the model predicts no session as
  walking (sensitivity = 0), so it is not summarised. The per-fold F1 values are in
  `loso_fold_results.csv`.
* MCC is set to 0 in those same folds, where it is otherwise undefined.

**Logistic Regression**

| Train → Test | Accuracy | Sensitivity | Specificity | MCC |
| --- | --- | --- | --- | --- |
| Open → Open | 0.98 ± 0.04 [0.92, 1.00] | 0.96 ± 0.09 [0.85, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.96 ± 0.08 [0.86, 1.00] |
| Open → Foil | 0.52 ± 0.04 [0.46, 0.58] | 1.00 ± 0.00 [1.00, 1.00] | 0.04 ± 0.09 [0.00, 0.15] | 0.07 ± 0.15 [-0.12, 0.25] |
| Open → No Foil | 0.78 ± 0.23 [0.50, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.56 ± 0.46 [0.00, 1.00] | 0.60 ± 0.43 [0.06, 1.00] |
| Foil → Open | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| Foil → Foil | 0.80 ± 0.19 [0.57, 1.00] | 0.76 ± 0.33 [0.35, 1.00] | 0.84 ± 0.36 [0.40, 1.00] | 0.66 ± 0.31 [0.28, 1.00] |
| Foil → No Foil | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| No Foil → Open | 0.60 ± 0.14 [0.42, 0.78] | 0.20 ± 0.28 [0.00, 0.55] | 1.00 ± 0.00 [1.00, 1.00] | 0.23 ± 0.32 [-0.17, 0.63] |
| No Foil → Foil | 0.60 ± 0.14 [0.42, 0.78] | 0.96 ± 0.09 [0.85, 1.00] | 0.24 ± 0.36 [0.00, 0.68] | 0.22 ± 0.30 [-0.16, 0.60] |
| No Foil → No Foil | 0.74 ± 0.23 [0.45, 1.00] | 0.68 ± 0.41 [0.17, 1.00] | 0.80 ± 0.45 [0.24, 1.00] | 0.49 ± 0.47 [-0.09, 1.00] |

**Decision Tree**

| Train → Test | Accuracy | Sensitivity | Specificity | MCC |
| --- | --- | --- | --- | --- |
| Open → Open | 0.98 ± 0.04 [0.92, 1.00] | 0.96 ± 0.09 [0.85, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.96 ± 0.08 [0.86, 1.00] |
| Open → Foil | 0.52 ± 0.04 [0.46, 0.58] | 1.00 ± 0.00 [1.00, 1.00] | 0.04 ± 0.09 [0.00, 0.15] | 0.07 ± 0.15 [-0.12, 0.25] |
| Open → No Foil | 0.78 ± 0.23 [0.50, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.56 ± 0.46 [0.00, 1.00] | 0.60 ± 0.43 [0.06, 1.00] |
| Foil → Open | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| Foil → Foil | 0.76 ± 0.19 [0.52, 1.00] | 0.76 ± 0.33 [0.35, 1.00] | 0.76 ± 0.43 [0.22, 1.00] | 0.56 ± 0.37 [0.09, 1.00] |
| Foil → No Foil | 0.54 ± 0.09 [0.43, 0.65] | 0.08 ± 0.18 [0.00, 0.30] | 1.00 ± 0.00 [1.00, 1.00] | 0.10 ± 0.22 [-0.18, 0.38] |
| No Foil → Open | 0.80 ± 0.27 [0.46, 1.00] | 0.60 ± 0.55 [0.00, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.60 ± 0.55 [-0.08, 1.00] |
| No Foil → Foil | 0.60 ± 0.17 [0.38, 0.82] | 0.96 ± 0.09 [0.85, 1.00] | 0.24 ± 0.43 [0.00, 0.78] | 0.23 ± 0.36 [-0.21, 0.67] |
| No Foil → No Foil | 0.60 ± 0.10 [0.48, 0.72] | 0.68 ± 0.46 [0.11, 1.00] | 0.52 ± 0.46 [0.00, 1.00] | 0.27 ± 0.25 [-0.05, 0.58] |

**Random Forest**

| Train → Test | Accuracy | Sensitivity | Specificity | MCC |
| --- | --- | --- | --- | --- |
| Open → Open | 0.98 ± 0.04 [0.92, 1.00] | 0.96 ± 0.09 [0.85, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.96 ± 0.08 [0.86, 1.00] |
| Open → Foil | 0.52 ± 0.04 [0.46, 0.58] | 1.00 ± 0.00 [1.00, 1.00] | 0.04 ± 0.09 [0.00, 0.15] | 0.07 ± 0.15 [-0.12, 0.25] |
| Open → No Foil | 0.78 ± 0.23 [0.50, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.56 ± 0.46 [0.00, 1.00] | 0.60 ± 0.43 [0.06, 1.00] |
| Foil → Open | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| Foil → Foil | 0.78 ± 0.13 [0.62, 0.94] | 0.72 ± 0.33 [0.30, 1.00] | 0.84 ± 0.26 [0.52, 1.00] | 0.62 ± 0.21 [0.36, 0.88] |
| Foil → No Foil | 0.52 ± 0.15 [0.34, 0.70] | 0.12 ± 0.18 [0.00, 0.34] | 0.92 ± 0.18 [0.70, 1.00] | 0.07 ± 0.38 [-0.41, 0.54] |
| No Foil → Open | 0.66 ± 0.15 [0.47, 0.85] | 0.32 ± 0.30 [0.00, 0.70] | 1.00 ± 0.00 [1.00, 1.00] | 0.40 ± 0.30 [0.03, 0.77] |
| No Foil → Foil | 0.74 ± 0.23 [0.45, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.48 ± 0.46 [0.00, 1.00] | 0.49 ± 0.47 [-0.09, 1.00] |
| No Foil → No Foil | 0.70 ± 0.12 [0.55, 0.85] | 0.64 ± 0.26 [0.32, 0.96] | 0.76 ± 0.43 [0.22, 1.00] | 0.45 ± 0.26 [0.13, 0.77] |

**SVM**

| Train → Test | Accuracy | Sensitivity | Specificity | MCC |
| --- | --- | --- | --- | --- |
| Open → Open | 0.98 ± 0.04 [0.92, 1.00] | 0.96 ± 0.09 [0.85, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.96 ± 0.08 [0.86, 1.00] |
| Open → Foil | 0.54 ± 0.05 [0.47, 0.61] | 1.00 ± 0.00 [1.00, 1.00] | 0.08 ± 0.11 [0.00, 0.22] | 0.13 ± 0.18 [-0.09, 0.36] |
| Open → No Foil | 0.78 ± 0.23 [0.50, 1.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.56 ± 0.46 [0.00, 1.00] | 0.60 ± 0.43 [0.06, 1.00] |
| Foil → Open | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| Foil → Foil | 0.80 ± 0.19 [0.57, 1.00] | 0.76 ± 0.33 [0.35, 1.00] | 0.84 ± 0.36 [0.40, 1.00] | 0.66 ± 0.31 [0.28, 1.00] |
| Foil → No Foil | 0.50 ± 0.00 [0.50, 0.50] | 0.00 ± 0.00 [0.00, 0.00] | 1.00 ± 0.00 [1.00, 1.00] | 0.00 ± 0.00 [0.00, 0.00] |
| No Foil → Open | 0.54 ± 0.09 [0.43, 0.65] | 0.08 ± 0.18 [0.00, 0.30] | 1.00 ± 0.00 [1.00, 1.00] | 0.10 ± 0.22 [-0.18, 0.38] |
| No Foil → Foil | 0.64 ± 0.15 [0.45, 0.83] | 0.96 ± 0.09 [0.85, 1.00] | 0.32 ± 0.36 [0.00, 0.77] | 0.32 ± 0.31 [-0.07, 0.71] |
| No Foil → No Foil | 0.72 ± 0.22 [0.45, 0.99] | 0.64 ± 0.41 [0.13, 1.00] | 0.80 ± 0.45 [0.24, 1.00] | 0.46 ± 0.44 [-0.09, 1.00] |

---

## 6. Ablation study

To check whether the performance comes from the selected feature, every non-empty
subset of the four engineered features was evaluated: 4 single features, 6 pairs,
4 triples and all 4 features (15 feature sets in total). Each set used the same LOSO
folds, IQR fences, min–max scaling and models as the main evaluation, with no feature
selection inside the ablation.

* **within** = mean accuracy over the 3 matched pairs (Open→Open, Foil→Foil,
  No Foil→No Foil) × 5 folds
* **cross** = mean accuracy over the 6 cross-environment pairs × 5 folds

The `var_roc_mean_mag` row is the predictor chosen by the Spearman vote and used in
the paper. Its results are identical to `loso_summary.csv`. The full metrics per
training/testing environment are in `ablation_summary.csv`, and per fold in
`ablation_fold_results.csv`.

| # | Feature set | LR within | DT within | RF within | SVM within | LR cross | DT cross | RF cross | SVM cross |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `var_roc_mean_mag` *(selected)* | 0.84 | 0.78 | 0.82 | 0.83 | 0.58 | 0.62 | 0.62 | 0.58 |
| 1 | `var_roc_mean_pha` | 0.77 | 0.74 | 0.75 | 0.76 | 0.60 | 0.63 | 0.65 | 0.60 |
| 1 | `var_roc_std_mag` | 0.85 | 0.87 | 0.85 | 0.87 | 0.57 | 0.64 | 0.66 | 0.56 |
| 1 | `var_roc_pha_coh` | 0.79 | 0.79 | 0.83 | 0.79 | 0.64 | 0.65 | 0.68 | 0.64 |
| 2 | `var_roc_mean_mag` + `var_roc_mean_pha` | 0.79 | 0.75 | 0.79 | 0.84 | 0.59 | 0.62 | 0.64 | 0.58 |
| 2 | `var_roc_mean_mag` + `var_roc_std_mag` | 0.89 | 0.85 | 0.90 | 0.86 | 0.63 | 0.63 | 0.63 | 0.59 |
| 2 | `var_roc_mean_mag` + `var_roc_pha_coh` | 0.79 | 0.72 | 0.78 | 0.84 | 0.59 | 0.62 | 0.64 | 0.58 |
| 2 | `var_roc_mean_pha` + `var_roc_std_mag` | 0.86 | 0.84 | 0.81 | 0.87 | 0.58 | 0.64 | 0.66 | 0.57 |
| 2 | `var_roc_mean_pha` + `var_roc_pha_coh` | 0.75 | 0.76 | 0.75 | 0.78 | 0.63 | 0.66 | 0.65 | 0.62 |
| 2 | `var_roc_std_mag` + `var_roc_pha_coh` | 0.83 | 0.82 | 0.80 | 0.88 | 0.64 | 0.60 | 0.67 | 0.58 |
| 3 | `var_roc_mean_mag` + `var_roc_mean_pha` + `var_roc_std_mag` | 0.84 | 0.85 | 0.86 | 0.89 | 0.61 | 0.63 | 0.62 | 0.58 |
| 3 | `var_roc_mean_mag` + `var_roc_mean_pha` + `var_roc_pha_coh` | 0.79 | 0.72 | 0.79 | 0.86 | 0.60 | 0.63 | 0.63 | 0.58 |
| 3 | `var_roc_mean_mag` + `var_roc_std_mag` + `var_roc_pha_coh` | 0.82 | 0.83 | 0.85 | 0.89 | 0.63 | 0.61 | 0.64 | 0.59 |
| 3 | `var_roc_mean_pha` + `var_roc_std_mag` + `var_roc_pha_coh` | 0.81 | 0.83 | 0.80 | 0.87 | 0.60 | 0.64 | 0.65 | 0.59 |
| 4 | `var_roc_mean_mag` + `var_roc_mean_pha` + `var_roc_std_mag` + `var_roc_pha_coh` | 0.81 | 0.83 | 0.84 | 0.89 | 0.61 | 0.61 | 0.64 | 0.58 |

LR = Logistic Regression, DT = Decision Tree, RF = Random Forest, SVM = Support Vector Machine.

---

## 7. Confusion matrices

`confusion_matrices/` has one plot for each of the 180 LOSO folds, named

```
cm_<Model>_<trainEnv>_<testEnv>_<heldOutSubject>.png
e.g. cm_RandomForest_nofoil_open_abel.png
```

where `<Model>` is `LogisticRegression`, `DecisionTree`, `RandomForest` or `SVM`
and the environments are `open`, `foil` or `nofoil`. Each matrix covers the held-out
subject's 10 sessions (5 standing, 5 walking) in the test environment.
