############################################################
# Peritoneal Metastasis Prediction in Gastric Cancer
# IMPORTANT:
# - The paper did not publish the exact scripts.
############################################################

############################
# 0) PACKAGES
############################
required_packages <- c(
  "readxl", "dplyr", "stringr", "forcats", "ggplot2", "tidyr",
  "mice", "tidymodels", "rpart", "ranger",
  "kernlab", "kknn", "klaR", "vip", "pROC"
)

to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(tidyr)
library(mice)
library(tidymodels)
library(rpart)
library(ranger)
library(kernlab)
library(kknn)
library(klaR)
library(vip)
library(pROC)

set.seed(1234)
options(yardstick.event_first = FALSE)

############################################################
# 01_preprocessing.R
############################################################

############################
# 1.1) FILE PATH
############################
file_path <- "MLG7GastricCancer.xlsx"

############################
# 1.2) LOAD DATA
############################
df_raw <- read_excel(file_path, sheet = 1)

cat("Dataset dimensions:", dim(df_raw), "\n")
cat("Column names:\n")
print(names(df_raw))

############################
# 1.3) CLEANING HELPERS
############################
clean_numeric_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "N/A", "null", "NULL")] <- NA
  x <- gsub("^<\\s*", "", x)
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(x))
}

############################
# 1.4) SELECT ANALYSIS VARIABLES
############################
df <- df_raw %>%
  dplyr::select(
    Gender,
    Age,
    BMI,
    ASA_Score,
    TMN_Stage,  # dataset uses TMN_Stage
    Platelet_Lymphocyte_Ratio,
    Neutrophil_Lymphocyte_Ratio,
    Hemoglobin_PreOp,
    Platelet_Count,
    Albumin_Level,
    Neutrophil_Count,
    Lymphocyte_Count,
    Monocyte_Count,
    WBC_Count,
    Histology_Type,
    Borrmann_Type,
    Pathology_Type,
    Invasion_Depth,
    Invasion_Depth_Grouped,
    Lymph_Nodes_Cleaned,
    Lymph_Nodes_Positive,
    Lymphatic_Invasion,
    Node_Dissection_Extent,
    Peritoneal_Metastasis
  )

############################
# 1.5) FIX DATA TYPES
############################
numeric_cols <- c(
  "Age", "BMI", "Platelet_Lymphocyte_Ratio", "Neutrophil_Lymphocyte_Ratio",
  "Hemoglobin_PreOp", "Platelet_Count", "Albumin_Level",
  "Neutrophil_Count", "Lymphocyte_Count", "Monocyte_Count", "WBC_Count",
  "Lymph_Nodes_Cleaned", "Lymph_Nodes_Positive"
)

for (col in numeric_cols) {
  df[[col]] <- clean_numeric_text(df[[col]])
}

# core categorical variables
df$Gender <- as.factor(df$Gender)
df$ASA_Score <- as.factor(df$ASA_Score)
df$TMN_Stage <- as.factor(df$TMN_Stage)

# coded clinical variables as factors
factor_code_cols <- c(
  "Histology_Type", "Borrmann_Type", "Pathology_Type", "Invasion_Depth",
  "Invasion_Depth_Grouped", "Lymphatic_Invasion", "Node_Dissection_Extent"
)

for (col in factor_code_cols) {
  df[[col]] <- as.factor(df[[col]])
}

# outcome variable
# assumes 0 = No metastasis, 1 = Yes metastasis
df$Peritoneal_Metastasis <- factor(
  df$Peritoneal_Metastasis,
  levels = c(0, 1),
  labels = c("No", "Yes")
)

############################
# 1.6) QUICK DATA CHECK
############################
cat("\nMissing values per column:\n")
print(sort(colSums(is.na(df)), decreasing = TRUE))

cat("\nOutcome counts:\n")
print(table(df$Peritoneal_Metastasis, useNA = "ifany"))

############################
# 1.7) MULTIPLE IMPUTATION (MICE)
############################
df_mice <- df

meth <- make.method(df_mice)
pred <- make.predictorMatrix(df_mice)

# Do not impute the outcome itself
meth["Peritoneal_Metastasis"] <- ""
pred["Peritoneal_Metastasis", ] <- 0
pred[, "Peritoneal_Metastasis"] <- 1
pred["Peritoneal_Metastasis", "Peritoneal_Metastasis"] <- 0

imp <- mice(
  df_mice,
  m = 5,
  maxit = 10,
  method = meth,
  predictorMatrix = pred,
  seed = 1234,
  printFlag = FALSE
)

# Use first completed dataset for analysis
df_complete <- complete(imp, 1)

# Restore factor types after imputation
df_complete$Gender <- as.factor(df_complete$Gender)
df_complete$ASA_Score <- as.factor(df_complete$ASA_Score)
df_complete$TMN_Stage <- as.factor(df_complete$TMN_Stage)
df_complete$Peritoneal_Metastasis <- factor(
  df_complete$Peritoneal_Metastasis,
  levels = c("No", "Yes")
)

for (col in factor_code_cols) {
  df_complete[[col]] <- as.factor(df_complete[[col]])
}

############################
# 1.8) TRAIN / TEST SPLIT
############################
split_obj <- initial_split(df_complete, prop = 0.80, strata = Peritoneal_Metastasis)
train_data <- training(split_obj)
test_data  <- testing(split_obj)

############################
# 1.9) PREPROCESSING RECIPE
############################
# paper mentions normalization before training
rec <- recipe(Peritoneal_Metastasis ~ ., data = train_data) %>%
  step_zv(all_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%
  step_novel(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_numeric_predictors())

############################################################
# 02_baseline_stats.R
############################################################

run_baseline_stats <- function(data, outcome = "Peritoneal_Metastasis") {
  results <- list()
  vars <- setdiff(names(data), outcome)

  for (v in vars) {
    x <- data[[v]]
    y <- data[[outcome]]

    if (is.numeric(x)) {
      tmp <- data.frame(x = x, y = y)
      tmp <- tmp[complete.cases(tmp), ]
      pval <- tryCatch(t.test(x ~ y, data = tmp)$p.value, error = function(e) NA)
      test_type <- "t-test"
    } else {
      tmp <- table(x, y)
      pval <- tryCatch(chisq.test(tmp)$p.value, error = function(e) NA)
      test_type <- "chi-square"
    }

    results[[v]] <- data.frame(
      Variable = v,
      Test = test_type,
      P_value = pval
    )
  }

  bind_rows(results) %>% arrange(P_value)
}

baseline_results <- run_baseline_stats(df_complete)
write.csv(baseline_results, "baseline_stats_results.csv", row.names = FALSE)

############################################################
# 03_correlation.R
############################################################

corr_df <- df_complete %>%
  dplyr::select(
    Age, BMI, Platelet_Lymphocyte_Ratio, Neutrophil_Lymphocyte_Ratio,
    Hemoglobin_PreOp, Platelet_Count, Albumin_Level,
    Neutrophil_Count, Lymphocyte_Count, Monocyte_Count,
    WBC_Count, Lymph_Nodes_Cleaned, Lymph_Nodes_Positive
  )

cor_matrix <- cor(corr_df, use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_matrix, "correlation_matrix.csv")

############################################################
# 04_train_models.R
############################################################

############################
# 4.1) MODEL DEFINITIONS
############################

# Logistic Regression
lr_spec <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

# Decision Tree
dt_spec <- decision_tree(
  mode = "classification",
  cost_complexity = tune(),
  tree_depth = tune(),
  min_n = tune()
) %>%
  set_engine("rpart")

# Random Forest
rf_spec <- rand_forest(
  mode = "classification",
  trees = 500,
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "permutation", probability = TRUE)

# Support Vector Machine
svm_spec <- svm_rbf(
  mode = "classification",
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_engine("kernlab")

# K-Nearest Neighbors
knn_spec <- nearest_neighbor(
  mode = "classification",
  neighbors = tune(),
  weight_func = tune(),
  dist_power = tune()
) %>%
  set_engine("kknn")

# Naive Bayes
nb_spec <- naive_Bayes(
  mode = "classification",
  smoothness = tune(),
  Laplace = tune()
) %>%
  set_engine("klaR")

############################
# 4.2) WORKFLOWS
############################
wf_lr  <- workflow() %>% add_recipe(rec) %>% add_model(lr_spec)
wf_dt  <- workflow() %>% add_recipe(rec) %>% add_model(dt_spec)
wf_rf  <- workflow() %>% add_recipe(rec) %>% add_model(rf_spec)
wf_svm <- workflow() %>% add_recipe(rec) %>% add_model(svm_spec)
wf_knn <- workflow() %>% add_recipe(rec) %>% add_model(knn_spec)
wf_nb  <- workflow() %>% add_recipe(rec) %>% add_model(nb_spec)

############################################################
# 05_cross_validation_tuning.R
############################################################

############################
# 5.1) CROSS-VALIDATION
############################
folds <- vfold_cv(train_data, v = 5, strata = Peritoneal_Metastasis)

############################
# 5.2) METRICS FOR TUNING
############################
my_metrics <- metric_set(roc_auc, accuracy, sens, spec)
ctrl_grid <- control_grid(save_pred = TRUE, save_workflow = TRUE)

############################
# 5.3) TUNING GRIDS
############################
dt_grid <- grid_regular(
  cost_complexity(),
  tree_depth(),
  min_n(),
  levels = 3
)

rf_grid <- grid_regular(
  mtry(range = c(2L, min(10L, ncol(train_data) - 1L))),
  min_n(),
  levels = 3
)

svm_grid <- grid_regular(
  cost(),
  rbf_sigma(),
  levels = 3
)

knn_grid <- grid_regular(
  neighbors(range = c(3L, 25L)),
  dist_power(),
  levels = 3
)

nb_grid <- grid_regular(
  smoothness(),
  Laplace(),
  levels = 3
)

############################
# 5.4) FIT / TUNE MODELS
############################

# Logistic Regression
res_lr <- fit_resamples(
  wf_lr,
  resamples = folds,
  metrics = my_metrics,
  control = control_resamples(save_pred = TRUE)
)

# Decision Tree
res_dt <- tune_grid(
  wf_dt,
  resamples = folds,
  grid = dt_grid,
  metrics = my_metrics,
  control = ctrl_grid
)

# Random Forest
res_rf <- tune_grid(
  wf_rf,
  resamples = folds,
  grid = rf_grid,
  metrics = my_metrics,
  control = ctrl_grid
)

# SVM
res_svm <- tune_grid(
  wf_svm,
  resamples = folds,
  grid = svm_grid,
  metrics = my_metrics,
  control = ctrl_grid
)

# KNN
res_knn <- tune_grid(
  wf_knn,
  resamples = folds,
  grid = knn_grid,
  metrics = my_metrics,
  control = ctrl_grid
)

# Naive Bayes
res_nb <- tune_grid(
  wf_nb,
  resamples = folds,
  grid = nb_grid,
  metrics = my_metrics,
  control = ctrl_grid
)

############################
# 5.5) FINALIZE BEST MODELS
############################
best_dt  <- select_best(res_dt, metric = "roc_auc")
best_rf  <- select_best(res_rf, metric = "roc_auc")
best_svm <- select_best(res_svm, metric = "roc_auc")
best_knn <- select_best(res_knn, metric = "roc_auc")
best_nb  <- select_best(res_nb, metric = "roc_auc")

final_wf_lr  <- wf_lr
final_wf_dt  <- finalize_workflow(wf_dt, best_dt)
final_wf_rf  <- finalize_workflow(wf_rf, best_rf)
final_wf_svm <- finalize_workflow(wf_svm, best_svm)
final_wf_knn <- finalize_workflow(wf_knn, best_knn)
final_wf_nb  <- finalize_workflow(wf_nb, best_nb)

############################
# 5.6) FINAL FIT ON TRAINING DATA
############################
fit_lr  <- fit(final_wf_lr,  data = train_data)
fit_dt  <- fit(final_wf_dt,  data = train_data)
fit_rf  <- fit(final_wf_rf,  data = train_data)
fit_svm <- fit(final_wf_svm, data = train_data)
fit_knn <- fit(final_wf_knn, data = train_data)
fit_nb  <- fit(final_wf_nb,  data = train_data)

############################################################
# 06_model_evaluation.R
############################################################

calc_all_metrics <- function(fitted_workflow, test_df, truth_col = "Peritoneal_Metastasis") {
  probs <- predict(fitted_workflow, test_df, type = "prob")
  preds <- predict(fitted_workflow, test_df, type = "class")

  out <- bind_cols(
    test_df %>% dplyr::select(all_of(truth_col)),
    probs,
    preds
  )

  truth <- out[[truth_col]]
  prob_yes <- out$.pred_Yes
  pred_class <- out$.pred_class

  tibble(
    Accuracy    = accuracy_vec(truth, pred_class),
    Sensitivity = sens_vec(truth, pred_class, event_level = "second"),
    Specificity = spec_vec(truth, pred_class, event_level = "second"),
    PPV         = ppv_vec(truth, pred_class, event_level = "second"),
    NPV         = npv_vec(truth, pred_class, event_level = "second"),
    MSE         = mean((as.numeric(truth == "Yes") - prob_yes)^2, na.rm = TRUE),
    AUC         = roc_auc_vec(truth, prob_yes, event_level = "second")
  )
}

############################
# 6.1) TEST SET PERFORMANCE
############################
perf_lr  <- calc_all_metrics(fit_lr,  test_data) %>% mutate(Model = "LR")
perf_dt  <- calc_all_metrics(fit_dt,  test_data) %>% mutate(Model = "DT")
perf_rf  <- calc_all_metrics(fit_rf,  test_data) %>% mutate(Model = "RF")
perf_svm <- calc_all_metrics(fit_svm, test_data) %>% mutate(Model = "SVM")
perf_knn <- calc_all_metrics(fit_knn, test_data) %>% mutate(Model = "KNN")
perf_nb  <- calc_all_metrics(fit_nb,  test_data) %>% mutate(Model = "NB")

performance_table <- bind_rows(
  perf_dt, perf_knn, perf_svm, perf_nb, perf_rf, perf_lr
) %>%
  dplyr::select(Model, Accuracy, Sensitivity, Specificity, PPV, NPV, MSE, AUC) %>%
  arrange(desc(Accuracy))

print(performance_table)
write.csv(performance_table, "model_performance_testset.csv", row.names = FALSE)

############################
# 6.2) CLASSIFIER COMPARISON PLOT
############################
perf_long <- performance_table %>%
  pivot_longer(
    cols = c(Accuracy, Sensitivity, Specificity, PPV, NPV),
    names_to = "Metric",
    values_to = "Value"
  )

p1 <- ggplot(perf_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_col(position = position_dodge()) +
  labs(
    title = "Comparison of Classifier Performance",
    x = NULL,
    y = "Metric Value"
  ) +
  theme_minimal(base_size = 12)

ggsave("classifier_performance_comparison.png", p1, width = 10, height = 6, dpi = 300)

############################
# 6.3) SAVE CV METRICS
############################
collect_metrics(res_lr)  %>% write.csv("cv_metrics_lr.csv", row.names = FALSE)
collect_metrics(res_dt)  %>% write.csv("cv_metrics_dt.csv", row.names = FALSE)
collect_metrics(res_rf)  %>% write.csv("cv_metrics_rf.csv", row.names = FALSE)
collect_metrics(res_svm) %>% write.csv("cv_metrics_svm.csv", row.names = FALSE)
collect_metrics(res_knn) %>% write.csv("cv_metrics_knn.csv", row.names = FALSE)
collect_metrics(res_nb)  %>% write.csv("cv_metrics_nb.csv", row.names = FALSE)

############################################################
# 07_feature_importance.R
############################################################

rf_fit_obj <- extract_fit_parsnip(fit_rf)$fit

rf_importance <- data.frame(
  Feature = names(rf_fit_obj$variable.importance),
  Importance = as.numeric(rf_fit_obj$variable.importance)
) %>%
  arrange(desc(Importance))

print(rf_importance)
write.csv(rf_importance, "rf_feature_importance.csv", row.names = FALSE)

p2 <- rf_importance %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Random Forest Feature Importance",
    x = "Feature",
    y = "Permutation Importance"
  ) +
  theme_minimal(base_size = 12)

ggsave("rf_feature_importance.png", p2, width = 9, height = 7, dpi = 300)

############################################################
# FINAL SAVE / SESSION INFO
############################################################
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

cat("\nDone. Files created:\n")
cat("- baseline_stats_results.csv\n")
cat("- correlation_matrix.csv\n")
cat("- model_performance_testset.csv\n")
cat("- classifier_performance_comparison.png\n")
cat("- rf_feature_importance.csv\n")
cat("- rf_feature_importance.png\n")
cat("- cv_metrics_lr.csv\n")
cat("- cv_metrics_dt.csv\n")
cat("- cv_metrics_rf.csv\n")
cat("- cv_metrics_svm.csv\n")
cat("- cv_metrics_knn.csv\n")
cat("- cv_metrics_nb.csv\n")
cat("- sessionInfo.txt\n")
