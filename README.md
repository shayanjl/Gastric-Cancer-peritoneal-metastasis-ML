# Peritoneal Metastasis Prediction in Gastric Cancer


This repository contains an **R-based reconstructed machine learning analysis pipeline** developed to reflect the main analytical workflow of the study:

**Peritoneal Metastasis Prediction in Gastric Cancer: A Machine Learning Approach**

## Overview

The goal of this repository is to provide a clear and organized implementation of the main analysis steps described in the study, including:

- data preprocessing
- baseline statistical analysis
- correlation analysis
- machine learning model training
- cross-validation and hyperparameter tuning
- model evaluation
- feature importance analysis

The repository is intended to support transparency, reproducibility, and future extension of the analytical workflow.

## Important Note

This repository **does not contain the exact final scripts used in the manuscript submission**.

The code provided here is a **carefully reconstructed public version** of the workflow based on the study design, dataset structure, and reported analytical methods. Some implementation details, preprocessing decisions, and publication-specific materials are intentionally not identical to the final internal analysis version.

As the study is currently under publication process, the **full original scripts and finalized publication-matched materials will be released after formal publication of the article**.

## What is included

This repository currently includes a structured implementation of the following main analytical components:

1. **Preprocessing**  
   Reading the dataset, selecting variables, cleaning values, handling missing data, type conversion, imputation, and normalization.

2. **Baseline Statistical Analysis**  
   Group comparisons using statistical tests for clinical and laboratory variables.

3. **Correlation Analysis**  
   Assessment of relationships among selected continuous variables.

4. **Model Training**  
   Training multiple machine learning classifiers for prediction of peritoneal metastasis.

5. **Cross-Validation and Hyperparameter Tuning**  
   Five-fold cross-validation and model tuning for improved robustness.

6. **Model Evaluation**  
   Evaluation using classification metrics including accuracy, sensitivity, specificity, positive predictive value, negative predictive value, mean squared error, and area under the ROC curve.

7. **Feature Importance Analysis**  
   Variable importance analysis, particularly from the Random Forest model.

## Models included

The reconstructed pipeline includes the following classifiers:

- Logistic Regression
- Decision Tree
- Random Forest
- Support Vector Machine
- K-Nearest Neighbors
- Naive Bayes

## Reproducibility statement

This repository is intended to reflect the general computational strategy of the study. However, because it is a reconstructed public version rather than the exact publication codebase, some numerical outputs may differ from the final manuscript results.

The definitive publication-aligned code and associated reproducibility materials will be made available after publication.

## Repository status

**Current status:** pre-publication public reconstruction

This repository is shared to document the analytical framework while reserving the exact final publication scripts and certain project-specific materials until completion of the publication process.

## Future release

After publication, this repository will be updated with:

- the final publication-aligned code
- clarified implementation details where appropriate
- additional documentation for exact reproducibility
- any approved publication-linked analysis materials

## Citation

If you use or reference this repository, please cite the associated article once published.


## Public dataset

The dataset supporting this study is publicly available in the **BioStudies** database under accession:

**S-EPMC5383064**

Dataset page:  
https://www.ebi.ac.uk/biostudies/studies?query=S-EPMC5383064

## Contact

**Shayan Jalali**  
Department of Health Sciences (DISS), University of Eastern Piedmont (UPO), Novara, Italy  
Email: shayanjalali.bioinformatics@gmail.com
