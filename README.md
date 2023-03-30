# superconductors

This repository contains codes related to the publication "Optimal composition based material descriptors for machine learning: A case study on superconductor classification" (add link arXiv).

In particular:

* Folder ```Classification``` contains code for training/validation/testing of all classifiers used in this work (ETCs, QEGs, Naive Bayesian), for assessing their performances (```Classifiers_performances.ipynb```), and the output file with such performances (```classifiers_comparison.xlsx```).
* Folder ```Mixed features optimization``` contains all ```.m``` files for finding the optimized mixed features with multi-objective optimization as linear or product combination of the original features extracted by means of Matminer. Specifically, the file ```MAIN.m``` has to be run, deciding how many features to mix and how many mixed features to have in output (1 or 2); on the contrary, ```MainSingle.m``` has to be run for finding mixed features with single-objective optimization. Folder ```Pareto fronts``` contains already computed Pareto fronts for the examples shown in this work.
* Folder ```Regression & invariance``` contains the file ```ETR&SHAP.ipynb``` for training/validating/testing the ETR model, with the SHAP analysis to rank the input features, together with such ranking ```SHAP_for_ETR_metallic_mean.xlsx``` and with the code for the search of invariant groups ```DNN&invariant_groups.ipynb```.
* File ```Coefficients_mixed_variables.xlsx``` contains the coefficients for mixing the first 30 or 52 original features extracted by means of Matminer in the ranking obtained with SHAP.
* File ```Database_construction.ipynb``` contains the code for cleaning the original SuperCon database.
* File ```predictions_on_MPj_materials.xlsx``` contains the probability predictions of the best two classifiers of this work (ETC-vanilla, ETC-SMOTE), the best QEG-based classifier (QEG 2D-mixed lin) and of the GEV classifier (for $T_{\rm{c}}\geq 35~\textup{K}$) over the $\sim$ 40,000 materials in MaterialsProject and not in SuperCon; furthermore, the class prediction 1/0 is also provided, considering the probability threshold maximizing the $F_{1, \textup{max}}$ score.
