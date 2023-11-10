# PRJNA883411

The goals of the exploartory analysis were to characterize cell types and states in this scRNA-seq dataset of human female prefrontal cortex samples. 

## Key steps included:

 - Performing quality control on raw count data from Cell Ranger preprocessing of single-nuclei samples
 - Estimate ambient RNA fractions in nuclei across different samples and correct expression matrices
 - Indentify doublets and exclude them from the analysis
 - Deriving thresholds for filtering out low quality nuclei
 - Transform count matrices and select highly variable genes for further analysis
 - Running dimensionality reduction and unsupervised clustering of high-quality nuclei
 - Reconciling clusters at multiple resolution levels into a clustering tree to define robust cell groups
 - Identifying differentially expressed genes between defined cell groups using UMIs count matrices

The final analysis aims to characterize cell types and transcriptional states in the prefrontal cortex of females with MDD compared to healthy female controls.
