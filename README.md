# CMS Classification
The CMS classification can be done using the CMSCaller Function shown in the CMSPipe.R Script.

# Enrichment scoring
The scoring of enrichment of various signatures in Bulk samples(Microarray/Rnaseq) of a dataset is done using the ssGSEA function in the GSEApy module in Python. The Signatures.xlsx file contains the various signatures we used- Epithelial, Mesenchymal, FAO, Glycolysis, OXPHOS, PD-L1, and the CMS subtypes 1 to 4. The enrichment scores for single-cell RNAseq data are calculated using the AUCell package in R.

# Entropy Calculation
The Entropy of a subset of genes (Epithelial, mesenchymal, Metabolic pathways, etc) is taken, and the heterogeneity/variability of their expression in many cells of a dataset is calculated using the Shannon entropy formula.
