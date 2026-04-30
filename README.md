# Single-Cell Transcriptomic Profiling of Progressive Mucinous Histiocytosis (PMH)

## 📌 Project Overview
This repository contains the computational pipeline and bioinformatics analysis for single-cell RNA sequencing (scRNA-seq) data derived from Progressive Mucinous Histiocytosis (PMH). 

The primary goal of this project is to map the cellular landscape of PMH, identify rare pathogenic cell states, and uncover potential disease-driving mechanisms through differential expression and pathway analysis.

## 🧬 Computational Workflow
This analysis was built using R and Scanpy and follows a standardized single-cell pipeline:

1. **Quality Control & Pre-processing:** Filtering of low-quality cells, mitochondrial content removal, and doublet detection using [e.g., DoubletFinder / scDblFinder].
2. **Normalization & Feature Selection:** [e.g., SCTransform or standard LogNormalize].
3. **Dimensionality Reduction & Integration:** PCA and UMAP generation. Batch effects were corrected using [e.g., Harmony / scVI / Seurat Integration].
4. **Clustering & Cell Type Annotation:** Graph-based clustering and manual annotation using canonical marker genes.
5. **Downstream Analysis:** * Differential Gene Expression (DGE) 
   * Pathway Enrichment Analysis (Gene Ontology / KEGG)

