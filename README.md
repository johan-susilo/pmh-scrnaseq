# Single-Cell Transcriptomic Profiling of Progressive Mucinous Histiocytosis (PMH)

This repository contains an automated, end-to-end Snakemake pipeline for the computational analysis of single-cell RNA sequencing (scRNA-seq) data. The pipeline is designed to process raw count matrices through quality control, integration, annotation, differential expression, trajectory inference, and cell-cell communication.

## 📁 Repository Structure

```text
.
├── config/
│   ├── config.yaml          # Main pipeline configuration (parameters, cutoffs, features)
│   └── input.csv            # Sample metadata (Sample ID, Condition/Group mappings)
├── workflow/
│   ├── Snakefile            # The Snakemake workflow definition
│   ├── envs/
│   │   └── scRNA_env.yaml   # Conda/Mamba environment dependencies
│   └── scripts/             # R scripts executed by Snakemake
│       ├── 00_utils.R               # Shared helpers, colors, and ggplot themes
│       ├── 01_preprocessing.R       # QC, DoubletFinder, Normalization, Harmony Integration
│       ├── 02_global_annotation.R   # SingleR, scCATCH, CelliD, and Consensus Voting
│       ├── 02b_global_dge.R         # Global-level Differential Gene Expression (DESeq2)
│       ├── 02c_global_go.R          # Global-level Pathway Analysis (clusterProfiler)
│       ├── 03_subset_clusters.R     # Isolating specific cell types (e.g., Fibroblasts)
│       ├── 04_detail_annotation.R   # High-resolution sub-clustering and Module Scoring
│       ├── 05_dge.R                 # Sub-cluster level DGE (Pseudobulk DESeq2)
│       ├── 06_go.R                  # Sub-cluster level GO/KEGG pathway enrichment
│       ├── 07_pseudotime.R          # Trajectory inference using Monocle3
│       ├── 08_cellchat_analysis.R   # Multi-way cross-talk for specific subsets
│       └── 08b_global_cellchat.R    # Macro-level cell-cell communication networks
└── README.md
```

---

## ⚙️ Setup & Installation

**1. Clone the repository and navigate to the directory:**
```bash
git clone <your-repo-url> skin_snakemake
cd skin_snakemake
```

**2. Install dependencies using Conda/Mamba:**
The pipeline relies on a tightly controlled Conda environment containing R 4.3+, Seurat v5, Monocle3, CellChat, and Bioconductor packages.
```bash
mamba env create -f workflow/envs/scRNA_env.yaml
```

*Note: Snakemake can also handle environment creation automatically if you run it with the `--use-conda` flag.*

---

## 🛠️ Configuration

Before running the pipeline, configure your inputs in the `config/` directory.

### 1. `config/input.csv`
Defines the samples to be processed. The pipeline expects three columns:
*   `sample_names`: The exact folder name of the raw data (e.g., CellRanger output) in your `data_dir`.
*   `ident1`: Broad condition (e.g., `Healthy_Fresh`, `Acute_PMH`).
*   `ident2`: Specific sample identifier (e.g., `HTY131`, `AC259`).

### 2. `config/config.yaml`
The central control hub for the pipeline. Key sections include:
*   **I/O:** `outdir` and `data_dir`.
*   **Preprocessing:** Mitochondrial cutoffs, doublet rates, PCs to calculate, and resolution.
*   **Subsetting:** Which cell types (e.g., `fibroblasts`, `macrophages`) to isolate for deep-dive analysis (`03` through `07`).
*   **DGE & Pathways:** LogFC and p-value cutoffs, plus the metadata column to compare (`Condition`).
*   **CellChat:** Defines pairs/groups of cell types to evaluate for communication networks.

---

## 🚀 Running the Pipeline

Activate your base Snakemake environment (or the `seurat_env` if Snakemake is installed there), then execute the workflow.

**Dry Run (Check what will be executed without actually running it):**
```bash
snakemake -n
```

**Run Locally (e.g., using 8 cores):**
```bash
snakemake --cores 8
```

**Run specific steps (Targeting rules):**
You can tell Snakemake to run up to a specific rule. For example, to just run preprocessing and global annotation:
```bash
snakemake --cores 8 global_annotation
```

---

## 📊 Pipeline Stages & Outputs

By default, everything is saved into the `results/` folder (or whatever `outdir` is set to in `config.yaml`).

| Stage | Script | Description | Primary Output |
| :--- | :--- | :--- | :--- |
| **1. QC & Preprocessing** | `01_preprocessing.R` | Reads 10x matrices, filters dead cells/doublets, normalizes (SCT/Log), and integrates via Harmony. | `01_preprocessing/TN.combined_dim30.rds` |
| **2. Global Annotation** | `02_global_annotation.R` | Automated consensus annotation mapping cells to broad lineages. | `02_global_annotation/res_*/TN.combined_annotated.rds` |
| **3. Global DGE/GO** | `02b_` & `02c_` | Compares Healthy vs PMH globally (Pseudobulk). | `02_global_annotation/res_*/pathways/` |
| **4. Subsetting** | `03_subset_clusters.R` | Isolates specific cell types defined in `config.yaml` into independent objects. | `03_subsets/<celltype>/00_data/*_subset_processed.rds` |
| **5. Detailed Annotation** | `04_detail_annotation.R` | Applies high-res module scoring (e.g., M1/M2 macrophages, F1-F8 fibroblasts). | `03_subsets/<celltype>/00_data/*_detailed_annotated.rds` |
| **6. Subset DGE/GO** | `05_dge.R` & `06_go.R` | Pseudobulk DESeq2 and clusterProfiler enrichment per sub-cluster. | `03_subsets/<celltype>/04_pathways/` |
| **7. Trajectory** | `07_pseudotime.R` | Monocle3 pseudotime ordering and trajectory-dependent gene mapping. | `03_subsets/<celltype>/05_pseudotime/` |
| **8. Communication** | `08_cellchat*.R` | Maps receptor-ligand interactions dynamically across conditions (Healthy vs PMH). | `04_cellchat/` |

## 📝 Logging and Debugging
Every rule in the Snakefile outputs logs to `results/logs/`. If a step fails, check the respective `.log` file in that directory for standard R error messages.