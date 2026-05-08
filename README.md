# 🧠 RNA-Seq Analysis: Gene Expression in Autism Spectrum Disorder (ASD)
 
> Transcriptomic profiling of human brain cortex samples to uncover novel genes, regulators, and pathways associated with ASD.
 
---
 
## 📌 Overview
 
This project reproduces and expands upon the analysis from:
 
> **"Comprehensive Analysis of RNA-Seq Gene Expression Profiling of Brain Transcriptomes Reveals Novel Genes, Regulators, and Pathways in Autism Spectrum Disorder"**
> [PMC7603078](https://pmc.ncbi.nlm.nih.gov/articles/PMC7603078/)
 
The original study aimed to profile gene expression signatures in the cerebral cortex of ASD patients using two publicly available RNA-seq datasets. 
My analysis focused specifically on **identifying differentially expressed genes and classifying them by biological role** - separating them into **Neuronal**, **Inflammatory**, and **Mixed** categories.
 
---
 
## 🗂️ Dataset
 
**Source:** [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) - transcriptomics public database
 
**Dataset:** [GSE64018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64018)
 
| Parameter | Details |
|-----------|---------|
| Samples | 12 ASD + 12 healthy controls |
| Matching | Paired by age and sex |
| Tissue | Post-mortem brain — Superior Temporal Gyrus (STG) |
| Brain Areas | Brodmann areas BA41/42/22 |
| Dissection | All cortical layers, gray matter preserved |
| Dimensions | 63,152 features × 24 samples |
 
---
 
## 🔬 Analysis Pipeline
 
### 1. 📊 Differential Expression — DESeq2
Applied the **DESeq2** framework to identify statistically significant differentially expressed genes (DEGs) between ASD and control samples. Normalized read counts and applied shrinkage estimators to improve fold-change accuracy.
 
### 2. 📉 PCA — All Samples
Principal Component Analysis (PCA) performed on all 24 samples (ASD + controls) to assess **overall transcriptomic variance**, detect batch effects, and confirm sample grouping.
 
### 3. 📉 PCA — ASD Samples Only
A second PCA was conducted **exclusively on the 12 ASD samples** to investigate internal heterogeneity within the ASD group and identify potential transcriptomic subgroups.
 
### 4. 🔍 Notable Gene Extraction
Highlighted genes with the most striking expression profiles — considering both statistical significance (adjusted p-value) and biological magnitude (log2 fold-change). These were prioritized for downstream interpretation.
 
### 5. 🗃️ Cluster Structure Analysis
Explored the organization of gene co-expression clusters, examining how DEGs group together and what shared biological processes may be driving their co-regulation.
 
### 6. 🏷️ Gene Classification by Marker Type
Classified differentially expressed genes into three functional categories:
 
| Category | Description |
|----------|-------------|
| 🟦 **Neuronal** | Genes associated with synaptic function, neuronal signaling, axon guidance, and neural development |
| 🟥 **Inflammatory** | Genes linked to immune activation, microglial response, cytokine signaling, and neuroinflammation |
| 🟨 **Mixed** | Genes with roles in both neuronal and inflammatory contexts, or with pleiotropic functions |
 
### 7. ⬆️⬇️ Activated vs. Repressed Genes
Identified and separated **upregulated** and **downregulated** genes in ASD relative to controls, with a focus on the most biologically relevant candidates in each direction.
 
### 8. 🧬 GO Enrichment Analysis
Performed **Gene Ontology (GO) Enrichment Analysis** across Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) categories to identify overrepresented functional terms among the DEGs.
 
### 9. 🎯 Top Driver Analysis
Investigated the **top driver genes** — those contributing most to the transcriptomic variance and cluster separation — to understand which molecular players may be central to ASD-related dysregulation.
 
---
 
## 🧰 Tools & Technologies
 
- **R** - DESeq2, ggplot2, clusterProfiler, pheatmap
- **Bioconductor** - Bioinformatics package ecosystem
- **GEO / NCBI** - Data acquisition
- **Gene Ontology (GO)** - Functional annotation and enrichment
---
 
## 📁 Repository Structure
 
```
📦 rnaseq-asd-analysis/
├── data/                  # Raw and processed count matrices
├── scripts/               # R scripts for each analysis step
│   ├── transcriptomic_Analysis.R
├── results/               # Output tables and figures
├── figures/               # Publication-ready plots
└── README.md
```
 
---
 
## 📚 Reference
 
Xiong K, Green S, Cao Q, et al. *Comprehensive Analysis of RNA-Seq Gene Expression Profiling of Brain Transcriptomes Reveals Novel Genes, Regulators, and Pathways in Autism Spectrum Disorder.* Front Psychiatry. 2020;11:587465. [doi:10.3389/fpsyt.2020.587465](https://doi.org/10.3389/fpsyt.2020.587465)
 
---
 
*Analysis performed as part of a personal bioinformatics portfolio.*
 
