# *Dissecting the intra- and intertumoral heterogeneity of adrenocortical carcinoma by single-cell multi-omics analyses*

This page recorded the codes and data used and mentioned in [*xxx*](XXX). And you could downloaded this paper by clicking [here](pdf/XXX)

![image-20231220102432789](D:\xiangyu.ubuntu\code\log-summery\code.public\Dissecting-heterogeneity-of-adrenocortical-carcinoma-by-single-cell-multi-omics-analyses\README.assets\image-20231220102432789.png)

Adrenocortical carcinoma (ACC) is a rare but aggressive malignancy originating in the adrenal cortex, characterized by significant intra- and intertumoral heterogeneity. In this study, we employed a combination of whole-genome sequencing, single-cell RNA sequencing, T cell receptor sequencing, spatial transcriptome sequencing, and multiple fluorescent staining to construct a comprehensive multi-omics landscape of ACC. Our findings demonstrated that, although all tumor cells exhibited a "confused cell identity" feature, distinct subpopulations were present in each ACC patient. Among these, the *LDLR*+ and *MKI67*+ subpopulations expressed elevated levels of C1A signature genes and established robust communication with endothelial cells through *NECTIN*-*PVR* and *NOTCH2*-*JAG* pairs, both significantly associated with poor prognosis. Interestingly, the *PCDH15*+ subpopulations displayed moderate levels of both C1A signature genes and extracellular matrix related genes. On the other hand, the *HLA-B*+ and subpopulation exhibited the highest levels of antigen-processing related genes, a characteristic not previously reported. The amalgamation of these distinct subpopulations facilitated the stratification of ACC patients into subtypes with varying prognoses. Patients (Group II) dominated by *LDLR*+ and *MKI67*+ subpopulations exhibited the worst prognosis, while those (Group III) with a higher proportion of *PCDH15*+ subpopulation demonstrated the most favorable outcomes. Significantly, patients (Group I) with a higher prevalence of the *HLA-B*+ subpopulation also showed notably increased *LAG3*+ memory CD8 T cells, indicating their potential suitability for immunotherapy. Despite responding averagely to conventional treatment, this subgroup was predicted to be responsive to immunotherapy. Notably, one patient from Group I exhibited a positive response to camrelizumab treatment after relapse with mitotane, leading to a successful outcome, while two patients from Group II did not respond favorably. Our study offers insights into the single-cell level heterogeneity of ACC and provides valuable implications for precision treatments for this disease.

To effectively demonstrate our step-by-step analysis of single-cell RNA sequencing (scRNA-seq), single-cell T-cell receptor sequencing (scTCR-seq), and single-cell spatial transcriptomics (scSpatial), we have meticulously compiled and stored detailed procedural information in Markdown files. This extensive documentation encompasses records of quality control measures, batch effect reduction, dimensionality reduction, cell clustering, pseudotime construction, identification of dynamically expressed genes, and pathway enrichment.

# **1. Codes of analyzing and visualization**

**Introduction to Our Script Compilation for ACC Analysis**

In our comprehensive analysis of adrenocortical carcinoma (ACC), we have organized our scripts into 10 distinct chapters. Each chapter focuses on a specific aspect of the analysis, allowing for detailed examination and understanding of the different facets of ACC. Below is a guide to where you can find each part in the respective chapters:

- **[Chapter 0](Chapter0.md): Pre-processing of Single-cell Sequencing Data**
  - Content: Initial processing steps for scRNA-seq, scTCR-seq, and scSpatial datasets, including data cleaning, normalization, and quality control.
- **[Chapter 1](Chapter1.md): scRNA-seq Analysis of Normal Adrenal and ACC Landscape**
  - Content: Scripts for analyzing scRNA-seq data to compare the cellular landscapes of normal adrenal glands and ACC tumors.
- **[Chapter 2](Chapter2.md): Developmental Analysis of Normal Adrenal Cortex**
  - Content: Computational workflows for studying the developmental processes in the normal adrenal cortex.
- **[Chapter 3](Chapter3.md): Identification of 'Confused Cell Identity' in ACC**
  - Content: Code for analyzing the mixed or confused cell identity phenomenon within ACC cells.
- **[Chapter 4](Chapter4.md): ACC Intra-Tumor Heterogeneity Analysis**
  - Content: Scripts for dissecting the heterogeneity within ACC tumors at the cellular level.
- **[Chapter 5](Chapter5.md): Molecular Classification of ACC**
  - Content: Code for categorizing ACC into different molecular subtypes based on single-cell data.
- **[Chapter 6](Chapter6.md): Intra-tumoral Cell-Cell Interaction Analysis**
  - Content: Computational methods for exploring cell-cell interactions within ACC tumors.
- **[Chapter 7](Chapter7.md): Tumor-Endothelial Cell Interaction in ACC**
  - Content: Scripts focused on the interactions between tumor cells and endothelial cells in the ACC microenvironment.
- **[Chapter 8](Chapter8.md): Tumor-Myeloid Cell Interaction in ACC**
  - Content: Analysis code for studying interactions between tumor cells and myeloid cells in ACC.
- **[Chapter 9](Chapter9.md): Tumor-T Cell Interaction in ACC**
  - Content: Scripts for investigating the interactions between tumor cells and T cells within the ACC context.
- **[Chapter 1](Chapter10.md)0: WES Analysis Pipeline**
  - Content: Development and implementation of a pipeline for Whole Exome Sequencing (WES) analysis, focusing on genomic alterations in ACC tumors.

Each chapter contains detailed scripts, methodologies, and analyses pertinent to the specific aspect of ACC it addresses. This structured approach allows researchers to navigate through our comprehensive analysis with ease, enhancing the understanding and study of ACC.

# **2. Raw data download**

- **Description**: This section includes all the raw FASTQ files from our study. These files are crucial for in-depth data analysis and understanding the sequencing results from single-cell RNA, T-cell receptor, and spatial transcriptomics.

- **Download**: You can access and download these files from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
[4.0K]  .
├── [4.0K]  scRNA
│   ├── [ 44G]  ACC10_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 40G]  ACC10_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 36G]  ACC1_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 33G]  ACC1_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 27G]  ACC2_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 25G]  ACC2_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 29G]  ACC3_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 26G]  ACC3_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 29G]  ACC4_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 27G]  ACC4_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 45G]  ACC5_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 35G]  ACC5_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 53G]  ACC6_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 42G]  ACC6_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 56G]  ACC7_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 44G]  ACC7_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 47G]  ACC8_RNA_S1_L001_R1_001.fastq.gz
│   ├── [ 36G]  ACC8_RNA_S1_L001_R2_001.fastq.gz
│   ├── [ 52G]  ACC9_RNA_S1_L001_R1_001.fastq.gz
│   └── [ 43G]  ACC9_RNA_S1_L001_R2_001.fastq.gz
├── [4.0K]  scTCR
│   ├── [4.6G]  ACC1_TCR_S1_L001_R1_001.fastq.gz
│   ├── [4.2G]  ACC1_TCR_S1_L001_R2_001.fastq.gz
│   ├── [3.3G]  ACC2_TCR_S1_L001_R1_001.fastq.gz
│   ├── [2.8G]  ACC2_TCR_S1_L001_R2_001.fastq.gz
│   ├── [3.2G]  ACC3_TCR_S1_L001_R1_001.fastq.gz
│   ├── [3.0G]  ACC3_TCR_S1_L001_R2_001.fastq.gz
│   ├── [3.6G]  ACC4_TCR_S1_L001_R1_001.fastq.gz
│   └── [3.1G]  ACC4_TCR_S1_L001_R2_001.fastq.gz
└── [4.0K]  spatial
    ├── [ 24G]  ACC11_S1_L001_R1_001.fastq.gz
    ├── [ 25G]  ACC11_S1_L001_R2_001.fastq.gz
    ├── [ 24G]  ACC12_S1_L001_R1_001.fastq.gz
    ├── [ 23G]  ACC12_S1_L001_R2_001.fastq.gz
    ├── [ 26G]  ACC13_S1_L001_R1_001.fastq.gz
    ├── [ 24G]  ACC13_S1_L001_R2_001.fastq.gz
    ├── [ 36G]  ACC14_S1_L001_R1_001.fastq.gz
    ├── [ 34G]  ACC14_S1_L001_R2_001.fastq.gz
    ├── [ 10G]  Normal1_spatial_S1_L001_R1_001.fastq.gz
    ├── [ 28G]  Normal1_spatial_S1_L001_R2_001.fastq.gz
    ├── [9.6G]  Normal2_spatial_S1_L001_R1_001.fastq.gz
    ├── [ 28G]  Normal2_spatial_S1_L001_R2_001.fastq.gz
    ├── [ 12G]  Normal3_spatial_S1_L001_R1_001.fastq.gz
    └── [ 34G]  Normal3_spatial_S1_L001_R2_001.fastq.gz
3 directories, 42 files
```

# **3. Processed Data Download**

## 3.1. CellRanger and SpaceRanger Output

- **Description**: This section includes the output files from Cell Ranger and Space Ranger, essential for the initial data processing and analysis of single-cell RNA, T-cell receptor, and spatial transcriptomics data.
- **Download**: These files are available for access and download from the [GEO database](https://chat.openai.com/c/link-to-download).

Below is a detailed annotation of the file structure and contents:

```shell
plaintextCopy code[4.0K]  .
├── [4.0K]  cellranger_output
│   ├── [424M]  ACC10_RNA.tar.gz
│   ├── [2.2G]  ACC11_spatial.tar.gz
│   ├── [1.9G]  ACC12_spatial.tar.gz
│   ├── [2.1G]  ACC13_spatial.tar.gz
│   ├── [2.2G]  ACC14_spatial.tar.gz
│   ├── [281M]  ACC1_RNA.tar.gz
│   ├── [855K]  ACC1_TCR.tar.gz
│   ├── [195M]  ACC2_RNA.tar.gz
│   ├── [1.9M]  ACC2_TCR.tar.gz
│   ├── [148M]  ACC3_RNA.tar.gz
│   ├── [3.1M]  ACC3_TCR.tar.gz
│   ├── [174M]  ACC4_RNA.tar.gz
│   ├── [2.0M]  ACC4_TCR.tar.gz
│   ├── [315M]  ACC5_RNA.tar.gz
│   ├── [359M]  ACC6_RNA.tar.gz
│   ├── [322M]  ACC7_RNA.tar.gz
│   ├── [225M]  ACC8_RNA.tar.gz
│   ├── [408M]  ACC9_RNA.tar.gz
│   ├── [115M]  Normal1_RNA.tar.gz
│   ├── [1.1G]  Normal1_spatial.tar.gz
│   ├── [142M]  Normal2_RNA.tar.gz
│   ├── [1.1G]  Normal2_spatial.tar.gz
│   ├── [322M]  Normal3_RNA.tar.gz
│   └── [1.0G]  Normal3_spatial.tar.gz
```

- Contents
  - Each `_RNA.tar.gz` file includes the filtered_feature_bc_matrix output and loupe file from the Cell Ranger count model.
  - Each `_TCR.tar.gz` file contains the filtered_contig_annotations.csv, clonotypes.csv, and loupe.vloupe files generated from the Cell Ranger TCR model.
  - Each `_spatial.tar.gz` file includes the filtered_feature_bc_matrix output and loupe.cloupe file from the Space Ranger count model.

## 3.2. R Data Files Generated in This Study

- **Description**: All R data files (.rds) related to single-cell RNA sequencing (scRNA-seq), single-cell T-cell receptor sequencing (scTCR-seq), single-cell spatial transcriptomics (scSpatial), and The Cancer Genome Atlas Adrenocortical Carcinoma (TCGA-ACC) data are available. These files encompass a comprehensive range of analyses and findings from our study.
- **Download**: You can download these files from [Zenodo Zenodo10416598](https://zenodo.org/uploads/10416598).

Here's an annotation for each file to give you :

~~~shell
tree -lh
[4.0K]  .
├── [ 16M]  CCI.cor_genes.GSEA.rds
CCI.cor_genes.GSEA.rds: Gene Set Enrichment Analysis (GSEA) results for genes correlated with Confused Cell Identity (CCI) in ACC.

├── [791K]  CCI.cor_genes.rds
CCI.cor_genes.rds: List or data frame of genes correlated with Confused Cell Identity in ACC.

├── [534K]  TCGA.ACC_GSVA.Differences.rds
TCGA.ACC_GSVA.Differences.rds: Differential Gene Set Variation Analysis (GSVA) scores for ACC samples from The Cancer Genome Atlas (TCGA) dataset.

├── [ 16M]  TCGA.ACC_GSVA.scores.rds
TCGA.ACC_GSVA.scores.rds: GSVA scores for ACC samples from TCGA.

├── [6.0K]  TCGA.ACC_clinical.classify.rds
TCGA.ACC_clinical.classify.rds: Clinical classification data for ACC patients from TCGA.

├── [306M]  TCGA.ACC_decon.egm.rds
TCGA.ACC_decon.egm.rds: Results from deconvolution analysis of gene expression matrix in ACC from TCGA.

├── [ 11M]  TCGA.ACC_exp.rds
TCGA.ACC_exp.rds: Expression data for ACC samples from TCGA.

├── [ 11M]  TCGA.ACC_exp_log.rds
TCGA.ACC_exp_log.rds: Log-transformed expression data for ACC samples from TCGA.

├── [1.6G]  scRNA.ACC.and.Normal.adrenal.merge.only.Cortex.rds
scRNA.ACC.and.Normal.adrenal.merge.only.Cortex.rds: Merged scRNA-seq data focusing only on the cortex of ACC and normal adrenal samples.

├── [1.0G]  scRNA.ACC.and.Normal.adrenal.merge.only.Cortex_downsample_exp.rds
scRNA.ACC.and.Normal.adrenal.merge.only.Cortex_downsample_exp.rds: Downsampled expression data from the cortex of merged ACC and normal adrenal samples.

├── [490K]  scRNA.ACC.and.Normal.adrenal.merge.only_Endo.PVR.pos_vs_OTS.genes.rds
scRNA.ACC.and.Normal.adrenal.merge.only_Endo.PVR.pos_vs_OTS.genes.rds: Gene expression comparison in endothelial cells (Endo) expressing PVR in ACC and normal samples.

├── [ 26K]  scRNA.ACC.and.Normal.adrenal.merge.only_Endo.PVR.pos_vs_OTS_GO.rds
scRNA.ACC.and.Normal.adrenal.merge.only_Endo.PVR.pos_vs_OTS_GO.rds: Gene Ontology analysis for Endo cells expressing PVR in ACC versus other cell types (OTS).

├── [ 24M]  scRNA.ACC.and.Normal.adrenal.merge.only_Endo.rds
scRNA.ACC.and.Normal.adrenal.merge.only_Endo.rds: scRNA-seq data of endothelial cells from merged ACC and normal adrenal samples.

├── [135M]  scRNA.ACC.and.Normal.adrenal.merge.only_Myeloid.rds
scRNA.ACC.and.Normal.adrenal.merge.only_Myeloid.rds: scRNA-seq data of myeloid cells from merged ACC and normal adrenal samples.

├── [ 19M]  scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.LAG3.CD8Tm_vs_OTS.T.GSEA.rds
scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.LAG3.CD8Tm_vs_OTS.T.GSEA.rds: GSEA results for NKT and LAG3+ CD8Tm cells in ACC compared to other T cell types.

├── [515K]  scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.LAG3.CD8Tm_vs_OTS.T.genes.rds
scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.LAG3.CD8Tm_vs_OTS.T.genes.rds: Gene expression data for NKT and LAG3+ CD8Tm cells in ACC compared to other T cell types.

├── [345M]  scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.rds
scRNA.ACC.and.Normal.adrenal.merge.only_NKT_T.rds: scRNA-seq data of NKT cells from merged ACC and normal adrenal samples.

├── [2.0G]  scRNA.ACC.and.Normal.adrenal.merge.seurat.rds
scRNA.ACC.and.Normal.adrenal.merge.seurat.rds: Seurat object containing merged scRNA-seq data from ACC and normal adrenal samples.

├── [698M]  scRNA.ACC.cellchat.rds
scRNA.ACC.cellchat.rds: CellChat analysis results for cell-cell communication in ACC samples.

├── [389K]  scRNA.ACC.only.Cortex_DEGs.rds
scRNA.ACC.only.Cortex_DEGs.rds: Differentially expressed genes in the cortex of ACC samples.

├── [ 24K]  scRNA.ACC.only.Cortex_KEGG.rds
scRNA.ACC.only.Cortex_KEGG.rds: KEGG pathway analysis results for the cortex of ACC samples.

├── [ 78M]  scRNA.ACC.only.Cortex_downsamples.rds
scRNA.ACC.only.Cortex_downsamples.rds: Downsampled scRNA-seq data from the cortex of ACC samples.

├── [1.5G]  scRNA.ACC.only.Cortex_harmony.rds
scRNA.ACC.only.Cortex_harmony.rds: Harmony algorithm results for integrating scRNA-seq data from the cortex of ACC samples.

├── [2.0G]  scRNA.ACC.seurat.rds
scRNA.ACC.seurat.rds: Seurat object containing scRNA-seq data for ACC samples.

├── [ 11K]  scRNA.Normal.adrenal.only.Cortex.and.Medulla_signatures.rds
scRNA.Normal.adrenal.only.Cortex.and.Medulla_signatures.rds: Gene expression signatures for the cortex and medulla of normal adrenal glands.

├── [1.1G]  scRNA.Normal.adrenal.only.Cortex_URD.rds
scRNA.Normal.adrenal.only.Cortex_URD.rds: URD analysis results for the cortex of normal adrenal samples.

├── [604M]  scRNA.Normal.adrenal.only.Cortex_downsample_exp.rds
scRNA.Normal.adrenal.only.Cortex_downsample_exp.rds: Downsampled expression data from the cortex of normal adrenal samples.

├── [247M]  scRNA.Normal.adrenal.only.Cortex_harmony.rds
scRNA.Normal.adrenal.only.Cortex_harmony.rds: Harmony algorithm results for the cortex of normal adrenal samples.

├── [ 23K]  scRNA.Normal.adrenal.only.Cortex_signatures.rds
scRNA.Normal.adrenal.only.Cortex_signatures.rds: Gene expression signatures for the cortex of normal adrenal glands.

├── [491M]  scRNA.Normal.adrenal.seurat.rds
scRNA.Normal.adrenal.seurat.rds: Seurat object containing scRNA-seq data for normal adrenal samples.

├── [2.0M]  scRNA.inferCNV.scores.rds
scRNA.inferCNV.scores.rds: inferCNV analysis results, likely for detecting copy number variations in scRNA-seq data.

├── [263M]  scSpatial.ACC.seurat.rds
scSpatial.ACC.seurat.rds: Seurat object containing single-cell spatial transcriptomics data for ACC samples.

└── [161M]  scSpatial.Normal.adrenal.seurat.rds
scSpatial.Normal.adrenal.seurat.rds: Seurat object containing single-cell spatial transcriptomics data for normal adrenal samples.

0 directories, 33 files
~~~

Each file seems to contain specific data subsets or analysis results, crucial for a comprehensive understanding of ACC and normal adrenal tissues at the single-cell level.

# **Citation**

Our paper has been published in [*XXX Journal*](https://chat.openai.com/c/xxxx). For further reference and details, you can access the publication at the provided link.

The raw data supporting the findings of this study can be downloaded from the following repositories:

- **GEO Database**: Access our dataset by visiting [GSEXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXX). This link will take you directly to the dataset's page.
- **Zenodo**: Additional data files are available on Zenodo. Download them at [Zenodo10416598](https://zenodo.org/record/10416598).