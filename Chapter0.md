# code of data pre-analysis

## 1. **Single-cell and Single-nucleus Transcriptome Data Analysis**

For our analysis, we utilized the cellranger count model (version 5.0.1) to align sequencing files generated by Single Cell 5' PE chemistry with the human reference genome GRCh38-2020-A. For single-nucleus samples processed by Single Cell 3' v3 chemistry, the `--include-introns` parameter was added during alignment. Loupe files were generated using the reanalyze model of cellranger (version 5.0.1).

Quality control, clustering, dimension reduction, cell type annotation, and visualization were performed using Seurat (version 3.2.3). Cells with fewer than 200 detectable genes, more than 8000 detectable genes, or with mitochondrial gene content exceeding 15% were considered poor-quality and removed. Similarly, genes detected in fewer than three cells were also excluded from further analysis. Following these quality control measures, a total of 154,428 cells were harvested across various samples, with specific counts for each normal and ACC sample.

For feature selection, the `FindVariableFeatures` function was used to identify 4,000 highly variable genes using the vst model. This step was essential for data scaling and subsequent clustering and dimension reduction. Twenty-five principal component analysis (PCA) factors were employed for embedding calculations, and twenty for neighbors’ estimation. Dimensionality reduction was performed using openTSNE, t-SNE, UMAP, and force atlas2 (fa2), following previously reported methods. Major populations and subpopulations were identified using classical signatures, and top markers were calculated using the `FindAllMarkers` function with a Wilcoxon test.

To analyze both normal adrenal cortex and tumor subpopulations, we integrated the harmony algorithm (version 0.1.0) to mitigate batch effects arising from different platforms and library construction chemistries. This was crucial for identifying common molecular features in ACC. Twenty harmony factors were used for clustering and dimensionality reduction.

The code for cellranger processing is as follows:

~~~shel
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
vim scRNA_v1.sh

mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC1_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC1_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC1_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC1_RNA \
--sample=ACC1_RNA \
--nosecondary
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC2_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC2_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC2_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC2_RNA \
--sample=ACC2_RNA \
--nosecondary
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC3_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC3_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC3_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC3_RNA \
--sample=ACC3_RNA \
--nosecondary
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC4_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC4_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC4_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC4_RNA \
--sample=ACC4_RNA \
--nosecondary
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal1_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal1_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=Normal1_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA_Normal1/N2/ \
--sample=Normal1_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal2_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal2_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=Normal2_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA_Normal1/X0805/ \
--sample=Normal2_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal3_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/Normal3_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=Normal3_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA_Normal1/ACC_scRNA_Normal_v2/ \
--sample=OES215872150-L23-01 \
--nosecondary --include-introns

cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
vim scRNA_v2.sh
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC5_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC5_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC5_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/HT2021-10647-1/S0_3_0_5/ACC5_RNA \
--sample=ACC5_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC6_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC6_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC6_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/HT2021-10647-1/S1_6_1_10/ACC6_RNA \
--sample=ACC6_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC7_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC7_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC7_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/S730ai/ACC7_RNA \
--sample=ACC7_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC8_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC8_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC8_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/S730aipang/ACC8_RNA \
--sample=ACC8_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC9_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC9_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC9_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/S850_1/ACC9_RNA \
--sample=ACC9_RNA \
--nosecondary --include-introns
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC10_RNA
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC10_RNA
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger count \
--id=ACC10_RNA \
--localcores 20 \
--transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/scNuclei/S851_2_851_4/ACC10_RNA \
--sample=ACC10_RNA \
--nosecondary --include-introns

cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
bash scRNA_v1.sh
bash scRNA_v2.sh
~~~



## **2. Single-cell immune profiling processing and analysis.** 

The vdj model of cellranger (v5.0.1) was used to align the sequencing files generated by Single Cell V(D)J with the human reference genome of vdj_GRCh38_alts_ensembl-5.0.0. The loupe files were generated by cellranger (v5.0.1) with reanalyze function. The filtered contig annotations and clonotypes files were used to generate the TCR α, TCR β-chain and CDR3 amino acid sequence. T cells without TCR annotation, B cells and NK cells were removed for clonal proportion analysis in lymphocyte lineage. T cells with the same TCR α, TCR β-chain and clonotype were regarded as clonal cells. TCR clonotype compromising more than 100 T cells were selected as the top enriched clonotype. The T cells with high proportion and low-diversity TCR clonotype were regarded as specific anti-tumor T cells. 

The code for cellranger processing is as follows:

~~~shell
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
vim scTCR_v1.sh

mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC1_TCR
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC1_TCR
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger vdj \
--id=ACC1_TCR \
--reference=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC1_TCR \
--sample=ACC1_TCR \
--localcores=20
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC2_TCR
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC2_TCR
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger vdj \
--id=ACC2_TCR \
--reference=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC2_TCR \
--sample=ACC2_TCR \
--localcores=20
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC3_TCR
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC3_TCR
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger vdj \
--id=ACC3_TCR \
--reference=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC3_TCR \
--sample=ACC3_TCR \
--localcores=20
mkdir /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC4_TCR
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project/ACC4_TCR
/mnt/data/user_data/xiangyu/programme/cellranger-5.0.1/cellranger vdj \
--id=ACC4_TCR \
--reference=/mnt/data/user_data/xiangyu/programme/genome_index/cellranger_ref_10x/human/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
--fastqs=/mnt/data/sequencedata/scRNA/ACC_scRNA/ACC4_TCR \
--sample=ACC4_TCR \
--localcores=20

cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project
bash scTCR_v1.sh
~~~



## 3. **Alignment and annotation in spatial transcriptome.** 

The count model of spaceranger (v1.2.2) was used to align the sequencing files generated by Spatial 3' (v1) with the human reference genome of GRCh38-2020-A using parameter of --reorient-images. The area A1, B1, C1 and D1 in silde of V11D06-026 were used to construct the library of spatial transcriptome sample. 
 The Seurat (v3.2.3), Matrix (v1.4.1) and spacexr (v2.0.1) were implemented into clustering, cell type annotation and visualization. The FindVariableFeatures was implemented to identify 4,000 high-variable genes with vst model for data scaling. The 30 principal component analysis factors were used to embedding calculation and the 20 principal component analysis factors were used to neighbors’ estimation. The Reference function was used to generate the reference files by inputting the counts of scRNA-seq data, cell annotation information and number of UMI with parameter of n_max_cells= 20000. The RCTD algorithm implemented in spacexr (v2.0.1) was used to deconvolute the potential cell type in each spot of spatial transcriptome data by inputting with spatial counts and parameter of doublet_mode = doublet. And the normalized weights of RCTD results were stored into meta.data in seurat object. The second cell type were regraded as the final annotation results in each spot. 

~~~r
vim scSpatial_v1.sh
export PATH=/mnt/data/user_data/xiangyu/programme/spaceranger-1.2.2:$PATH
mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC11_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC11_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/CP2021062400128/H101SC22033524/RSSQ00504/X101SC22033524-Z01/X101SC22033524-Z01-J062/01.RawData/唐云瑞-E1 \
    --sample=TYR \
    --id=ACC11_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/picture/KJ-1_01.tif \
    --slide=V11D06-026 \
    --area=A1 --localcores 20

mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC1_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC1_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/CP2021062400128/H101SC22033524/RSSQ00504/X101SC22033524-Z01/X101SC22033524-Z01-J062/01.RawData/汤丽萍-E2 \
    --sample=YLP \
    --id=ACC1_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/picture/KJ-1_02.tif \
    --slide=V11D06-026 \
    --area=C1 --localcores 20

mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC12_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC12_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/CP2021062400128/H101SC22033524/RSSQ00504/X101SC22033524-Z01/X101SC22033524-Z01-J062/01.RawData/田凤英-E3 \
    --sample=TFY \
    --id=ACC12_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/picture/KJ-1_03.tif \
    --slide=V11D06-026 \
    --area=B1 --localcores 20

mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC13_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/ACC13_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/CP2021062400128/H101SC22033524/RSSQ00504/X101SC22033524-Z01/X101SC22033524-Z01-J062/01.RawData/蒋婷-E4 \
    --sample=JT \
    --id=ACC13_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/picture/KJ-1_04.tif \
    --slide=V11D06-026 \
    --area=D1 --localcores 20

export PATH=/mnt/data/user_data/xiangyu/programme/spaceranger-1.2.2:$PATH
mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal1_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal1_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/N01_1/ \
    --sample=OES21582590B-01 \
    --id=Normal1_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/tif/N01_1.tif \
    --slide=V11D08-325 \
    --area=A1 --localcores 20

mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal2_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal2_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/S89_3/ \
    --sample=OES21582790B-01 \
    --id=Normal2_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/tif/S89_3.tif \
    --slide=V11D08-325 \
    --area=C1 --localcores 20

mkdir -p /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal3_spatial
cd /mnt/data/user_data/xiangyu/workshop/scRNA/ACC_project_spatial/Normal3_spatial
spaceranger count  --reorient-images \
    --transcriptome=/mnt/data/user_data/xiangyu/programme/genome_index/spatiaranger_10x/human/refdata-gex-GRCh38-2020-A \
    --fastqs=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/X0805_2/ \
    --sample=OES21582690B-01 \
    --id=Normal3_spatial \
    --image=/mnt/data/sequencedata/scSpatial_RNAseq/ACC_project/ACC_normal_scSpatial/DOE20226654-b1/tif/X0805_2.tif \
    --slide=V11D08-325 \
    --area=B1 --localcores 20
~~~




