

# WES Analysis Pipeline

- Content: Development and implementation of a pipeline for Whole Exome Sequencing (WES) analysis, focusing on genomic alterations in ACC tumors.

~~~java
mkdir -p /mnt/data/userdata/xiangyu/workshop/ACC/mpileup
cd /mnt/data/userdata/xiangyu/workshop/ACC/mpileup

vim sample_file2
TF_ACC
TF_PB
TL_ACC_1
TL_ACC_2
vim sample_file3
TY_ACC
TY_PB
YD_ACC
YD_Adj
vim sample_file4
ZJ_ACC
ZJ_Adj

EXON_BED=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/reference/hg19_Exome-Agilent_V6.bed
GENOME=/mnt/data/userdata/abao/reference/hg19_genome.fa
cat sample_file1 | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
echo mpileupping...
samtools mpileup -d 1000 -q 1 -Q 15 -A -f $GENOME /mnt/data/userdata/abao/project/4_whole_genome_sequencing/WGS_LSZ_230224/${sample}_BQSR.bam -l $EXON_BED > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup
echo varscanning...
java -Xmx100G -jar /mnt/data/userdata/xiangyu/programme/varscan-2.4.2/VarScan.v2.4.2.jar mpileup2cns \
/mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.05 --output-vcf 1 > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.varscan.vcf ;
done

EXON_BED=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/reference/hg19_Exome-Agilent_V6.bed
GENOME=/mnt/data/userdata/abao/reference/hg19_genome.fa
cat sample_file2 | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
echo mpileupping...
samtools mpileup -d 1000 -q 1 -Q 15 -A -f $GENOME /mnt/data/userdata/abao/project/4_whole_genome_sequencing/WGS_LSZ_230224/${sample}_BQSR.bam -l $EXON_BED > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup
echo varscanning...
java -Xmx100G -jar /mnt/data/userdata/xiangyu/programme/varscan-2.4.2/VarScan.v2.4.2.jar mpileup2cns \
/mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.05 --output-vcf 1 > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.varscan.vcf ;
done

EXON_BED=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/reference/hg19_Exome-Agilent_V6.bed
GENOME=/mnt/data/userdata/abao/reference/hg19_genome.fa
cat sample_file3 | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
echo mpileupping...
samtools mpileup -d 1000 -q 1 -Q 15 -A -f $GENOME /mnt/data/userdata/abao/project/4_whole_genome_sequencing/WGS_LSZ_230224/${sample}_BQSR.bam -l $EXON_BED > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup
echo varscanning...
java -Xmx100G -jar /mnt/data/userdata/xiangyu/programme/varscan-2.4.2/VarScan.v2.4.2.jar mpileup2cns \
/mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.05 --output-vcf 1 > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.varscan.vcf ;
done

EXON_BED=/mnt/data/userdata/abao/project/4_whole_genome_sequencing/reference/hg19_Exome-Agilent_V6.bed
GENOME=/mnt/data/userdata/abao/reference/hg19_genome.fa
cat sample_file4 | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
echo mpileupping...
samtools mpileup -d 1000 -q 1 -Q 15 -A -f $GENOME /mnt/data/userdata/abao/project/4_whole_genome_sequencing/WGS_LSZ_230224/${sample}_BQSR.bam -l $EXON_BED > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup
echo varscanning...
java -Xmx100G -jar /mnt/data/userdata/xiangyu/programme/varscan-2.4.2/VarScan.v2.4.2.jar mpileup2cns \
/mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.05 --output-vcf 1 > /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$sample.varscan.vcf ;
done

cd /mnt/data/userdata/xiangyu/workshop/ACC/mpileup
for id in *.varscan.vcf; do
sample=$(basename $id .varscan.vcf)
bcftools view -i 'FILTER="PASS"' $id > $sample.varscan.PASS.vcf ;
done

sed 's/Sample1/JT_ACC_1/g' JT_ACC_1.varscan.PASS.vcf > JT_ACC_1.varscan.PASS.new.vcf
sed 's/Sample1/JT_ACC_2/g' JT_ACC_2.varscan.PASS.vcf > JT_ACC_2.varscan.PASS.new.vcf
sed 's/Sample1/JT_PB_1/g' JT_PB_1.varscan.PASS.vcf > JT_PB_1.varscan.PASS.new.vcf
sed 's/Sample1/LL_Adj/g' LL_Adj.varscan.PASS.vcf > LL_Adj.varscan.PASS.new.vcf
sed 's/Sample1/TFY_ACC/g' TFY_ACC.varscan.PASS.vcf > TFY_ACC.varscan.PASS.new.vcf
sed 's/Sample1/TFY_PB/g' TFY_PB.varscan.PASS.vcf > TFY_PB.varscan.PASS.new.vcf
sed 's/Sample1/TLP_ACC_1/g' TLP_ACC_1.varscan.PASS.vcf > TLP_ACC_1.varscan.PASS.new.vcf
sed 's/Sample1/TLP_ACC_2/g' TLP_ACC_2.varscan.PASS.vcf > TLP_ACC_2.varscan.PASS.new.vcf
sed 's/Sample1/TYR_ACC/g' TYR_ACC.varscan.PASS.vcf > TYR_ACC.varscan.PASS.new.vcf
sed 's/Sample1/TYR_PB/g' TYR_PB.varscan.PASS.vcf > TYR_PB.varscan.PASS.new.vcf
sed 's/Sample1/YDB_ACC/g' YDB_ACC.varscan.PASS.vcf > YDB_ACC.varscan.PASS.new.vcf
sed 's/Sample1/YDB_Adj/g' YDB_Adj.varscan.PASS.vcf > YDB_Adj.varscan.PASS.new.vcf
sed 's/Sample1/ZJM_ACC/g' ZJM_ACC.varscan.PASS.vcf > ZJM_ACC.varscan.PASS.new.vcf
sed 's/Sample1/ZJM_Adj/g' ZJM_Adj.varscan.PASS.vcf > ZJM_Adj.varscan.PASS.new.vcf

cat sample_file1 | while read id ; do
arr=($id)
vcf_name=${arr[0]}
tumor_sam=${arr[0]}
echo $vcf_name
echo $tumor_sam
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.varscan.PASS.new.vcf \
--output-maf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.vep.maf \
--tumor-id $tumor_sam \
--vcf-tumor-id $tumor_sam \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_hg19_index/hg19.fa \
--filter-vcf /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species homo_sapiens/ --ncbi-build GRCh37 ;
done

cat sample_file2 | while read id ; do
arr=($id)
vcf_name=${arr[0]}
tumor_sam=${arr[0]}
echo $vcf_name
echo $tumor_sam
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.varscan.PASS.new.vcf \
--output-maf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.vep.maf \
--tumor-id $tumor_sam \
--vcf-tumor-id $tumor_sam \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_hg19_index/hg19.fa \
--filter-vcf /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species homo_sapiens/ --ncbi-build GRCh37 ;
done

cat sample_file3 | while read id ; do
arr=($id)
vcf_name=${arr[0]}
tumor_sam=${arr[0]}
echo $vcf_name
echo $tumor_sam
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.varscan.PASS.new.vcf \
--output-maf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.vep.maf \
--tumor-id $tumor_sam \
--vcf-tumor-id $tumor_sam \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_hg19_index/hg19.fa \
--filter-vcf /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species homo_sapiens/ --ncbi-build GRCh37 ;
done

cat sample_file4 | while read id ; do
arr=($id)
vcf_name=${arr[0]}
tumor_sam=${arr[0]}
echo $vcf_name
echo $tumor_sam
perl /mnt/data/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.varscan.PASS.new.vcf \
--output-maf /mnt/data/userdata/xiangyu/workshop/ACC/mpileup/$vcf_name.vep.maf \
--tumor-id $tumor_sam \
--vcf-tumor-id $tumor_sam \
--ref-fasta /mnt/data/userdata/xiangyu/programme/genome_index/bwa_index/bwa_hg19_index/hg19.fa \
--filter-vcf /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
--vep-path /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all \
--vep-data /mnt/data/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 --species homo_sapiens/ --ncbi-build GRCh37 ;
done
~~~