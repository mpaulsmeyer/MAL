#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 8
#SBATCH --mem=64g
#SBATCH -N 1
#SBATCH -J samtools
#########################################################################################
###################### Variant Calling Pipeline  ########################################
#########################################################################################
module load SAMtools/1.10-IGB-gcc-8.2.0

cwd="/home/MAL" #working directory
ref="/home/ref" #ref directory
log="/home/MAL/SNP/samtools/variant.log" #optional log
depth=6 #read depth per SNP for MAL
max=250 #max read depth for MAL
depth=150 #A3 RNAseq min depth (50 times 3 for each rep)
max=50000 #A3 RNAseq max count

#Need to index the reference genome for the next steps
echo "index B73" | tee -a "${log}"
samtools faidx "${ref}"/Zm-B73-REFERENCE-NAM-5.0.fa
echo "index Mo17" | tee -a "${log}"
samtools faidx "${ref}"/Zm-Mo17-REFERENCE-CAU-1.0.fa
module purge

module load BCFtools/1.12-IGB-gcc-8.2.0
#max is depth saved, B is do not readjust alignments, Ou is a way to pipe commands
#-mv means variants only and rare allele caller
echo "Starting variant pipeline for C Block" | tee -a "${log}"
bcftools mpileup -Ou --skip-indels -d "${max}" --threads $SLURM_NTASKS \
	-f "${ref}"/Zm-B73-REFERENCE-NAM-5.0.fa "${cwd}"/Aligned/MAL2/*.sort.bam | \
bcftools call --threads $SLURM_NTASKS -mv > "${cwd}"/SNP/samtools/MAL2.raw.SNP.vcf

echo "Starting variant pipeline for Mo17 MAL" | tee -a "${log}"
bcftools mpileup -Ou --skip-indels -d "${max}" --threads $SLURM_NTASKS \
	-f "${ref}"/Zm-Mo17-REFERENCE-CAU-1.0.fa "${cwd}"/Aligned/MAL1/*.sort.bam | \
bcftools call --threads $SLURM_NTASKS -mv > "${cwd}"/SNP/samtools/MAL1.raw.SNP.vcf

#none with depth <  $depth, qual < 30 or indels
echo "filtering final VCF" | tee -a "${log}"
bcftools filter -i '%QUAL>=30 & DP>=6 & TYPE="snp"' \
	"${cwd}"/SNP/samtools/MAL2.raw.SNP.vcf > "${cwd}"/SNP/samtools/MAL2.filtered.SNP.vcf
bcftools filter -i '%QUAL>=30 & DP>=6 & TYPE="snp"' \
	"${cwd}"/SNP/samtools/MAL1.raw.SNP.vcf > "${cwd}"/SNP/samtools/MAL1.filtered.SNP.vcf
module purge

#optional step to impute missing data
# module load beagle/5.1-Java-1.8.0_152
# beagle --java-options "-Xmx56g" gt="${cwd}"/SNP/samtools/CBlock.filtered.SNP.vcf \
# out="${cwd}"/SNP/samtools/CBlock.imputed.b1.SNP.vcf gp=true window=10 overlap=5 ne=2
# beagle --java-options "-Xmx56g" gt="${cwd}"/SNP/samtools/CBlock.filtered.SNP.vcf \
# out="${cwd}"/SNP/samtools/CBlock.imputed.b2.SNP.vcf gp=true
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/Mo17.filtered.SNP.vcf \
# out="${cwd}"/SNP/samtools/Mo17.imputed.b1.SNP.vcf gp=true window=10 overlap=5 ne=2
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/Mo17.filtered.SNP.vcf \
# out="${cwd}"/SNP/samtools/Mo17.imputed.b2.SNP.vcf gp=true
