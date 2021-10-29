#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 24
#SBATCH --mem=192g
#SBATCH -N 1
#SBATCH -J GATK
# ----------------Load Modules--------------------
module load GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1

chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)

echo "Begin Compiling MAL2 Genomic Database"
geno="/home/MAL/SNP/GATK/MAL2"
ref="/home/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
for int in $(echo ${chr[@]}); do
    gatk --java-options "-Xmx186g" \
	GenomicsDBImport \
	--genomicsdb-workspace-path "${geno}"/genomicsdb \
	--sample-name-map "${geno}"/MAL2_gatk_map \
	--batch-size 0 \
	-L "${int}" \
	--reader-threads $SLURM_NTASKS
    gatk --java-options "-Xmx186g" GenotypeGVCFs \
	-R "${ref}" \
	-V gendb://"${geno}"/genomicsdb \
	-O "${geno}"/MAL2."${int}".combined.raw.vcf.gz
    rm -fr "${geno}"/genomicsdb;
done

echo "Begin Compiling Mo17 Genomic Database"
geno="/home/MAL/SNP/GATK/MAL1"
ref="/home/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa"
for int in $(echo ${chr[@]}); do
    gatk --java-options "-Xmx186g" \
	GenomicsDBImport \
	--genomicsdb-workspace-path "${geno}"/genomicsdb \
	--sample-name-map "${geno}"/MAL1_gatk_map \
	--batch-size 0 \
	-L "${int}" \
	--reader-threads $SLURM_NTASKS
    gatk --java-options "-Xmx186g" GenotypeGVCFs \
        -R "${ref}" \
        -V gendb://"${geno}"/genomicsdb \
        -O "${geno}"/MAL1."${int}".combined.raw.vcf.gz
    rm -fr "${geno}"/genomicsdb;
done
