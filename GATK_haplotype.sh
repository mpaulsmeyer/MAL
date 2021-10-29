#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=192g
#SBATCH -N 1
#SBATCH -J GATK
# ----------------Load Modules--------------------  
module load GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1
module load BCFtools/1.9-IGB-gcc-4.9.4

cwd="/home/MAL" #working directory
chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)

echo "Begin GATK Pipeline for MAL2"
ref="/home/ref/Zm-B73-REFERENCE-NAM-5.0.fa" #B73 ref v5 genome
geno="$cwd/SNP/GATK/MAL2"

for f in $(ls -1 "${cwd}"/Aligned/MAL2/ | cut -d "." -f 1 | sort | uniq); do 
    timestamp=$(date +%H:%M:%S); echo "Calling SNPs For $f at $timestamp"
    for int in $(echo ${chr[@]}); do
        gatk --java-options "-Xmx16g" HaplotypeCaller \
            -R "${ref}" \
            -I "${cwd}"/Aligned/MAL2/"${f}".aligned.sort.bam \
            -O "${geno}"/"${int}".c.vcf.gz \
            -ERC GVCF \
            -L "${int}" \
            #neat trick to parallelize single core processes. Saves hours
            --native-pair-hmm-threads 6 & done; #& done runs on all cores
        wait
	echo "${f}	${f}.raw.g.vcf.gz" >> "${geno}"/MAL2_gatk_map.txt
    bcftools concat -o "${geno}"/"${f}".raw.g.vcf "${geno}"/chr1.c.vcf.gz \
    	"${geno}"/chr2.c.vcf.gz "${geno}"/chr3.c.vcf.gz "${geno}"/chr4.c.vcf.gz \
    	"${geno}"/chr5.c.vcf.gz "${geno}"/chr6.c.vcf.gz "${geno}"/chr7.c.vcf.gz \
    	"${geno}"/chr8.c.vcf.gz "${geno}"/chr9.c.vcf.gz "${geno}"/chr10.c.vcf.gz
    rm "${geno}"/*.c.gz "${geno}"/*.tbi
    gzip "${geno}"/"${f}".raw.g.vcf #optional
    timestamp=$(date +%H:%M:%S); echo "Ended Calling SNPs For $f at $timestamp";
done

echo "Begin GATK Pipeline for MAL1"
ref="/home/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa" #Mo17 genome
geno="$cwd/SNP/GATK/MAL1"

for f in $(ls -1 "${cwd}"/Aligned/MAL1/ | cut -d "." -f 1 | sort | uniq); do
    timestamp=$(date +%H:%M:%S); echo "Calling SNPs For $f at $timestamp"
    for int in $(echo ${chr[@]}); do
        gatk --java-options "-Xmx16g" HaplotypeCaller \
            -R "${ref}" \
	    	-I "${cwd}"/Aligned/MAL1/"${f}".aligned.sort.bam \
            -O "${geno}"/"${int}".c.vcf.gz \
            -ERC GVCF \
            -L "${int}" \
            --native-pair-hmm-threads 6 & done;
	wait
	echo "${f}	${f}.raw.g.vcf.gz" >> "${geno}"/MAL1_gatk_map.txt
    bcftools concat -o "${geno}"/"${f}".raw.g.vcf "${geno}"/chr1.c.vcf.gz \
    	"${geno}"/chr2.c.vcf.gz "${geno}"/chr3.c.vcf.gz "${geno}"/chr4.c.vcf.gz \
    	"${geno}"/chr5.c.vcf.gz "${geno}"/chr6.c.vcf.gz "${geno}"/chr7.c.vcf.gz \
    	"${geno}"/chr8.c.vcf.gz "${geno}"/chr9.c.vcf.gz "${geno}"/chr10.c.vcf.gz
    rm "${geno}"/*.c.gz "${geno}"/*.tbi
    gzip "${geno}"/"${f}".raw.g.vcf #optional
    timestamp=$(date +%H:%M:%S); echo "Ended Calling SNPs For $f at $timestamp";
done

#You can split this next part into a different job and use 24 CPUs to be more efficient
echo "Begin Compiling MAL2 Genomic Database"
geno="$cwd/SNP/GATK/MAL2"
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
geno="$cwd/SNP/GATK/MAL1"
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
