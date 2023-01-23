#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 16
#SBATCH --mem=64g
#SBATCH -N 1
#SBATCH --mail-user=paulsme2@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J MP_GATK
# ----------------Load Modules--------------------
module load Stacks/2.54-IGB-gcc-8.2.0
module load BWA/0.7.17-IGB-gcc-8.2.0
module load SAMtools/1.10-IGB-gcc-8.2.0

cwd="/home/n-z/paulsme2/MAL"
home="/home/n-z/paulsme2"
ref="/home/n-z/paulsme2/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

#cat KeyFile2019MAL.txt | cut -f 2- | grep "^3" | cut -f 2-3 | grep -v "D" > Mo17_barcodes.txt
#cat KeyFile2019MAL.txt | cut -f 2- | grep "^4" | cut -f 2-3 | grep -v "D" > CBlock_barcodes.txt
#cat CBlock_barcodes.txt | grep -v "C083" | grep -v "C093" | grep -v "C095" | grep -v "C096" > CBlock_barcodes_clean.txt

echo "Demultiplex C Block"
process_radtags -p "${cwd}"/raw/Lane2 -o "${cwd}"/samples -b "${cwd}"/CBlock_barcodes_clean.txt \
-e pstI --len_limit 50 --retain_header -r -c -q

echo "Demultiplex Mo17"
process_radtags -p "${cwd}"/raw/Lane1 -o "${cwd}"/samples -b "${cwd}"/Mo17_barcodes.txt \
-e pstI --len_limit 50 --retain_header -r -c -q

#created index already
#bwa index -p "${home}"/db/bwa/B73v5 -a bwtsw "${ref}"

echo "Begin C Block Alignment"
for f in $(ls -1 "${cwd}"/samples/C* | sed -E 's/\/home\/n-z\/paulsme2\/MAL\/samples\/(C[0-9]{3}).fq.gz/\1/' | cut -d "." -f 1 | sort | uniq); do 
	echo "Aligning sample ${f}"
	RGID=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 13-21)
	RGLB=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 23)
	RGPU=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 9:23)
	bwa mem -t $SLURM_NTASKS -R "@RG\tID:${RGID}\tPU:${RGPU}\tSM:${f}\tPL:ILLUMINA\tLB:${RGLB}" \
	"${home}"/db/bwa/B73v5 "${cwd}"/"${f}".fq.gz | \
	samtools view -bh -o "${cwd}"/Aligned/CBlock/"${f}".aligned.bam
	samtools index -@ $SLURM_NTASKS-1 "${cwd}"/Aligned/CBlock/"${f}".aligned.bam
	samtools sort -@ $SLURM_NTASKS-1 -o "${cwd}"/Aligned/CBlock/"${f}".aligned.bam "${cwd}"/Aligned/CBlock/"${f}".aligned.bam;
done

ref="/home/n-z/paulsme2/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa"
#created index already
#bwa index -p "${home}"/db/bwa/Mo17 -a bwtsw "${ref}"

echo "Begin Mo17 Alignment"
for f in $(ls -1 "${cwd}"/samples/2* | sed -E 's/\/home\/n-z\/paulsme2\/MAL\/samples\/(2[0-9]{2}-[0-9]+).fq.gz/\1/' | cut -d "." -f 1 | sort | uniq)); do 
	echo "Aligning sample ${f}"
	RGID=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 13-21)
	RGLB=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 23)
	RGPU=$(zcat "${cwd}"/"${f}".fq.gz | head -n 1 | cut -c 9:23)
	bwa mem -t $SLURM_NTASKS -R "@RG\tID:${RGID}\tPU:${RGPU}\tSM:${f}\tPL:ILLUMINA\tLB:${RGLB}" \
	"${home}"/db/bwa/Mo17 "${cwd}"/"${f}".fq.gz | \
	samtools view -bh -o "${cwd}"/Aligned/Mo17/"${f}".aligned.bam
	samtools index -@ $SLURM_NTASKS-1 "${cwd}"/Aligned/Mo17/"${f}".aligned.bam
	samtools sort -@ $SLURM_NTASKS-1 -o "${cwd}"/Aligned/Mo17/"${f}".aligned.bam "${cwd}"/Aligned/Mo17/"${f}".aligned.bam;
done

########################################################################################################################
#!/bin/bash                                                                                                                         
# ----------------SLURM Parameters----------------                                                                                  
#SBATCH -p normal                                                                                                                   
#SBATCH -n 64                                                                                                                       
#SBATCH --mem=172g                                                                                                                  
#SBATCH -N 1                                                                                                                        
#SBATCH --mail-user=paulsme2@illinois.edu                                                                                           
#SBATCH --mail-type=ALL                                                                                                             
#SBATCH -J MP_GATK                                                                                                                  
# ----------------Load Modules--------------------  
#6 min per sample
module load GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1
module load BCFtools/1.9-IGB-gcc-4.9.4

mkdir /scratch/paulsme2
cwd="/home/n-z/paulsme2/MAL"
chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)

echo "Begin GATK Pipeline for C Block"
ref="/home/n-z/paulsme2/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
geno="/home/n-z/paulsme2/MAL/SNP/GATK/CBlock"

for f in $(ls -1 "${cwd}"/samples/C* | sed -E 's/\/home\/n-z\/paulsme2\/MAL\/samples\/(C[0-9]{3}).fq.gz/\1/' | cut -d "." -f 1 | sort | uniq); do 
    timestamp=$(date +%H:%M:%S); echo "Calling SNPs For $f at $timestamp"
    for int in $(echo ${chr[@]}); do
        gatk --java-options "-Xmx16g" HaplotypeCaller \
            -R "${ref}" \
            -I "${cwd}"/Aligned/CBlock/"${f}".aligned.sort.bam \
            -O "${geno}"/"${int}".c.vcf \
            -ERC GVCF \
            -L "${int}" \
            --native-pair-hmm-threads 6 & done;
        wait
    echo "${f}	MAL/SNP/GATK/CBlock/${f}.raw.g.vcf" >> "${geno}"/CBlock_gatk_map
    bcftools concat -o "${geno}"/"${f}".raw.g.vcf "${geno}"/chr1.c.vcf.gz "${geno}"/chr2.c.vcf \
         "${geno}"/chr3.c.vcf.gz "${geno}"/chr4.c.vcf "${geno}"/chr5.c.vcf "${geno}"/chr6.c.vcf \
         "${geno}"/chr7.c.vcf.gz "${geno}"/chr8.c.vcf "${geno}"/chr9.c.vcf "${geno}"/chr10.c.vcf
    rm "${geno}"/*.c "${geno}"/*.tbi
    gatk IndexFeatureFile -F "${geno}"/"${f}".raw.g.vcf
    timestamp=$(date +%H:%M:%S); echo "Ended Calling SNPs For $f at $timestamp";
done

echo "Begin GATK Pipeline for Mo17"
ref="/home/n-z/paulsme2/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa"
geno="/home/n-z/paulsme2/MAL/SNP/GATK/Mo17"

for f in $(ls -1 "${cwd}"/Aligned/Mo17/ | sed -E 's/\/home\/n-z\/paulsme2\/MAL\/samples\/(2[0-9]{2}-[0-9]+).fq.gz/\1/' | cut -d "." -f 1 | sort | uniq); do
    timestamp=$(date +%H:%M:%S); echo "Calling SNPs For $f at $timestamp"
    for int in $(echo ${chr[@]}); do
        gatk --java-options "-Xmx16g" HaplotypeCaller \
            -R "${ref}" \
            -I "${cwd}"/Aligned/Mo17/"${f}".aligned.sort.bam \
            -O "${geno}"/"${int}".c.vcf \
            -ERC GVCF \
            -L "${int}" \
            --native-pair-hmm-threads 6 & done;
        wait
    echo "${f}	MAL/SNP/GATK/Mo17/${f}.raw.g.vcf" >> "${geno}"/Mo17_gatk_map
    bcftools concat -o "${geno}"/"${f}".raw.g.vcf "${geno}"/chr1.c.vcf "${geno}"/chr2.c.vcf \
         "${geno}"/chr3.c.vcf "${geno}"/chr4.c.vcf "${geno}"/chr5.c.vcf "${geno}"/chr6.c.vcf \
         "${geno}"/chr7.c.vcf "${geno}"/chr8.c.vcf "${geno}"/chr9.c.vcf "${geno}"/chr10.c.vcf
    rm "${geno}"/*.c "${geno}"/*.tbi
    gatk IndexFeatureFile -F "${geno}"/"${f}".raw.g.vcf
    timestamp=$(date +%H:%M:%S); echo "Ended Calling SNPs For $f at $timestamp";
done

#grep "Calling SNPs" slurm-7250975.out
#ls -l MAL/SNP/GATK/CBlock/
#less slurm-7250975.out
########################################################################################################################
#cat MAL/SNP/GATK/CBlock/CBlock_gatk_map.txt | sed -E 's/(C[0-9]{3})	C[0-9]{3}.raw.g.vcf.gz/\1	/home/n-z/paulsme2/MAL/SNP/GATK/CBlock/\1.raw.g.vcf.gz/' > CBlock_gatk_map2.txt
#cat MAL/SNP/GATK/Mo17/Mo17_gatk_map.txt | sed -E 's/(2[0-9]{3}-[0-9]+)	2[0-9]{3}-[0-9]+\.raw.g.vcf.gz/\1	\/home\/n-z\/paulsme2\/MAL\/SNP\/GATK\/Mo17\/\1.raw.g.vcf.gz/' > Mo17_gatk_map2.txt

#renaming gzipped files to regular files
for f in $(cat MAL/SNP/GATK/CBlock/CBlock_gatk_map | cut -f 2 | sed -E "s/MAL\/SNP\/GATK\/CBlock\/(C[0-9]{3}).raw.g.vcf.gz/MAL\/SNP\/GATK\/CBlock\/\1.raw.g.vcf/"); do mv "${f}".gz $f; done
for f in $(cat MAL/SNP/GATK/Mo17/Mo17_gatk_map | cut -f 2 | sed -E "s/MAL\/SNP\/GATK\/Mo17\/(2[0-9-]+).raw.g.vcf.gz/MAL\/SNP\/GATK\/Mo17\/\1.raw.g.vcf/"); do mv "${f}".gz $f; done

for f in $(ls -1 MAL/SNP/GATK/CBlock/*.g.vcf); do echo $f; gatk IndexFeatureFile -F $f; done
for f in $(cat MAL/SNP/GATK/Mo17/Mo17_gatk_map | cut -f 2); do echo $f; gatk IndexFeatureFile -F $f; done

for f in $(ls -1 *.combined.raw.vcf.gz | sed -E "s/Mo17.chr([0-9]+).combined.raw.vcf.gz/Mo17.chr\1.combined.raw.vcf/"); do echo $f; mv "${f}".gz $f; done
rm *.tbi

########################################################################################################################
#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 24
#SBATCH --mem=192g
#SBATCH -N 1
#SBATCH --mail-user=paulsme2@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J MP_GATK
# ----------------Load Modules--------------------
#10 min per chromsome, not sure if 24 threads was necessary
module load GATK/4.0.9.0-IGB-gcc-4.9.4-Java-1.8.0_152-Python-3.6.1

mkdir /scratach/paulsme2
chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)

echo "Begin Compiling CBlock Genomic Database"
geno="/home/n-z/paulsme2/MAL/SNP/GATK/CBlock"
ref="/home/n-z/paulsme2/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
for int in $(echo ${chr[@]}); do
    gatk --java-options "-Xmx186g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path "${geno}"/genomicsdb \
        --sample-name-map "${geno}"/CBlock_gatk_map \
        --batch-size 0 \
        -L "${int}" \
        --reader-threads $SLURM_NTASKS
    gatk --java-options "-Xmx186g" GenotypeGVCFs \
        -R "${ref}" \
        -V gendb://"${geno}"/genomicsdb \
        -O "${geno}"/CBlock."${int}".combined.raw.vcf
    rm -fr "${geno}"/genomicsdb;
done

echo "Begin Compiling Mo17 Genomic Database"
geno="/home/n-z/paulsme2/MAL/SNP/GATK/Mo17"
ref="/home/n-z/paulsme2/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa"
for int in $(echo ${chr[@]}); do
    gatk --java-options "-Xmx186g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path "${geno}"/genomicsdb \
        --sample-name-map "${geno}"/Mo17_gatk_map \
        --batch-size 0 \
        -L "${int}" \
        --reader-threads $SLURM_NTASKS
    gatk --java-options "-Xmx186g" GenotypeGVCFs \
        -R "${ref}" \
        -V gendb://"${geno}"/genomicsdb \
        -O "${geno}"/Mo17."${int}".combined.raw.vcf
    rm -fr "${geno}"/genomicsdb;
done

module purge
module load BCFtools/1.12-IGB-gcc-8.2.0
geno="/home/n-z/paulsme2/MAL/SNP/GATK/CBlock"
for f in $(ls -1 "${geno}"/*raw.vcf); do echo $f; bcftools index $f; done
bcftools concat -O z -o "${geno}"/CBlock.combined.raw.vcf.gz "${geno}"/CBlock.chr1.combined.raw.vcf "${geno}"/CBlock.chr2.combined.raw.vcf "${geno}"/CBlock.chr3.combined.raw.vcf "${geno}"/CBlock.chr4.combined.raw.vcf \
"${geno}"/CBlock.chr5.combined.raw.vcf "${geno}"/CBlock.chr6.combined.raw.vcf "${geno}"/CBlock.chr7.combined.raw.vcf "${geno}"/CBlock.chr8.combined.raw.vcf "${geno}"/CBlock.chr9.combined.raw.vcf "${geno}"/CBlock.chr10.combined.raw.vcf
bcftools filter -i '%QUAL>=30 & INFO/DP>=6 & TYPE="snp"' -O z -o "${geno}"/CBlock.combined.filtered.vcf.gz "${geno}"/CBlock.combined.raw.sort.vcf.gz
bcftools index "${geno}"/CBlock.combined.filtered.vcf.gz
bcftools sort -m 64G --temp-dir /scratch/paulsme2 -O z -o "${geno}"/CBlock.combined.filtered.sort.vcf.gz "${geno}"/CBlock.combined.filtered.vcf.gz

geno="/home/n-z/paulsme2/MAL/SNP/GATK/Mo17"
for f in $(ls -1 "${geno}"/*raw.vcf); do echo $f; bcftools index $f; done
bcftools concat -O z -o "${geno}"/Mo17.combined.raw.vcf.gz "${geno}"/Mo17.chr1.combined.raw.vcf "${geno}"/Mo17.chr2.combined.raw.vcf "${geno}"/Mo17.chr3.combined.raw.vcf "${geno}"/Mo17.chr4.combined.raw.vcf \
"${geno}"/Mo17.chr5.combined.raw.vcf "${geno}"/Mo17.chr6.combined.raw.vcf "${geno}"/Mo17.chr7.combined.raw.vcf "${geno}"/Mo17.chr8.combined.raw.vcf "${geno}"/Mo17.chr9.combined.raw.vcf "${geno}"/Mo17.chr10.combined.raw.vcf
bcftools filter -i '%QUAL>=30 & INFO/DP>=6 & TYPE="snp"' -O z -o "${geno}"/Mo17.combined.filtered.vcf.gz "${geno}"/Mo17.combined.raw.sort.vcf.gz
bcftools index "${geno}"/Mo17.combined.filtered.vcf.gz
bcftools sort -m 64G --temp-dir /scratch/paulsme2 -O z -o "${geno}"/Mo17.combined.filtered.sort.vcf.gz "${geno}"/Mo17.combined.filtered.vcf.gz

########################################################################################################################
module purge
module load SAMtools/1.10-IGB-gcc-8.2.0

cwd="/home/n-z/paulsme2/MAL"
ref="/home/n-z/paulsme2/ref"
log="/home/n-z/paulsme2/MAL/SNP/samtools/variant.log"
depth=6
max=250

echo "index B73" | tee -a "${log}"
samtools faidx "${ref}"/Zm-B73-REFERENCE-NAM-5.0.fa
echo "index Mo17" | tee -a "${log}"
samtools faidx "${ref}"/Zm-Mo17-REFERENCE-CAU-1.0.fa

module purge

module load BCFtools/1.12-IGB-gcc-8.2.0
#max is depth saved, B is do not readjust alignments, Ou is a way to pipe commands
#-mv means variants only and rare allele caller
echo "Starting variant pipeline for C Block" | tee -a "${log}"
bcftools mpileup -Ou --skip-indels -d "${max}" --threads $SLURM_NTASKS -f "${ref}"/Zm-B73-REFERENCE-NAM-5.0.fa "${cwd}"/Aligned/CBlock/*sort.bam | \
bcftools call --threads $SLURM_NTASKS -mv > "${cwd}"/SNP/samtools/CBlock.raw.SNP.vcf

echo "Starting variant pipeline for Mo17 MAL" | tee -a "${log}"
bcftools mpileup -Ou --skip-indels -d "${max}" --threads $SLURM_NTASKS -f "${ref}"/Zm-Mo17-REFERENCE-CAU-1.0.fa "${cwd}"/Aligned/Mo17/*sort.bam | \
bcftools call --threads $SLURM_NTASKS -mv > "${cwd}"/SNP/samtools/Mo17.raw.SNP.vcf

#none with depth <  $depth, qual < 30 or indels
echo "filtering final VCF" | tee -a "${log}"
bcftools filter -i '%QUAL>=30 & DP>=6 & TYPE="snp"' "${cwd}"/SNP/samtools/CBlock.raw.SNP.vcf > "${cwd}"/SNP/samtools/CBlock.filtered.SNP.vcf
bcftools filter -i '%QUAL>=30 & DP>=6 & TYPE="snp"' "${cwd}"/SNP/samtools/Mo17.raw.SNP.vcf > "${cwd}"/SNP/samtools/Mo17.filtered.SNP.vcf
module purge

# module load beagle/5.1-Java-1.8.0_152
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/CBlock.filtered.SNP.vcf out="${cwd}"/SNP/samtools/CBlock.imputed.b1.SNP.vcf gp=true window=10 overlap=5 ne=2
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/CBlock.filtered.SNP.vcf out="${cwd}"/SNP/samtools/CBlock.imputed.b2.SNP.vcf gp=true
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/Mo17.filtered.SNP.vcf out="${cwd}"/SNP/samtools/Mo17.imputed.b1.SNP.vcf gp=true window=10 overlap=5 ne=2
# beagle -Xmx60g gt="${cwd}"/SNP/samtools/Mo17.filtered.SNP.vcf out="${cwd}"/SNP/samtools/Mo17.imputed.b2.SNP.vcf gp=true
