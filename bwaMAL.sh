#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 16
#SBATCH --mem=64g
#SBATCH -N 1
#SBATCH -J bwa
# ----------------Load Modules--------------------
module load Stacks/2.54-IGB-gcc-8.2.0
module load BWA/0.7.17-IGB-gcc-8.2.0
module load SAMtools/1.10-IGB-gcc-8.2.0

cwd="/home/MAL" #working directory
ref="/home/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
bwa index -p "${home}"/db/bwa/B73v5 -a bwtsw "${ref}"

#barcodes are in a file with sample\tbarcode on each line
echo "Demultiplex MAL2"
process_radtags -p "${cwd}"/raw/Lane2 -o "${cwd}"/samples -b "${cwd}"/MAL2_barcodes.txt \
	-e pstI --len_limit 50 --retain_header -r -c -q

echo "Demultiplex MAL1"
process_radtags -p "${cwd}"/raw/Lane1 -o "${cwd}"/samples -b "${cwd}"/MAL1_barcodes.txt \
	-e pstI --len_limit 50 --retain_header -r -c -q

echo "Begin C Block Alignment"
for f in $(ls -1 "${cwd}"/Aligned/MAL2/ | cut -d "." -f 1 | sort | uniq); do 
    echo "Aligning sample ${f}"
    RGID=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 13-21)
    RGLB=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 23)
    RGPU=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 9-23)
    bwa mem -t $SLURM_NTASKS \
    	-R "@RG\tID:${RGID}\tPU:${RGPU}\tSM:${f}\tPL:ILLUMINA\tLB:${RGLB}" \
		"${ref}"/db/bwa/B73v5 "${cwd}"/samples/"${f}".fq.gz | \
	samtools view -bh -o "${cwd}"/Aligned/MAL2/"${f}".aligned.bam
    samtools sort -@ $SLURM_NTASKS-1 -o "${cwd}"/Aligned/MAL2/"${f}".aligned.sort.bam \
    	"${cwd}"/Aligned/MAL2/"${f}".aligned.bam
    samtools index -@ $SLURM_NTASKS-1 "${cwd}"/Aligned/MAL2/"${f}".aligned.sort.bam;
done

ref="/home/ref/Zm-Mo17-REFERENCE-CAU-1.0.fa"
bwa index -p "${ref}"/db/bwa/MAL1 -a bwtsw "${ref}"

echo "Begin MAL1 Alignment"
for f in $(ls -1 "${cwd}"/Aligned/MAL1/ | cut -d "." -f 1 | sort | uniq); do 
    echo "Aligning sample ${f}"
    RGID=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 13-21)
    RGLB=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 23)
    RGPU=$(zcat "${cwd}"/samples/"${f}".fq.gz | head -n 1 | cut -c 9-23)
    bwa mem -t $SLURM_NTASKS \
    	-R "@RG\tID:${RGID}\tPU:${RGPU}\tSM:${f}\tPL:ILLUMINA\tLB:${RGLB}" \
		"${ref}"/db/bwa/MAL1 "${cwd}"/samples/"${f}".fq.gz | \
	samtools view -bh -o "${cwd}"/Aligned/MAL1/"${f}".aligned.bam
    samtools sort -@ $SLURM_NTASKS-1 -o "${cwd}"/Aligned/MAL1/"${f}".aligned.sort.bam \
    	"${cwd}"/Aligned/MAL1/"${f}".aligned.bam
    samtools index -@ $SLURM_NTASKS-1 "${cwd}"/Aligned/MAL1/"${f}".aligned.sort.bam;
done
