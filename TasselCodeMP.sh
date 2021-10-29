#!/bin/bash
#Download the latest reference genome from https://download.maizegdb.org/
#Download the latest version of Bowtie from http://bowtie-bio.sourceforge.net/index.shtml
#Download TASSEL 5.0 standalone from https://www.maizegenetics.net/tassel
#Add you sequence to "mySequence"
#Make empty directories for db, maizeTags, logfiles, hdf5, and keys

#Define the directory where run_pipeline.py and all directories are stored
cwd="/home/MAL" 
cd "${cwd}"

"${cwd}"/run_pipeline.pl -Xmx16g -GBSSeqToTagDBPlugin -i "${cwd}"/mySequence/Lane1 \
	-db "${cwd}"/db/maizeTags.db \
	-e PstI -k "${cwd}"/keys/MAL1.txt \
	-kmerLength 64 -minKmerL 20 -batchSize 1 -c 1 \
	-mxKmerNum 100000000 -endPlugin | tee "${cwd}"/logfiles/logfile1.txt
"${cwd}"/run_pipeline.pl -Xmx16g -TagExportToFastqPlugin -db "${cwd}"/db/maizeTags.db \
	-o "${cwd}"/mergedTags/maizeTags.fa.gz \
	-c 1 -endPlugin | tee "${cwd}"/logfiles/logfile2.txt
"${cwd}"/Bowtie/bowtie2 -p 4 --very-sensitive \
	-U "${cwd}"/mergedTags/maizeTags.fa.gz \
	-x "${cwd}"/build/Mo17build \
	-S "${cwd}"/mergedTags/aligned.MAL1.sam  | tee "${cwd}"/logfiles/logfile4.txt
"${cwd}"/run_pipeline.pl -Xmx16g -SAMToGBSdbPlugin -db "${cwd}"/db/maizeTags.db \
	-i "${cwd}"/mergedTags/aligned.MAL1.sam \
	-aProp 0.0 -aLen 0 -endPlugin | tee "${cwd}"/logfiles/logfile5.txt
${cwd}/run_pipeline.pl -Xmx16g -DiscoverySNPCallerPluginV2 -db "${cwd}"/db/maizeTags.db \
	-sC 1 -eC 10 -mnLCov 0.001 -mnMAF 0.001 \
	-endPlugin | tee "${cwd}"/logfiles/logfile6.txt
"${cwd}"/run_pipeline.pl -Xmx16g -SNPQualityProfilerPlugin -db "${cwd}"/db/maizeTags.db \
	-statFile "outputStats.MAL1.txt" 
	-deleteOldData true -endPlugin | tee "${cwd}"/logfiles/logfileMo7.txt
${cwd}/run_pipeline.pl -Xmx16g -ProductionSNPCallerPluginV2 -db "${cwd}"/db/maizeTags.db \
	-e PstI -i "${cwd}"/mySequence/Lane1 -k "${cwd}"/keys/MAL1.txt 
	-o "${cwd}"/hdf5/MAL1.genotypes.h5  
	-batchSize 1 -mnQS 30 -endPlugin | tee "${cwd}"/logfiles/logfile8.txt