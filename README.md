# MAL

MALVariantCalling.SH is the pipeline I used to align the GBS reads and call variants

"raw" datasets are the output from Samtools filtering

"unsmoothed" datasets are the result of thinning markers from the dataset to be at least 10 kb apart and removing taxa that had less than 100k reads.

The final dataset was subjected to a smoothing function to remove erroneous heterozygous SNPs.
