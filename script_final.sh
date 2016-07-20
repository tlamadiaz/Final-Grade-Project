#!/bin/bash

set -euo pipefail

coding_regions_ucsc="coding_regions.bed.gz"

# print header
echo -ne "CHR\tSIZE_chromosome\tSIZE_coding_regions\t" > results.txt
echo -ne "SNVs\tts_tv_SNVs\t" >> results.txt
echo -ne "SNPs\tts_tv_SNPs\t" >> results.txt
echo -ne "cSNVs\tts_tv_cSNVs\t" >> results.txt
echo -ne "SNVs_frequency\tcSNVs_frequency\tSNPs_frequency\tcSNPs_frequency\n" >> results.txt


#loop for the analysis of 23 chromosomes
for num in {1..22} X Y; do

	# file names
	whole_chromosome_vcf="ALL.chr${num}.*.genotypes.vcf.gz"
	coding_regions_bed="chr${num}_coding_regions.bed"
	SNV_stats="chr${num}_SNVs_stats.txt"
	SNPs_stats="chr${num}_SNPs_stats.txt"
	SNVs_coding_regions_stats="chr${num}_SNVs_coding_regions_stats.txt"
	SNPs_coding_regions_stats="chr${num}_SNPs_coding_regions_stats.txt"

	# get regions of interest
	gunzip -c $coding_regions_ucsc \
	| grep -P "^chr${num}\t" \
	| sed 's/^chr//' \
	| cut -f1-3 \
	> $coding_regions_bed

	# chromosome number
	echo -n "chromosome number: "
	CHR=$(echo -n "${num}")
	echo $CHR

	# chromosome size
	echo -n "chr${num} size "
	SIZE_chromosome=$(bcftools view -h $whole_chromosome_vcf | grep "^##contig=<ID=${num}," | sed 's/.*length=\([0-9]*\).*/\1/')
	echo $SIZE_chromosome

	# coding regions size
	echo -n "chr${num} exons size: "
	SIZE_coding_regions=$(gunzip -c coding_regions.bed.gz | grep -P "^chr${num}\t" | awk '{s+= $3 - $2} END {print s}')
	echo $SIZE_coding_regions

	# statistics whole chromosome
	if [ ! -f $SNV_stats ]; then
		bcftools stats $whole_chromosome_vcf > $SNV_stats
	fi
	echo -n "chr${num} SNPs: "
	SNVs=$(cat $SNV_stats| awk '{if ($5=="SNPs:"){print $NF}}')
	echo $SNVs
	echo -n "chr${num} ts/tv: "
	ts_tv_SNVs=$(cat $SNV_stats | grep -A1 "ts/tv" | tail -1 | cut -f5)
	echo $ts_tv_SNVs

	# SNPs' statistics of whole chromosome
	if [ ! -f $SNPs_stats ]; then
		bcftools view -i 'AF>=0.01' $whole_chromosome_vcf | bcftools stats > $SNPs_stats;
	fi
	echo -n "chr${num} SNPs: "
	SNPs=$(cat $SNPs_stats | awk '{if ($5=="SNPs:"){print $NF}}')
	echo $SNPs
	echo -n "chr${num} ts/tv SNPs: "
	ts_tv_SNPs=$(cat $SNPs_stats | grep -A1 "ts/tv" | tail -1 | cut -f5)
	echo $ts_tv_SNPs

	# SNVs statistics of coding regions
	if [ ! -f $SNVs_coding_regions_stats ]; then
		bcftools stats -R $coding_regions_bed $whole_chromosome_vcf > $SNVs_coding_regions_stats
	fi
	echo -n "chr${num} SNVs in coding regions: "
	cSNVs=$(cat $SNVs_coding_regions_stats | awk '{if ($5=="SNPs:"){print $NF}}')
	echo $cSNVs
	echo -n "chr${num} ts/tv cSNVs: "
	ts_tv_cSNVs=$(cat $SNVs_coding_regions_stats | grep -A1 "ts/tv" | tail -1 | cut -f5)
	echo $ts_tv_cSNVs

	# SNPs statistics of coding regions's
	if [ ! -f $SNPs_coding_regions_stats ]; then
		bcftools view -i 'AF>0.01' -R $coding_regions_bed $whole_chromosome_vcf | bcftools stats > $SNPs_coding_regions_stats
	fi
	echo -n "chr${num} SNPs in coding regions: "
	cSNP=$(cat  | awk '{if ($5=="SNPs:"){print $NF}}')
	echo $cSNP
	echo -n "chr${num} ts/tv cSNPs: "
	ts_tv_cSNPs=$(cat $SNPs_coding_regions_stats | grep -A1 "ts/tv" | tail -1 | cut -f5)
	echo $ts_tv_cSNPs

	# Measures of variability
	echo -n "chr${num} cSNVs frequency: "
	cSNVs_frequency=$(echo "scale=2; 100*$cSNVs/$SIZE_coding_regions" | bc -l)
	echo $cSNVs_frequency

	echo -n "chr${num} SNVs frequency: "
	SNVs_frequency=$(echo "scale=2; 100*$SNVs/$SIZE_chromosome" | bc -l)
	echo $SNVs_frequency

	echo -n "chr${num} SNPs_frequency: "
	SNPs_frequency=$(echo "scale=2; 100*$SNPS_variants/$SIZE_chromosome" | bc -l)
	echo $SNPs_frequency

	echo -n "chr${num} cSNPs frequency: "
	cSNPs_frequency=$(echo "scale=2; 100*$cSNPs/$SIZE_coding_regions" | bc -l)
	echo $cSNPs_frequency

	# print all chromosome information
	echo -ne "$CHR\t$SIZE_chromosome\t$SIZE_coding_regions\t" >> results.txt
	echo -ne "$SNVs\t$ts_tv_SNVs\t" >> results.txt
	echo -ne "$SNPs\t$ts_tv_SNPs\t" >> results.txt
	echo -ne "$cSNVs\t$ts_tv_cSNVs\t" >> results.txt
	echo -ne "$SNVs_frequency\t$cSNVs_frequency\t$SNPs_frequency\t$cSNPs_frequency\n" >> results.txt
done
