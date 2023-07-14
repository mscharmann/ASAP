configfile: "config.yaml"

genomefile = config["genomefile"]
samples_units_fqs_map = config["samples_units_fqs_map"]
pedigree_file = config["pedigree_file"]
gff = config["gff"]

import pandas as pd
import gzip
import scipy

samples_units_fqs = pd.read_table(samples_units_fqs_map, dtype=str).set_index(
	["sample", "unit", "fq1", "fq2"], drop=False)

SAMPLES = list( set(samples_units_fqs["sample"]) )

pedigree_map = pd.read_table(pedigree_file, dtype=str).set_index(
	["progeny"], drop=False)

PROGENIES = list( set(pedigree_map["progeny"]) )


def get_sample_bams(wildcards):
	"""Get all aligned reads of given sample."""
	return expand(
		"mapped_reads_per_unit/{sample}-{unit}.sorted.bam",
		sample=wildcards.sample,
		unit=samples_units_fqs.loc[wildcards.sample].unit,
	)


def get_fastq_sample_unit(wildcards):
	"""Get fastq files of given sample-unit."""
	fastqs = samples_units_fqs.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]
	if fastqs.fq2.isnull().values.any():
		return [ fastqs.fq1.item() ]
	return [ fastqs.fq1.item(), fastqs.fq2.item() ]


def get_fastq_sample_ALL(wildcards):
	"""Get list of fastq files of given sample including ALL units."""
	fastqs_pd = samples_units_fqs.loc[(wildcards.sample), ["fq1", "fq2"]]
	fastqs = set( fastqs_pd.fq1.tolist() + fastqs_pd.fq2.tolist() )
	fastqs_clean = list( {x for x in fastqs if pd.notna(x)} )
	return fastqs_clean


rule all:
	input:
		expand("results/{progeny}.gene_level_parental_allele_stats.txt", progeny=PROGENIES),
		"results/mapping_statistics_report.txt"
	run:
		import os

		os.system( "rm -r varcall*" )
		os.system( "rm -r mapped_reads_per_unit" )


rule bwa_idx:
	input:
	   genomefile
	output:
		"bwa_index/reference.fa.bwt",
		"bwa_index/reference.fa
	shell:
		"""
		if [[ ! $( grep ">" {input} ) =~ "|" ]]; then
			mkdir -p bwa_index
			cd bwa_index
			cp {input} ./reference.fa
			bwa index reference.fa
		else
			echo "refusing to run, fasta headers contain pipe '|' character, dying"
		fi
		"""

rule bwa_map:
	input:
		fa="bwa_index/reference.fa",
		gidx="bwa_index/reference.fa.bwt",
		reads=get_fastq_sample_unit
	output:
		temp("mapped_reads_per_unit/{sample}-{unit}.bam")
	threads: 16
	run:
		if len(input.reads) == 2: # paired-end!
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen):
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# sum of the bit flags: 2304 => filters against BOTH non-primary and supplementary alignments; verified with samtools flagstat
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} {input.reads[1]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2304 -q 1 -b -@ 2 - > {output}
				""")
		else: # single-end
			shell("""
				# filtering alignments to be primary (i.e. each read only once; if multiple locations equally possible than a random one is chosen):
				# -F 256 == -F 0x0100 == NOT not primary alignment
				# filtering alignments to be NOT supplementary (supplementary: sections of read map to discontinuous coordinates, e.g. across an inversion breakpoint..):
				# -F 2048 == -F 0x800 == NOT supplementary alignment
				# -F 4 read unmapped (0x4)
				# sum of the bit flags: 2308 => filters against non-primary and supplementary alignments and unmapped
				# filtering against multi-mapping alignments (which have MAPQ=0): -q 1
				bwa mem -t {threads} -a {input.fa} {input.reads[0]} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina" | samtools view -F 2308 -q 1 -b -@ 2 - > {output}
				""")


rule samtools_sort:
	input:
		"mapped_reads_per_unit/{sample}-{unit}.bam"
	output:
		temp( "mapped_reads_per_unit/{sample}-{unit}.sorted.bam" )
	shell:
		"""
		samtools sort -T mapped_reads_per_unit/{wildcards.sample}.{wildcards.unit} -O bam {input} > {output}
		"""

rule merge_bams_per_sample:
	input:
		bams=get_sample_bams
	output:
		"mapped_reads/{sample}.sorted.bam"
	threads: 4
	shell:
		"""
		samtools merge --threads {threads} {output} {input.bams}
		"""


rule samtools_index:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		"mapped_reads/{sample}.sorted.bam.bai"
	shell:
		"samtools index {input}"




# A checkpoint that shall trigger re-evaluation of the DAG
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
# we need this because the number of output files (i.e. chunks / region files for variant calling) is not known before execution!
checkpoint split_ref_for_varcall:
	input:
		ref=genomefile,
		indexes=expand("mapped_reads/{sample}.sorted.bam.bai", sample=SAMPLES),
		samples=expand("mapped_reads/{sample}.sorted.bam", sample=SAMPLES)
	output:
		directory("varcall_chunks")
	params:
		chunksize_Mb = config["varcall_chunksize"]
	shell:
		"""
		samtools faidx {input.ref}

		echo {input.samples} | sed 's/ /\\n/g' > bamlistforsplit
		# wget https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
		python scripts/split_ref_by_bai_datasize.py -r {input.ref}.fai -L bamlistforsplit --target-data-size {params.chunksize_Mb} > varcall_regionsALL.bed
		rm bamlistforsplit
		split --lines=1 varcall_regionsALL.bed varcall_regions_ --numeric-suffixes --suffix-length=6 --additional-suffix=.bed
		mkdir -p {output}
		for i in varcall_regions_* ; do mv $i {output} ; done
		"""

rule call_variants:
	input:
		ref=genomefile,
		regions="varcall_chunks/{i}.bed",
		gff=config["gff"]
	output:
		temp( "varcall_chunk_VCFs/{i}.bed.vcf.gz" )
	shell:
		"""
		# call variants ONLY in GFF regions:
		bedtools intersect -a {input.regions} -b {input.gff} | sort -k 1,1 -k2,2n | bedtools merge > {input.regions}.gff_intersected.bed

		# handle regions not overlapping with any genes, returning empty bed files for bcftools.
		if [ -s {input.regions}.gff_intersected.bed ] ; then
			# --max-depth 1000: use at most 1000 reads per input BAM file, apparently these are sampled RANDOMLY!?
			# --min-MQ 20: minimum mapping quality of an alignment, otherwise skip
			# --no-BAQ : do NOT re-calculate mapping quality (which involves re-aligning). Instead, will use MAPQ as stored in BAM file.
			# --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
			# --variants-only
			bcftools mpileup -Ou -f {input.ref} -R {input.regions}.gff_intersected.bed --bam-list <( ls mapped_reads/*.bam ) --max-depth 1000 --min-MQ 20 --min-BQ 15 --no-BAQ -a INFO/AD -a FORMAT/AD | bcftools call -m --variants-only --skip-variants indels -Ov | bgzip -c > {output}
		else
			# input region is empty, therefore run without a specific region but return the HEADER of the VCF only, then terminate. VCF with at least the header are required for downstream; empty files will break pipeline.
			# ensure that output directory exists.. and force return exist status 0, otherwise it will often break. Unclear why.
			mkdir -p varcall_chunk_VCFs
			bcftools mpileup -Ou -f {input.ref} --bam-list <( ls mapped_reads/*.bam ) --max-depth 1000 --min-MQ 20 --min-BQ 15 --no-BAQ -a INFO/AD -a FORMAT/AD | bcftools call -m --variants-only --skip-variants indels -Ov | awk '{{if ($1 == "#CHROM")  {{print ; exit;}} else print}}' | bgzip -c > {output} 2>&1 || true
		fi
		rm {input.regions}.gff_intersected.bed
		"""


def merge_vcfs_input(wildcards):
	checkpoint_output = checkpoints.split_ref_for_varcall.get(**wildcards).output[0]
	# print(checkpoint_output)
	vcfs_before_filter = sorted( expand("varcall_chunk_VCFs/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	# print (thef)
	vcfs_after_filter = sorted( expand("varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bed")).i) )
	return vcfs_after_filter



rule merge_filtered_vcfs:
	input:
		region_vcfs=merge_vcfs_input # refers to the function above which evaluates the checkpoint
	output:
		"results/variants.post_filter.vcf.gz"
	threads: 3
	shell:
		"""
		# MUST NOT USE bcftools concat: it cannot resolve POS that are non-monotonically icreasing (which ca happen at the interval boundaries)
		# the code below is only slightly modified from freebayes-parallel script: https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel
		# zcat input.region_vcfs | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		# zcat alone may complain about too many arguments, so better use find -exec :
		find varcall_chunk_VCFs_filtered/*.bed.vcf.gz -type f -exec zcat {{}} \\; | python $(which vcffirstheader) | vcfstreamsort -w 10000 | vcfuniq | bgzip -c > {output}
		sleep 20
		tabix {output}
		"""



rule VCF_filter_variants:
	input:
		gzvcf="varcall_chunk_VCFs/{i}.bed.vcf.gz"
	output:
		temp( "varcall_chunk_VCFs_filtered/{i}.bed.vcf.gz" )
	params:
		QUAL=config["VCF_QUAL"],
		MIN_DEPTH=config["MIN_DEPTH"]
	shell:
		"""
		## https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/#comment-1469
		## => "We use ... a very vanilla filter with only depth and quality (QUAL < 20, DP < 5)"
		## => impose a minimum of QUAL, and a minimum of DEPTH
		## 	NOT applying maximum DPETH cutoff
		## mac = minimum minor allele count: 2, because singletons are not relevant here (informative variants must occur in at least one parent and its progeny)
		## keep only bi-alleleic sites: --min-alleles 2 --max-alleles 2

		# prepare
		wd=DIR_{wildcards.i}
		mkdir -p varcall_chunk_VCFs_filtered/$wd
		cd varcall_chunk_VCFs_filtered/$wd

		# set filters

		# variants: || "vcftools --min-meanDP" is a bad idea here, because it is a site-filter, not a genotype-filter: some low-coverage samples in the set may cause mean depth across samples to be low, thus loosing sites that are good and credible in the other samples!
		vcftools --gzvcf ../../{input.gzvcf} --mac 2 --min-alleles 2 --max-alleles 2 --minQ {params.QUAL} --minDP {params.MIN_DEPTH} --recode --stdout | bgzip -c > tmp.1
		tabix tmp.1

		# filter for missingness: keep any site with at least two samples present
		# first find out how many samples there are:
		nsamples=$( bcftools query -l tmp.1 | wc -l )
		echo $nsamples
		total_minus_two=$((nsamples-2))
		echo $total_minus_two
		vcftools --gzvcf tmp.1 --max-missing-count $total_minus_two --mac 2 --recode --stdout | bcftools view --exclude-uncalled --trim-alt-alleles | bgzip -c > ../../{output}

		# cleanup
		cd ../
		rm -r $wd

		"""


rule count_reads_input_and_mapping:
	input:
		indexfile="mapped_reads/{sample}.sorted.bam.bai",
		bamfile="mapped_reads/{sample}.sorted.bam",
		reads=get_fastq_sample_ALL
	output:
		temp( "results/{sample}.raw_and_mapped_reads_report.txt" )
	threads: 3
	run:
		cmd = "samtools flagstat " + input.bamfile + " > " + input.bamfile + ".flagstat.txt"
		os.system( cmd )
		with open(input.bamfile +".flagstat.txt", "r") as I:
			total_mapped = int(I.readline().split()[0])
		os.system( "rm " +  input.bamfile +".flagstat.txt")
		# reconstruct the sample name...
		samplename = input.bamfile.split("/")[1].split(".sorted.bam")[0]
		# now count all the reads
		sum_of_read_fastq_lines = 0
		for f in list(input.reads):
			os.system( "pigz -p3 -dc {0} | wc -l > tmp.{1}.readcount".format(f,samplename) )
			with open( "tmp.{}.readcount".format(samplename), "r" ) as an_infile:
				sum_of_read_fastq_lines += int(an_infile.readline().strip())
			os.system( "rm tmp.{}.readcount".format(samplename) )
		n_reads = float(sum_of_read_fastq_lines)/4.0
		mapping_rate = float(total_mapped)/n_reads
		with open(str(output), "w") as O:
			O.write(samplename + "\t" + str(n_reads) + "\t" + str(total_mapped) + "\t" + str(round(mapping_rate,3)) + "\n")


rule collect_mapping_report:
	input:
		expand("results/{sample}.raw_and_mapped_reads_report.txt", sample=SAMPLES)
	output:
		"results/mapping_statistics_report.txt"
	shell:
		"""
		echo -e "counts are single reads, not pairs of reads" > {output}
		echo -e "sample\traw_reads\tmapped_reads\tmapping_rate" >> {output}
		cat results/*.raw_and_mapped_reads_report.txt >> {output}
		"""

rule get_raw_parental_counts:
	"""
	# this is the heart of parent-specific allele counting
	there are 9 possible m+f genotype combinations:
	m	f	comment
	0/0	0/0	never informative
	1/1	1/1	never informative
	0/1	0/1 never informative
	1/1	0/0	always informative
	0/0	1/1	always informative
	0/1	0/0	If progeny is individual: informative only if progeny het. If progeny is a (large) pool: informative with correction.
	0/1	1/1	If progeny is individual: informative only if progeny het. If progeny is a (large) pool: informative with correction.
	0/0	0/1	If progeny is individual: informative only if progeny het. If progeny is a (large) pool: informative with correction.
	1/1	0/1	If progeny is individual: informative only if progeny het. If progeny is a (large) pool: informative with correction.

	- important detail of scoring parent-of-origin expression at a certain type of sites:
		- one parent het, the other hom, progeny het.
			- if progeny is an individual: naive counts of the two alleles are directly the parental counts
			- if progeny is a pool, and assuming randomness at all steps from meioses to fertilization:
	 				the naive counts are not directly parental counts, because the heterozygote parent contributes
	 				BOTH alleles to the pool of progeny. Hence, apply a correction!
					- if mother is het: naive maternal count is a 2-fold underestimation
							corrected maternal count = naive maternal count * 2 | cap to the total count for consistency
							corrected paternal count = total count - corrected maternal count
					- if father is het: naive paternal count is a 2-fold underestimation
							corrected maternal count = naive maternal count - naive paternal count | cap to 0 for consistency
							corrected paternal count = total count - corrected maternal count
	"""
	input:
		vcf="results/variants.post_filter.vcf.gz"
	output:
		"results/{progeny}.raw_parental_counts.txt"
	run:
		parents_pd = pedigree_map.loc[(wildcards.progeny), ["mother", "father"]]
		type_of_sample = pedigree_map.loc[(wildcards.progeny), ["type_of_sample"]][0]

		outf = open( output[0] , "w")
		outf.write("\t".join(["chrom","pos","site_type","maternal_count","paternal_count"]) + "\n")

		with gzip.open( input.vcf ,"rt") as I:
			for line in I:
				if line.startswith("#CHROM"):
					header = line.strip("\n").split("\t")
					progeny_idx = header.index( wildcards.progeny )
					mother_idx = header.index( parents_pd.mother )
					father_idx = header.index( parents_pd.father )
				else:
					if not line.startswith("##"):
						fields = line.strip("\n").split("\t")
						chrom = fields[0]
						pos = fields[1]
						ref_allele = fields[3]
						alt_allele = fields[4]
						mother_gt = fields[mother_idx].split(":")[0]
						father_gt = fields[father_idx].split(":")[0]
						progeny_allelic_counts = [int(x) for x in fields[progeny_idx].split(":")[2].split(",")]
						# if the progeny is not heterozygous:## NO; MUST NOT DO THIS BECAUSE THIS excludes the extreme, perfectly allele-specific case!!
						if sum(progeny_allelic_counts) < 6: # the sum of counts must be at least 6, otherwise discard
							outline = [chrom, pos, "progeny_coverage_low", "NA", "NA"]
						else:
							if mother_gt == father_gt:# so this is 0/0,0/0 or 1/1,1/1 or 0/1,0/1
								# both parents hom for same alelle should already be impossible but, it can also be that both are het.
								# if parents are identical, site is not informative for this progeny
								outline = [chrom, pos, "parents_identical", "NA", "NA"]
							elif mother_gt == "0/0" and father_gt == "1/1":
								# mother and father each homozygous for different alleles, mother's allele here REF
								outline = [chrom, pos, "parents_reciprocally_hom", progeny_allelic_counts[0], progeny_allelic_counts[1]]
							elif mother_gt == "1/1" and father_gt == "0/0":
								# mother and father each homozygous for different alleles, mother's allele here ALT
								outline = [chrom, pos, "parents_reciprocally_hom", progeny_allelic_counts[1], progeny_allelic_counts[0]]
							# now dealing with scenario where only one parent is het. See above for the different handling of pools and individuals, the correction.
							elif type_of_sample == "individual":
								if progeny_allelic_counts[0] < 3 or progeny_allelic_counts[1] < 3:
									# in progeny individuals, must see both alleles for such sites to be informative
									outline = [chrom, pos, "one_parent_het_but_progeny_hom", "NA", "NA"]
								else:
									if mother_gt == "0/1" and father_gt == "1/1":
										# mother het and progeny het: mothers' allele is clear, here REF
										outline = [chrom, pos, "mother_het_father_hom", progeny_allelic_counts[0], progeny_allelic_counts[1]]
									elif mother_gt == "0/1" and father_gt == "0/0":
										# mother het and progeny het: mothers' allele is clear, here ALT
										outline = [chrom, pos, "mother_het_father_hom", progeny_allelic_counts[1], progeny_allelic_counts[0]]
									elif mother_gt == "0/0" and father_gt == "0/1":
										# father het and progeny het: mothers' allele is clear, here REF
										outline = [chrom, pos, "mother_hom_father_het", progeny_allelic_counts[0], progeny_allelic_counts[1]]
									elif mother_gt == "1/1" and father_gt == "0/1":
										# father het and progeny het: mothers' allele is clear, here ALT
										outline = [chrom, pos, "mother_hom_father_het", progeny_allelic_counts[1], progeny_allelic_counts[0]]
									else:
										outline = [chrom, pos, "unclear_problem", "NA", "NA"]
							elif type_of_sample == "pool": # now we consider informative even if progeny has only 1 allele, and apply a correction factor because both alleles must be present in the genomes of the progeny in the pool.
								if mother_gt == "0/1":
									if father_gt == "1/1":
										# mother het mothers' allele is clear, here RER
										mat_naive = int(progeny_allelic_counts[0])
										pat_naive = int(progeny_allelic_counts[1])
									elif father_gt == "0/0":
										# mother het mothers' allele is clear, here ALT
										mat_naive = int(progeny_allelic_counts[1])
										pat_naive = int(progeny_allelic_counts[0])
									mat_corrected = mat_naive*2
									if mat_corrected > mat_naive+pat_naive:
										mat_corrected = mat_naive+pat_naive
									pat_corrected = (mat_naive+pat_naive) - mat_corrected
									outline = [chrom, pos, "mother_het_father_hom_and_pool_correction", mat_corrected, pat_corrected]
								elif father_gt == "0/1":
									if mother_gt == "1/1":
										# mother het mothers' allele is clear, here REF
										mat_naive = int(progeny_allelic_counts[1])
										pat_naive = int(progeny_allelic_counts[0])
									elif mother_gt == "0/0":
										# mother het mothers' allele is clear, here ALT
										mat_naive = int(progeny_allelic_counts[0])
										pat_naive = int(progeny_allelic_counts[1])
									pat_corrected = pat_naive*2
									if pat_corrected > mat_naive+pat_naive:
										pat_corrected = mat_naive+pat_naive
									mat_corrected = (mat_naive+pat_naive) - pat_corrected
									outline = [chrom, pos, "mother_hom_father_het_and_pool_correction", mat_corrected, pat_corrected]
								else:
									outline = [chrom, pos, "unclear_problem", "NA", "NA"]
		#				print("###")
		#				print(mother_gt,father_gt,progeny_allelic_counts)
						outf.write( "\t".join([str(x) for x in outline]) + "\n")
		outf.close()


rule get_gene_level_parental_allele_counts:
	input:
		rawcounts="results/{progeny}.raw_parental_counts.txt",
		gff=config["gff"]
	output:
		temp( "results/{progeny}.gene_level_parental_allele_counts.txt" )
	shell:
		"""
		tail -n +2 {input.rawcounts} | awk '{{if($4 != "NA") print}}' | awk '{{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}}' | sort -k 1,1 -k2,2n > {input.rawcounts}.bedformat.bed
		cat {input.gff} | awk '{{if($3 == "gene") print}}' | bedtools sort > {input.rawcounts}.tmp.gff
		bedtools intersect -a {input.rawcounts}.tmp.gff -b {input.rawcounts}.bedformat.bed -c > {input.rawcounts}.snpcounts
		bedtools map -a {input.rawcounts}.tmp.gff -b {input.rawcounts}.bedformat.bed -c 4 -o mean > {input.rawcounts}.mat
		bedtools map -a {input.rawcounts}.tmp.gff -b {input.rawcounts}.bedformat.bed -c 5 -o mean > {input.rawcounts}.pat

		echo -e 'chrom\tstart\tend\tgene_name\tSNPs\tmean_maternal_count\tmean_paternal_count' > {output}
		paste <(cut -f1,4,5,9,10 {input.rawcounts}.snpcounts) <(cut -f10 {input.rawcounts}.mat) <(cut -f10 {input.rawcounts}.pat) | awk '{{ if($5 == "0") {{ print $1"\t"$2"\t"$3"\t"$4"\t0\t0\t0" }} else {{ print}} }}' >> {output}
		rm {input.rawcounts}.bedformat.bed {input.rawcounts}.tmp.gff {input.rawcounts}.snpcounts {input.rawcounts}.mat {input.rawcounts}.pat

		"""


rule get_gene_level_parental_allele_stats:
	input:
		counts="results/{progeny}.gene_level_parental_allele_counts.txt"
	output:
		"results/{progeny}.gene_level_parental_allele_stats.txt"
	run:
		infile = "results/progeny1.gene_level_parental_allele_stats.txt"
		expected_maternal_proportion = float( pedigree_map.loc[(wildcards.progeny), ["expected_maternal_proportion"]].iloc[0] )
		expected_paternal_proportion = 1.0-expected_maternal_proportion

		mydata = pd.read_table( input.counts , dtype=str)
		mydata[["mean_maternal_count", "mean_paternal_count"]] = mydata[["mean_maternal_count", "mean_paternal_count"]].apply(pd.to_numeric)
		mydata["rowsum"] = mydata["mean_maternal_count"] + mydata["mean_paternal_count"]
		mydata["maternal_proportion"] = mydata["mean_maternal_count"] / mydata["rowsum"]
		mydata["maternal_expect"] = mydata["rowsum"]*expected_maternal_proportion
		mydata["paternal_expect"] = mydata["rowsum"]*expected_paternal_proportion

		test_res = scipy.stats.chisquare(mydata[["mean_maternal_count","mean_paternal_count"]], mydata[["maternal_expect","paternal_expect"]], ddof=0, axis=1)

		# make the output cleaner
		mydata["mean_maternal_count"] = mydata["mean_maternal_count"].round(2)
		mydata["mean_paternal_count"] = mydata["mean_paternal_count"].round(2)
		mydata.drop(columns="rowsum", inplace=True)
		mydata["maternal_expect"] = mydata["maternal_expect"].round(2)
		mydata["paternal_expect"] = mydata["paternal_expect"].round(2)

		mydata["chisq_statistic"] = test_res[0].round(3)
		mydata["p_value"] = test_res[1]
		mydata["maternal_proportion"] = mydata["maternal_proportion"].round(4)

		# if there were no SNPs, do not report chisq statistic and p-value as 0.0, but as NA. same for descriptive stats.
		mydata["p_value"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)
		mydata["chisq_statistic"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)

		mydata["mean_maternal_count"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)
		mydata["mean_paternal_count"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)
		mydata["maternal_proportion"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)
		mydata["maternal_expect"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)
		mydata["paternal_expect"].mask(mydata["SNPs"] =="0" ,"NA", inplace=True)

		mydata.to_csv(output[0], sep='\t', header=True, index=False)
