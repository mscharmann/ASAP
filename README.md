# ASAP
Allele-Specific Abundance and imPrinting pipeline

inputs:
- reference genome
- GFF including features named "gene" (these could be actual genes, or anything else).
- sequencing reads for progenies and their parents (can be either RNA-seq or DNA-seq, single or paired-end data). One or multiple read file units per sample are possible.
- info about the samples: pedigree, pool/individual and expected maternal read proportion

outputs:
- mapped reads by HISAT2 in .BAM format
- read mapping statistics
- variants in .VCF format
- SNP-level parental allele counts for each progeny
- gene-level parental allele statistics for each progeny, including chi-squared test for deviation from expected maternal proportion

## 0. setup
create a conda environment
```
conda create --name imprinting
conda activate imprinting

conda install snakemake samtools seqtk bcftools bedtools tabix pigz hisat2 scipy "vcflib=1.0.3" -y
```
clone this repo
```
git clone https://github.com/mscharmann/ASAP
```

download this helper script to directory "scripts"
```
wget -P scripts/ https://gist.githubusercontent.com/travc/0c53df2c8eca81c3ebc36616869930ec/raw/eff3032ca7c955ca33bffd8758092e4006949c75/split_ref_by_bai_datasize.py
```

## 0. prepare
- adjust parameters and input files in config.yaml
- set up the input read files to samples map in "data/samples_units_readfiles.txt"
- set up "data/pedigree.txt" with info about the pedigree, type of sample and expected maternal read proportion

- the repo comes with some simulated test data (see below to repeat the simulation). Use this to test your installation, or remove it like so:
```
rm data/*.fa data/*.fq.gz
```

## 1. run pipeline
- work within the above conda environment
- locally do:
```
snakemake -j 48 --restart-times 3 --keep-going --rerun-incomplete
```
- else on clusters with job submission systems like SLURM or LSF, use an appropriate snakemake command

## N. Simulate test data
- fake GFFs are provided, made by hand: data/fake.*.gff
- python, gffread and wgsim (from samtools) are used to simulate genomes and sequencing data
```
import numpy as np

length_of_the_chrom = 1000000

refchrom = "".join( np.random.choice(["A","C","T","G"], size = length_of_the_chrom, replace = True) )

# and SNPs (3%)
haplo_with_all_SNPs = refchrom
snpsites = np.random.uniform(0, len(refchrom), size = int(0.03*len(refchrom)) )
print (len(snpsites))
for s in snpsites:
	s = int(s)
	isnuc = refchrom[s]
	newnuc = np.random.choice([ x for x in ["A","C","T","G"] if not x == isnuc])
	haplo_with_all_SNPs = haplo_with_all_SNPs[:s] + newnuc + haplo_with_all_SNPs[s+1:]


# now make four random haplotypes as mixes of these 'perfectly divergent' ones:
h1 = ""
h2 = ""
h3 = ""
h4 = ""
for i,j in zip(refchrom,haplo_with_all_SNPs):
	nucs = i+j
	draws = np.random.choice([0,1], size = 4, replace = True)
	h1 += nucs[draws[0]]
	h2 += nucs[draws[1]]
	h3 += nucs[draws[2]]
	h4 += nucs[draws[3]]


with open("data/fakegenome.ref.fa", "w") as O:
	O.write(">fake_chromosome" + "\n")
	O.write(refchrom + "\n")

with open("data/parent1.diploid.fa", "w") as O:
	O.write(">fake_chromosome_1" + "\n")
	O.write(h1 + "\n")
	O.write(">fake_chromosome_2" + "\n")
	O.write(h2 + "\n")

with open("data/parent2.diploid.fa", "w") as O:
	O.write(">fake_chromosome_1" + "\n")
	O.write(h3 + "\n")
	O.write(">fake_chromosome_2" + "\n")
	O.write(h4 + "\n")

# now simulate the progeny's genomes; their haplotypes are produced by cross-over between the parent's haplotypes
with open("progeny1.diploid.fa", "w") as O:
	O.write(">data/fake_chromosome_1" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h1[:co] + h2[co:] + "\n")
	O.write(">fake_chromosome_2" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h3[:co] + h4[co:] + "\n")

with open("data/progeny2.diploid.fa", "w") as O:
	O.write(">fake_chromosome_1" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h1[:co] + h2[co:] + "\n")
	O.write(">fake_chromosome_2" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h3[:co] + h4[co:] + "\n")

with open("data/progeny3.diploid.fa", "w") as O:
	O.write(">fake_chromosome_1" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h1[:co] + h2[co:] + "\n")
	O.write(">fake_chromosome_2" + "\n")
	co = np.random.choice(range(length_of_the_chrom), 1)[0]
	O.write(h3[:co] + h4[co:] + "\n")
```

get 20x WGS coverage from each of the parents:
```
wgsim -N 66667 -1 150 -2 150 data/parent1.diploid.fa parent1.1.fastq parent1.2.fastq
wgsim -N 66667 -1 150 -2 150 data/parent2.diploid.fa parent2.1.fastq parent2.2.fastq

gzip *.fastq
```
now simulate some RNA-seq data for the progeny. Lets say endosperms are triploid with 2:1 M:P expectation
- 18 genes total, 1000 reads per gene = 18000 reads total. triploid endosperm.
- 6 MEGs = 900 maternal reads, 100 paternal reads: 900*6=5400 reads mat, 100*6=600 reads pat.
- 6 PEGs = 400 maternal reads, 600 paternal reads: 400*6=2400 reads mat, 600*6=3600 reads pat.
- 6 unimprinted = 667 maternal reads, 333 paternal reads: 667*6=4002 reads mat, 333*6=1998 reads pat.

These are given in the GFF files supplied. See here for GFF specifications:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


now make reads for progeny 1
```
samtools faidx data/progeny1.diploid.fa fake_chromosome_1 | sed 's/fake_chromosome_1/fake_chromosome/g' > tmp1.fa
samtools faidx data/progeny1.diploid.fa fake_chromosome_2 | sed 's/fake_chromosome_2/fake_chromosome/g' > tmp2.fa
gffread -g tmp1.fa -x tmp1.MEGs.mrna.fa data/fake.MEGs.gff
gffread -g tmp2.fa -x tmp2.MEGs.mrna.fa data/fake.MEGs.gff
wgsim -N 5400 -1 150 -2 150 tmp1.MEGs.mrna.fa tmp1.MEGs.1.fastq tmp1.MEGs.2.fastq
wgsim -N 600 -1 150 -2 150 tmp2.MEGs.mrna.fa tmp2.MEGs.1.fastq tmp2.MEGs.2.fastq

gffread -g tmp1.fa -x tmp1.PEGs.mrna.fa data/fake.PEGs.gff
gffread -g tmp2.fa -x tmp2.PEGs.mrna.fa data/fake.PEGs.gff
wgsim -N 2400 -1 150 -2 150 tmp1.PEGs.mrna.fa tmp1.PEGs.1.fastq tmp1.PEGs.2.fastq
wgsim -N 3600 -1 150 -2 150 tmp2.PEGs.mrna.fa tmp2.PEGs.1.fastq tmp2.PEGs.2.fastq

gffread -g tmp1.fa -x tmp1.not_imprinted.mrna.fa data/fake.not_imprinted.gff
gffread -g tmp2.fa -x tmp2.not_imprinted.mrna.fa data/fake.not_imprinted.gff
wgsim -N 4002 -1 150 -2 150 tmp1.not_imprinted.mrna.fa tmp1.not_imprinted.1.fastq tmp1.not_imprinted.2.fastq
wgsim -N 1998 -1 150 -2 150 tmp2.not_imprinted.mrna.fa tmp2.not_imprinted.1.fastq tmp2.not_imprinted.2.fastq

rm progeny1.1.fastq progeny1.2.fastq
for i in tmp1.MEGs  tmp1.not_imprinted  tmp1.PEGs  tmp2.MEGs  tmp2.not_imprinted  tmp2.PEGs ; do
	cat $i.1.fastq >> progeny1.1.fastq
	cat $i.2.fastq >> progeny1.2.fastq
done

rm tmp*
```
and for progeny 2
```
samtools faidx data/progeny2.diploid.fa fake_chromosome_1 | sed 's/fake_chromosome_1/fake_chromosome/g' > tmp2.fa # here we turn mother/father around
samtools faidx data/progeny2.diploid.fa fake_chromosome_2 | sed 's/fake_chromosome_2/fake_chromosome/g' > tmp1.fa # here we turn mother/father around
gffread -g tmp1.fa -x tmp1.MEGs.mrna.fa data/fake.MEGs.gff
gffread -g tmp2.fa -x tmp2.MEGs.mrna.fa data/fake.MEGs.gff
wgsim -N 5400 -1 150 -2 150 tmp1.MEGs.mrna.fa tmp1.MEGs.1.fastq tmp1.MEGs.2.fastq
wgsim -N 600 -1 150 -2 150 tmp2.MEGs.mrna.fa tmp2.MEGs.1.fastq tmp2.MEGs.2.fastq

gffread -g tmp1.fa -x tmp1.PEGs.mrna.fa data/fake.PEGs.gff
gffread -g tmp2.fa -x tmp2.PEGs.mrna.fa data/fake.PEGs.gff
wgsim -N 2400 -1 150 -2 150 tmp1.PEGs.mrna.fa tmp1.PEGs.1.fastq tmp1.PEGs.2.fastq
wgsim -N 3600 -1 150 -2 150 tmp2.PEGs.mrna.fa tmp2.PEGs.1.fastq tmp2.PEGs.2.fastq

gffread -g tmp1.fa -x tmp1.not_imprinted.mrna.fa data/fake.not_imprinted.gff
gffread -g tmp2.fa -x tmp2.not_imprinted.mrna.fa data/fake.not_imprinted.gff
wgsim -N 4002 -1 150 -2 150 tmp1.not_imprinted.mrna.fa tmp1.not_imprinted.1.fastq tmp1.not_imprinted.2.fastq
wgsim -N 1998 -1 150 -2 150 tmp2.not_imprinted.mrna.fa tmp2.not_imprinted.1.fastq tmp2.not_imprinted.2.fastq

rm progeny2.1.fastq progeny2.2.fastq
for i in tmp1.MEGs  tmp1.not_imprinted  tmp1.PEGs  tmp2.MEGs  tmp2.not_imprinted  tmp2.PEGs ; do
	cat $i.1.fastq >> progeny2.1.fastq
	cat $i.2.fastq >> progeny2.2.fastq
done

rm tmp*
```

and finally for progeny 3
```
samtools faidx data/progeny3.diploid.fa fake_chromosome_1 | sed 's/fake_chromosome_1/fake_chromosome/g' > tmp1.fa # here we turn mother/father around
samtools faidx data/progeny3.diploid.fa fake_chromosome_2 | sed 's/fake_chromosome_2/fake_chromosome/g' > tmp2.fa # here we turn mother/father around
gffread -g tmp1.fa -x tmp1.MEGs.mrna.fa data/fake.MEGs.gff
gffread -g tmp2.fa -x tmp2.MEGs.mrna.fa data/fake.MEGs.gff
wgsim -N 5400 -1 150 -2 150 tmp1.MEGs.mrna.fa tmp1.MEGs.1.fastq tmp1.MEGs.2.fastq
wgsim -N 600 -1 150 -2 150 tmp2.MEGs.mrna.fa tmp2.MEGs.1.fastq tmp2.MEGs.2.fastq

gffread -g tmp1.fa -x tmp1.PEGs.mrna.fa data/fake.PEGs.gff
gffread -g tmp2.fa -x tmp2.PEGs.mrna.fa data/fake.PEGs.gff
wgsim -N 2400 -1 150 -2 150 tmp1.PEGs.mrna.fa tmp1.PEGs.1.fastq tmp1.PEGs.2.fastq
wgsim -N 3600 -1 150 -2 150 tmp2.PEGs.mrna.fa tmp2.PEGs.1.fastq tmp2.PEGs.2.fastq

gffread -g tmp1.fa -x tmp1.not_imprinted.mrna.fa data/fake.not_imprinted.gff
gffread -g tmp2.fa -x tmp2.not_imprinted.mrna.fa data/fake.not_imprinted.gff
wgsim -N 4002 -1 150 -2 150 tmp1.not_imprinted.mrna.fa tmp1.not_imprinted.1.fastq tmp1.not_imprinted.2.fastq
wgsim -N 1998 -1 150 -2 150 tmp2.not_imprinted.mrna.fa tmp2.not_imprinted.1.fastq tmp2.not_imprinted.2.fastq

# this one to be split in two units of read files, just for the sake of demonstration:
rm progeny3.qwrtzy123.1.fastq progeny3.qwrtzy123.2.fastq
for i in tmp1.MEGs  tmp1.not_imprinted  tmp1.PEGs ; do
	cat $i.1.fastq >> progeny3.qwrtzy123.1.fastq
	cat $i.2.fastq >> progeny3.qwrtzy123.2.fastq
done

rm progeny3.gklb486.1.fastq progeny3.gklb486.2.fastq
for i in tmp2.MEGs  tmp2.not_imprinted  tmp2.PEGs ; do
	cat $i.1.fastq >> progeny3.gklb486.1.fastq
	cat $i.2.fastq >> progeny3.gklb486.2.fastq
done

rm tmp*
```
in the end, clean up and move all simulated reads to data/
```
gzip *.fastq
mv *.fastq.gz data/
```