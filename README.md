# hetmers

Rust script for working with heterozygous k-mers in pooled sequencing data for population genetics

# Table of Contents

# Usage

```
Usage: hetmers [OPTIONS] --inputs <INPUTS>... --outputs <OUTPUTS>... --analyses <ANALYSES>... --minimums <MINIMUMS>... --coverages <COVERAGES>... --pools <POOLS>... --alphas <ALPHAS>... --betas <BETAS>...

Options:
  -i, --inputs <INPUTS>...        kmer count table file name
  -o, --outputs <OUTPUTS>...      prefix for output files
  -y, --analyses <ANALYSES>...    analysis type to run [possible values: hetmers, emp_freq, bayes_freq, fst, dxy, fit]
  -m, --minimums <MINIMUMS>...    minimum k-mer count
  -l, --alleles <ALLELES>         number of alleles in each hetmer [default: 2]
  -c, --coverages <COVERAGES>...  mean k-mer coverage
  -p, --pools <POOLS>...          pool size
  -a, --alphas <ALPHAS>...        shape parameter for prior distribution (bayes_freq analysis)
  -b, --betas <BETAS>...          shape parameter for prior distribution (bayes_freq_analysis)
  -h, --help                      Print help
  -V, --version                   Print version
```

# Tutorial

## 1. get pooled sequening data

### create pseudo pools

If you have sequencing data for individuals, you can create pseudo-pooled sequencing data by downsampling every individual to the same sequencing depth and then concatenating the downsampled fastq files.

## 2. count k-mers

You can use jellyfish or kmc or any other tool that counts k-mers and will dump the counts to a tab separated text file with two columns (canonical kmer sequence, kmer count)

It's typically a good idea at this point to examine the k-mer frequency spectrum.

## 3. use `hetmers`

You will want to know a few things:

* pool size

* minimum k-mer count to use (look at the k-mer frequency spectrum)

* average k-mer coverage ((L-k+1)/L)*(N*L/(G))

### Analysis of 1 population

`hetmers --inputs test.tsv --minimums 1 --alleles 2 --outputs yay --coverages 100 --pools 50 --alphas 1 --betas 1 --analyses hetmers`

### Analysis of 2 populations

#### Finding hetmers shared between two populations

#### Finding hetmers specific to one population

#### Fst

If you want to perform the same analysis on more than one population, then you can pass vectors to each argument.

`hetmers --inputs test.tsv test2.tsv --minimums 1 1 --alleles 2 --outputs yay yay2 --coverages 100 50 --pools 50 25 --alphas 1 1 --betas 1 1 --analyses fst`

### Analysis of >2 populations

# Examples

## Find putative SNPs in one pool

## Find putative SNPs specific to one pool (not shared with other pools)

## Measure Fst for hetmers shared between two pools

# To-do

- [x] empirical allele frequencies

- [x] bayesian allele states

- [x] loop over multiple populations separately

- [ ] bayesian pool size (is n constant for each k-mer)

- [ ] confidence intervals for allele frequencies

- [ ] credible intervals for allele frequencies

- [ ] Fst

- [ ] dxy

- [ ] Diversity

- [ ] Tajimas D for individual hetmers

- [ ] frequency increment test

- [ ] HWE?
