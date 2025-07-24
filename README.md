# hetmers

Rust script for working with heterozygous k-mers in pooled sequencing data for population genetics

# Usage

```
Usage: hetmers [OPTIONS] --inputs <INPUTS>... --outputs <OUTPUTS>... --minimums <MINIMUMS>... --coverages <COVERAGES>... --pools <POOLS>... --alphas <ALPHAS>... --betas <BETAS>... --sigmas <SIGMAS>...

Options:
  -i, --inputs <INPUTS>...        kmer count table file name
  -o, --outputs <OUTPUTS>...      prefix for output files
  -m, --minimums <MINIMUMS>...    analysis type to run minimum k-mer count
  -l, --alleles <ALLELES>         number of alleles in each hetmer [default: 2]
  -c, --coverages <COVERAGES>...  mean k-mer coverage
  -p, --pools <POOLS>...          pool size
  -a, --alphas <ALPHAS>...        shape parameter for prior distribution
  -b, --betas <BETAS>...          shape parameter for prior distribution
  -s, --sigmas <SIGMAS>...        thresholds for determining if k-mer has abnormal copy number
  -h, --help                      Print help
  -V, --version                   Print version
```

# Inputs

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
