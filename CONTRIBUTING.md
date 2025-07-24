
# Table of Contents

# Tutorial

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

- [x] get population specific hetmers

- [x] get hetmers shared across populations

- [x] write some tests

- [x] add some error checking code for inputs

- [x] check letters in k-mers

- [x] tag hetmers with really high total coverage

- [ ] verify that candidate hetmers (based on hash functions) actually differ

- [ ] include read length in k-mer coverage calculation?

- [ ] what if there's a duplicate k-mer?

- [ ] remove support for > 2 alleles

- [ ] maximum k-mer count?

- [ ] use hash from kmers_to_hetmers to find shared hetmers in two different samples

- [ ] bayesian pool size (is n constant for each k-mer)

- [ ] confidence intervals for allele frequencies

- [ ] credible intervals for allele frequencies

- [ ] Fst

- [ ] dxy

- [ ] Diversity

- [ ] Tajimas D for individual hetmers

- [ ] frequency increment test

- [ ] HWE?
