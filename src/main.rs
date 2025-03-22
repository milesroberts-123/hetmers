use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::Parser; // argument parser
//use bio::alphabets::dna;
use sha2::{Digest, Sha256};
use itertools::izip;
//use std::collections::hash_map::Keys

// Structure of input arguments
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// kmer count table file name
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    inputs: Vec<String>,

    /// prefix for output files
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    outputs: Vec<String>,

    /// analysis type to run
    #[arg(short = 'y', long, num_args = 1.., value_delimiter = ' ', required = true, value_parser = clap::builder::PossibleValuesParser::new(["hetmers", "emp_freq", "bayes_freq", "fst", "dxy", "fit"]))]
    analyses: Vec<String>,

    /// minimum k-mer count
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    minimums: Vec<usize>,

    /// number of alleles in each hetmer
    #[arg(short = 'l', long, default_value_t = 2)]
    alleles: usize,

    /// mean k-mer coverage
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    coverages: Vec<f64>,

    /// pool size
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    pools: Vec<i32>,

    /// shape parameter for prior distribution (bayes_freq analysis)
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    alphas: Vec<f64>,

    /// shape parameter for prior distribution (bayes_freq_analysis)
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    betas: Vec<f64>,
}


// Function to read in kmer count table
fn load_kmers(input: &String, minimum: usize) -> (Vec<String>, Vec<usize>) {
    println!("Loading k-mer count file {}...", input);
    let file = File::open(input).expect("Unable to open file");
    let reader = BufReader::new(file);

    let mut seqs = Vec::new();
    let mut counts = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            continue;
        }
        let seq = parts[0].to_string();
        let count: usize = parts[1].parse().expect("Invalid count value");
        if count >= minimum {
            seqs.push(seq);
            counts.push(count);
        }
    }
    return(seqs, counts);
}

// get kmer without central bp
fn extract_border(seqs: Vec<String>) -> Vec<String> {
    let k = seqs[0].len();
    println!("k is {}", k);
    let k_half = k / 2;

    let no_center_ks: Vec<String> = seqs.iter()
        .map(|s| format!("{}{}", &s[..k_half], &s[k_half + 1..]))
        .collect();

    return no_center_ks;
}

// reverse complement sequences
fn rev_comp(seqs: Vec<String>) -> Vec<String> {
   println!("Reverse complementing...");
    let complement = |c: char| match c {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c
    };
    let revseqs: Vec<String> = seqs.iter()
        .map(|s| s.chars().rev().map(complement).collect())
        .collect();
    return revseqs;
}

// hash sequences
fn hash_seqs(seqs: Vec<String>) -> Vec<u64>{
    println!("Hashing...");
    let hash_fn = |s: &String| {
        let mut hasher = Sha256::new();
        hasher.update(s.as_bytes());
        u64::from_be_bytes(hasher.finalize()[..8].try_into().unwrap())
    };

    let seq_hashes: Vec<u64> = seqs.iter().map(hash_fn).collect();

    return seq_hashes;
}

// minimum hash
fn min_hash(hash1: Vec<u64>, hash2: Vec<u64>) -> Vec<u64>{
    println!("Getting the minimum hash...");
    let min_hashes: Vec<u64> = hash1.iter()
        .zip(hash2.iter())
        .map(|(x, y)| *x.min(y))
        .collect();

    return min_hashes;
}

// group hashes
fn group_hashes(hashes: Vec<u64>) -> HashMap<u64, Vec<usize>>{
    println!("Grouping unique hashes into a dictionary...");
    let mut d: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, num) in hashes.iter().enumerate() {
        d.entry(*num).or_insert_with(Vec::new).push(i);
    }

    return d;
}

// filter hashes by number of alleles
fn filter_groups(input: HashMap<u64, Vec<usize>>, alleles: usize) -> HashMap<u64, Vec<usize>> {
    println!("Filtering hash groups by number of alleles...");
    let output = input.iter()
        .filter(|(_, v)| v.len() == alleles)
        .map(|(k, v)| (*k, v.clone()))
        .collect();

    return output;
}


// extract hetmer sequences and counts based on hashes
fn extract_hetmers(hashdict: HashMap<u64, Vec<usize>>, seqs: Vec<String>, counts: Vec<usize>) -> (Vec<String>, Vec<String>, Vec<u64>) {
    println!("Extracting counts and sequences...");
    let hetmer_seqs: Vec<String> = hashdict.values()
        .map(|indices| indices.iter().map(|&i| seqs[i].clone()).collect::<Vec<String>>().join(","))
        .collect();
    let hetmer_counts: Vec<String> = hashdict.values()
        .map(|indices| indices.iter().map(|&i| counts[i].to_string()).collect::<Vec<String>>().join(","))
        .collect();

    return (hetmer_seqs, hetmer_counts, hashdict.into_iter().map(|(id, _score)| id).collect());
}

// write vector to file
fn write_file(output: Vec<String>, prefix: String, suffix: String) {
    println!("Saving results...");
    let mut file = File::create(format!("{}_{}", prefix, suffix)).expect("Unable to create file");
    writeln!(file, "{}", output.join("\n")).expect("Unable to write to file");
}

// Collect functions into one 
fn kmers_to_hetmers(input: &String, output: String, minimum: usize, alleles: usize) -> (Vec<std::string::String>, Vec<std::string::String>, Vec<u64>){
    let kmers = load_kmers(input, minimum);
    //println!("{:?}", kmers);

    let borders = extract_border(kmers.0.clone());
    //println!("{:?}", borders);

    let revborders = rev_comp(borders.clone());
    //println!("{:?}", revborders);

    let hashbord = hash_seqs(borders);
    //println!("{:?}", hashbord);

    let hashrevbord = hash_seqs(revborders);
    //println!("{:?}", hashrevbord);

    let min_hashes = min_hash(hashbord, hashrevbord);
    //println!("{:?}", min_hashes);

    let grouped_hashes = group_hashes(min_hashes);
    //println!("{:?}", grouped_hashes);

    let filtered_groups = filter_groups(grouped_hashes, alleles);
    //println!("{:?}", filtered_groups);

    let hetmers = extract_hetmers(filtered_groups, kmers.0, kmers.1);
    //println!("{:?}", hetmers);

    write_file(hetmers.0.clone(), output.clone(), "seqs.csv".to_string());
    write_file(hetmers.1.clone(), output.clone(), "counts.csv".to_string());
    write_file(hetmers.2.clone().iter().map(|num| num.to_string()).collect(), output, "hashes.csv".to_string());

    return hetmers;
}

// empirical allele frequencies
fn counts_to_frequencies(count_pairs: Vec<String>, output: String) -> Vec<f64> {
    println!("Calculating frequencies...");
    let frequencies: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let min_num = num1.min(num2);
                    let max_num = num1.max(num2);
                    let sum = min_num + max_num;
                    if sum != 0.0 {
                        return Some(min_num / sum);
                    }
                }
            }
        None
        })
        .collect();

    // write frequencies to file
    write_file(frequencies.clone().into_iter().map(|s| s.to_string()).collect(), output, "empirical_freqs.csv".to_string());
    // done
    return frequencies;
}

// bayesian allele states
fn posterior(x: f64, z: f64, n: i32, c: f64, alpha: f64, beta: f64) -> Vec<f64> {
    let mut likelihood_times_prior = Vec::new();
    let mut total_probability = 0.0;

    for i in 1..n {
        //println!("i is {}", i);
        let p = i as f64 / n as f64;
        //println!("p is {}", p);
        //println!("c is {}", c);
        let lambdax = (p as f64) * c;
        //println!("lambdax is {}",lambdax);
        let lambday = (1.0-p as f64) * c;
        //println!("lambday is {}",lambday);
        
        let num = p.powf(x + alpha - 1.0) * (1.0-p).powf((z - x) + beta - 1.0);
        let denom = (lambdax.exp() - 1.0) * (lambday.exp() - 1.0);        
        let value = num/denom;
        //println!("Numerator is {}",num);
        //println!("Denominator is {}",denom);
        //println!("Value is {}",value);
        likelihood_times_prior.push(value);
        total_probability += value;
    }
        
    let post: Vec<f64> = likelihood_times_prior.iter().map(|&v| v / total_probability).collect();

    return post;    
}

fn highest_prob_index(probabilities: Vec<f64>) -> Option<usize> {
    probabilities.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal)).map(|(index, _)| index)
}

fn counts_to_bayes_state(count_pairs: Vec<String>, n: i32, c: f64, alpha: f64, beta: f64, output: String) -> Vec<usize>{
    println!("Calculating posterior...");
    let bayes_states: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let min_num = num1.min(num2);
                    let max_num = num1.max(num2);
                    let z = min_num + max_num;
                    let post = posterior(min_num,z,n,c,alpha,beta);
                    return highest_prob_index(post);
                } else { return Some(0);}
            } else { return Some(0); }
        })
        .collect();
    // write frequencies to file
    write_file(bayes_states.clone().into_iter().map(|s| s.to_string()).collect(), output, "bayes_states.csv".to_string());
    return bayes_states;
}

// fst
fn shared_hetmers(map1: Vec<u64>, map2: Vec<u64>) -> Vec<u64>{
    println!("Find hetmers shared between two maps...");
    let output = map1.into_iter()
        .filter(|k| map2.contains(k))
        .collect();

    return output;
}

fn population_specific_hetmers(map1: &HashMap<u64, Vec<usize>>, map2: &HashMap<u64, Vec<usize>>) -> HashMap<u64, Vec<usize>>{
    println!("Find hetmers specific to one map...");
    let output = map1.iter()
        .filter(|(k, _)| !map2.contains_key(k))
        .map(|(k, v)| (*k, v.clone()))
        .collect();

    return output;
}

// fn kmers_to_hetmers

// frequency increment test

// test functions
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_k_borders() {
        let test_vec = vec!["ATGCA".to_string(), "TTGAT".to_string(), "GGATA".to_string()];
        let result = extract_border(test_vec);
        let expected = vec!["ATCA".to_string(), "TTAT".to_string(), "GGTA".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn even_k_borders() {
        let test_vec = vec!["ATGCAT".to_string(), "TTGATC".to_string(), "GGATAA".to_string()];
        let result = extract_border(test_vec);
        let expected = vec!["ATGAT".to_string(), "TTGTC".to_string(), "GGAAA".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn short_rev_comp() {
        let test_vec = vec!["ATGCAT".to_string(), "TTGATC".to_string(), "GGATAA".to_string()];
        let result = rev_comp(test_vec);
        let expected = vec!["ATGCAT".to_string(), "GATCAA".to_string(), "TTATCC".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn small_min_hash() {
        let test_vec_1 = vec![59888, 1, 100];
        let test_vec_2 = vec![59887, 5000, 101];
        let result = min_hash(test_vec_1, test_vec_2);
        let expected = vec![59887, 1, 100];
        assert_eq!(result, expected);
    }

}

// mix it all together! :D
fn main() {
    //let args: Vec<String> = env::args().collect();
    let args = Args::parse(); 

    //if args.analyses.contains(&"hetmers".to_string()) { 
    // loop over individual populations
    for (input, output, minimum, coverage, pool, alpha, beta) in izip!(&args.inputs, &args.outputs, &args.minimums, args.coverages, args.pools, args.alphas, args.betas) {
        // find hetmers
        //println!("{}", "1");
        let hetmers = kmers_to_hetmers(input, output.clone(), *minimum, args.alleles);

        // empirical hetmer frequencies
        let empirical_freqs = counts_to_frequencies(hetmers.1.clone(), output.clone());
        //println!("{:?}", empirical_freqs);
    
        // bayesian hetmer frequencies
        //if args.analyses.contains(&"bayes".to_string()) { 
        let bayes_states = counts_to_bayes_state(hetmers.1, pool, coverage, alpha, beta, output.to_string());
        //}
        //println!("{:?}", bayes_states);
    }
    //}

   // analyses comparing two populations, or pairs of populations
   //let num_pops = args.inputs.len();

   //if (num_pops == 2) && args.analyses.contains(&"fst".to_string()){
       //for i in 0..(num_pops-1){
       //    for j in (i+1)..num_pops{
               //println!("{:?}", &args.inputs[0]);
               //println!("{:?}", &args.inputs[1]);
               // get hetmers that are specific to each population
               //let map1 = kmers_to_hetmers(&args.inputs[0], args.outputs[0].clone(), args.minimums[0], args.alleles.clone());
               //let map2 = kmers_to_hetmers(&args.inputs[1], args.outputs[1].clone(), args.minimums[1], args.alleles);

               // find hetmers shared between the two populations (intersection of hash vector)
               //let hetmers_in_common = shared_hetmers(map1.2,map2.2);
               //println!("{:?}", hetmers_in_common

               //let map1_specific_hetmers = population_specific_hetmers(&map1,&map2);
               //println!("{:?}", map1_specific_hetmers);
               //let map2_specific_hetmers = population_specific_hetmers(&map2,&map1);
               //println!("{:?}", map2_specific_hetmers);

               // calculate pairwise comparision
           //}
       //}
   //}

   // analyses comparing three or more populations
   //if num_pops > 2{
       //println!("{}", ":D");
   //}

   println!("{}", "Done!");
}
