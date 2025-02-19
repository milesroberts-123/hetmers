use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use clap::Parser; // argument parser
//use bio::alphabets::dna;
use sha2::{Digest, Sha256};

// Structure of input arguments
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    // kmer count table file name
    #[arg(short, long)]
    input: String,

    // minimum kmer count
    #[arg(short, long, default_value_t = 1)]
    minimum: usize,

    // number of alleles in each hetmer
    #[arg(short = 'l', long, default_value_t = 2)]
    alleles: usize,

    // kmer count table file name
    #[arg(short, long)]
    output: String,

    // mean k-mer coverage
    #[arg(short, long)]
    coverage: f64,

    // pool size
    #[arg(short, long)]
    pool: i32,

    // minimum kmer count
    #[arg(short, long, default_value_t = 1.0)]
    alpha: f64,

    // minimum kmer count
    #[arg(short, long, default_value_t = 1.0)]
    beta: f64,
}


// Function to read in kmer count table
fn load_kmers(input: String, minimum: usize) -> (Vec<String>, Vec<usize>) {
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
fn extract_hetmers(hashdict: HashMap<u64, Vec<usize>>, seqs: Vec<String>, counts: Vec<usize>) -> (Vec<String>, Vec<String>) {
    println!("Extracting counts and sequences...");
    let hetmer_seqs: Vec<String> = hashdict.values()
        .map(|indices| indices.iter().map(|&i| seqs[i].clone()).collect::<Vec<String>>().join(","))
        .collect();
    let hetmer_counts: Vec<String> = hashdict.values()
        .map(|indices| indices.iter().map(|&i| counts[i].to_string()).collect::<Vec<String>>().join(","))
        .collect();

    return (hetmer_seqs, hetmer_counts);
}

// write vector to file
fn write_file(output: Vec<String>, prefix: String, suffix: String) {
    println!("Saving results...");
    let mut file = File::create(format!("{}_{}", prefix, suffix)).expect("Unable to create file");
    writeln!(file, "{}", output.join("\n")).expect("Unable to write to file");
}

// Collect functions into one 
fn kmers_to_hetmers(input: String, output: String, minimum: usize, alleles: usize) -> (Vec<std::string::String>, Vec<std::string::String>){
    let kmers = load_kmers(input, minimum);
    println!("{:?}", kmers);

    let borders = extract_border(kmers.0.clone());
    println!("{:?}", borders);

    let revborders = rev_comp(borders.clone());
    println!("{:?}", revborders);

    let hashbord = hash_seqs(borders);
    println!("{:?}", hashbord);

    let hashrevbord = hash_seqs(revborders);
    println!("{:?}", hashrevbord);

    let min_hashes = min_hash(hashbord, hashrevbord);
    println!("{:?}", min_hashes);

    let grouped_hashes = group_hashes(min_hashes);
    println!("{:?}", grouped_hashes);

    let filtered_groups = filter_groups(grouped_hashes, alleles);
    println!("{:?}", filtered_groups);

    let hetmers = extract_hetmers(filtered_groups, kmers.0, kmers.1);
    println!("{:?}", hetmers);

    write_file(hetmers.0.clone(), output.clone(), "seqs.csv".to_string());
    write_file(hetmers.1.clone(), output, "counts.csv".to_string());

    return hetmers;
}

// empirical allele frequencies
fn counts_to_frequencies(count_pairs: Vec<String>) -> Vec<f64> {
    println!("Calculating frequencies...");
    let frequencies = count_pairs.iter()
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

fn counts_to_bayes_state(count_pairs: Vec<String>, n: i32, c: f64, alpha: f64, beta: f64) -> Vec<usize>{
    println!("Calculating posterior...");
    let bayes_states = count_pairs.iter()
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
    return bayes_states;
}

// fst

// frequency increment test

// mix it all together! :D
fn main() {
    //let args: Vec<String> = env::args().collect();
    let args = Args::parse(); 

    //for _ in 0..args.count {
    //    println!("Hello {}!", args.name);
    //}

    //let input = &args[1];  // Input file path
    //let alleles = &args[2]; // How many alleles to have for hetmers
    //let minimum = &args[2]; // Minimum k-mer count
    //let output_prefix = &args[3];  // Output file prefix

    //println!("Input file is {input}");
    //println!("allele type is {alleles}");
    //println!("Minimum hetmer count is {minimum}");
    //println!("Output prefix is {output_prefix}");

    //let minimum = minimum.parse::<usize>().unwrap();
    //let kmers = load_kmers(args.input, args.minimum);

    //println!("{:?}", kmers);

    //let borders = extract_border(kmers.0.clone());
    //println!("{:?}", borders);

    //let revborders = rev_comp(borders.clone());
    //println!("{:?}", revborders);

    //let hashbord = hash_seqs(borders);
    //println!("{:?}", hashbord);

    //let hashrevbord = hash_seqs(revborders);
    //println!("{:?}", hashrevbord);

    //let min_hashes = min_hash(hashbord, hashrevbord);
    //println!("{:?}", min_hashes);

    //let grouped_hashes = group_hashes(min_hashes);
    //println!("{:?}", grouped_hashes);

    //let filtered_groups = filter_groups(grouped_hashes, args.alleles);
    //println!("{:?}", filtered_groups);

    //let hetmers = extract_hetmers(filtered_groups, kmers.0, kmers.1);
    //println!("{:?}", hetmers);

    //write_file(hetmers.0, args.output.clone(), "seqs.csv".to_string());
    //write_file(hetmers.1, args.output, "counts.csv".to_string());

    let hetmers = kmers_to_hetmers(args.input, args.output, args.minimum, args.alleles);

    let empirical_freqs = counts_to_frequencies(hetmers.1.clone());
    println!("{:?}", empirical_freqs);

    //println!("{:?}", posterior(100.0, 200.0, 50, 100.0, 1.0, 1.0));
    //println!("{:?}", highest_prob_index(posterior(25.0, 300.0, 50, 300.0, 1.0, 1.0)));
    
    let bayes_states = counts_to_bayes_state(hetmers.1, args.pool, args.coverage, args.alpha, args.beta);
    println!("{:?}", bayes_states);
}
