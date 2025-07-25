use sha2::{Digest, Sha256};
use statrs::function::gamma::gamma;
use std::f64;
use std::collections::HashSet;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

mod input_checkers;

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
            println!("Skipping line that does not have two tab-separated columns:");
            println!("{}", line);
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
fn extract_border(seqs: &Vec<String>) -> Vec<String> {
    let k = seqs[0].len();
    println!("k is {}", k);
    let k_half = k / 2;

    let no_center_ks: Vec<String> = seqs.iter()
        .map(|s| format!("{}{}", &s[..k_half], &s[k_half + 1..]))
        .collect();

    return no_center_ks;
}

// reverse complement sequences
fn rev_comp(seqs: &Vec<String>) -> Vec<String> {
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
    let output = input.into_iter()
        .filter(|(_, v)| v.len() == alleles)
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
fn write_file(output: &Vec<String>, prefix: &String, suffix: &str) {
    println!("{}", format!("Saving results to {}_{}...", prefix, suffix));
    let mut file = File::create(format!("{}_{}", prefix, suffix)).expect("Unable to create file");
    writeln!(file, "{}", output.join("\n")).expect("Unable to write to file");
}

// empirical allele frequencies
fn counts_to_frequencies(count_pairs: &Vec<String>) -> Vec<String> {
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

    // reformatting frequency list
    let freq_strings: Vec<String> = frequencies.into_iter().map(|s| s.to_string()).collect();
    return freq_strings;
}

// Tag hetmers with really high total coverage (potentially due to paralogous sequences)
fn high_cov_hetmers(count_pairs: &Vec<String>, sigma: f64, n: i32, cov: f64) -> Vec<String> {
    println!("Checking for questionable hetmers...");
    let potential_filter: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let sum = num1 + num2;
                    let stderr = ((n as f64)*(cov as f64)).sqrt();
                    if sum > (sigma*stderr) as f64 {
                        return Some(1);
                    } else{
                        return Some(0);
                    }
                }
            }
        None
        })
        .collect();

    // reformatting frequency list
    let potent_strings: Vec<String> = potential_filter.into_iter().map(|s| s.to_string()).collect();
    return potent_strings;

}

// Compute truncation constant
fn truncation_constant(c: usize, lambda: f64) -> f64 {
    let sum: f64 = (0..c)
        .map(|x| {
            let numerator = (-lambda).exp() * lambda.powi(x as i32);
            let denominator = gamma((x + 1) as f64); // factorial(x)
            numerator / denominator
        })
        .sum();

    1.0 - sum
}

// Compute posterior where the minor k-mer count is variable
fn posterior_min_kmer_count(x: f64, z: f64, n: i32, cov: f64, c: usize, alpha: f64, beta: f64) -> usize {
    let mut likelihood_times_prior = Vec::new();
    let mut total_probability = 0.0;
    let max_minor_count = n/2; // minor allele can't have frequency above 1/2 by definition

    for i in 1..max_minor_count {
        let p = i as f64 / n as f64;
        let lambda_x = (i as f64) * (cov as f64);
        let lambda_y = ((n - i) as f64) * (cov as f64);

        let tx = truncation_constant(c, lambda_x);
        let ty = truncation_constant(c, lambda_y);

        let likelihood = p.powf(x as f64 + alpha - 1.0) * (1.0 - p).powf((z - x) as f64 + beta - 1.0);
        let likelihood_truncated = likelihood / (tx * ty);
        likelihood_times_prior.push(likelihood_truncated);
        total_probability += likelihood;
    }

    let posterior: Vec<f64> = likelihood_times_prior
        .iter()
        .map(|val| val / total_probability)
        .collect();

    let max_index = posterior
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(idx, _)| idx)
        .unwrap();

    max_index + c
}

// calculate posterior distribution for allele count
fn counts_to_bayes_state(count_pairs: &Vec<String>, n: i32, cov:f64, c: usize, alpha: f64, beta: f64) -> Vec<usize>{
    println!("Calculating posterior...");
    let bayes_states: Vec<_> = count_pairs.iter()
        .filter_map(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                if let (Ok(num1), Ok(num2)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    let x = num1.min(num2);
                    //let max_num = num1.max(num2);
                    //let z = x + max_num;
                    let z = num1 + num2;
                    let post = posterior_min_kmer_count(x,z,n,cov,c,alpha,beta);
                    return Some(post);
                    //return highest_prob_index(post);
                } else { return Some(0);}
            } else { return Some(0); }
        })
        .collect();
    return bayes_states;
}

// Collect functions into one 
pub fn kmers_to_hetmers(input: &String, output: &String, minimum: usize, alleles: usize, pool: i32, coverage: f64, alpha: f64, beta: f64, sigma: f64) -> (Vec<std::string::String>, Vec<std::string::String>, Vec<u64>) {
    // load k-mers
    let kmers = load_kmers(input, minimum);

    // input checks
    input_checkers::all_checks(&kmers.0);

    // remove central base from each k-mer
    let borders = extract_border(&kmers.0);

    // reverse complement borders
    let revborders = rev_comp(&borders);

    // get hash of borders
    let hashbord = hash_seqs(borders);
    let hashrevbord = hash_seqs(revborders);

    // compare forward and reverse hash and take the min
    let min_hashes = min_hash(hashbord, hashrevbord);

    // group the hashes into a dictionary
    let grouped_hashes = group_hashes(min_hashes);

    // remove hashes that had only one or more than 2 k-mers per group
    let filtered_groups = filter_groups(grouped_hashes, alleles);

    // extract sequences for each hash group
    let hetmers = extract_hetmers(filtered_groups, kmers.0, kmers.1);

    // empirical frequencies
    let empirical_frequencies = counts_to_frequencies(&hetmers.1);

    // bayesian allele states
    let bayes_states = counts_to_bayes_state(&hetmers.1, pool, coverage, minimum, alpha, beta);

    // check for hetmers with weirdly high coverage
    let check_these_hetmers = high_cov_hetmers(&hetmers.1, sigma, pool, coverage);

    // write output files
    write_file(&hetmers.0, output, "seqs.csv");
    write_file(&hetmers.1, output, "counts.csv");
    write_file(&hetmers.2.iter().map(|num| num.to_string()).collect(), output, "hashes.csv");
    write_file(&bayes_states.into_iter().map(|s| s.to_string()).collect(), output, "bayes_states.csv");
    write_file(&empirical_frequencies, output, "empirical_freqs.csv");
    write_file(&check_these_hetmers, output, "bad_hetmers.csv");

    return hetmers;
}


//
//
//
// TRASH FUNCTIONS
//
//

// functions to process hetmers with 2 base differences
fn drop_bases(s: String, positions: &Vec<usize>) -> String {

    let filtered_s: String = s
        .chars()
        .enumerate() // gives us (index, char) tuples
        .filter(|(index, _)| !positions.contains(index))
        .map(|(_, char)| char)
        .collect();
            
    filtered_s
}

fn drop_bases_all_kmers(strs: &Vec<String>, positions: &Vec<usize>) -> Vec<String> {
    strs.iter().map(|s| drop_bases(s.to_string(), positions)).collect()
}
    

    // check if the k-mers in a hetmer actually differ
    fn verify_muts_at_pos(strs: &Vec<String>, positions: &[usize]) -> bool {

        let byte_strs: Vec<&[u8]> = strs.iter().map(|s| s.as_bytes()).collect();

        positions.iter().all(|&pos| {
            let mut seen = std::collections::HashSet::new();
            for b in byte_strs.iter() {
                seen.insert(b[pos]);
            }
            seen.len() == byte_strs.len() // All must differ
        })
    }
    
    
// test functions
#[cfg(test)]
mod units {
    use super::*;

    #[test]
    fn odd_k_borders() {
        let test_vec = vec!["ATGCA".to_string(), "TTGAT".to_string(), "GGATA".to_string()];
        let result = extract_border(&test_vec);
        let expected = vec!["ATCA".to_string(), "TTAT".to_string(), "GGTA".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn even_k_borders() {
        let test_vec = vec!["ATGCAT".to_string(), "TTGATC".to_string(), "GGATAA".to_string()];
        let result = extract_border(&test_vec);
        let expected = vec!["ATGAT".to_string(), "TTGTC".to_string(), "GGAAA".to_string()];
        assert_eq!(result, expected);
    }

    #[test]
    fn short_rev_comp() {
        let test_vec = vec!["ATGCAT".to_string(), "TTGATC".to_string(), "GGATAA".to_string()];
        let result = rev_comp(&test_vec);
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

    #[test]
    fn a_few_freqs(){
        // Vec<String>, output: String
        let count_pairs = vec!["5,120".to_string(), "1,9".to_string(), "20,140".to_string(), "22,22".to_string()];
        let result = counts_to_frequencies(&count_pairs);
        let expected = vec!["0.04".to_string(), "0.1".to_string(), "0.125".to_string(), "0.5".to_string()];
        assert_eq!(result, expected);
    }
    
    #[test]
    fn check_single_muts(){
        let strs = vec!["GATTACA".to_string(),"GATCACA".to_string()];
        let positions = vec![3];
        let result = verify_muts_at_pos(&strs, &positions);
        let expected = true;
        assert_eq!(result, expected);
    }

    #[test]
    fn check_multi_muts(){
        let strs = vec!["GATTACA".to_string(), "GCTCACA".to_string()];
        let positions = vec![1,3];
        let result = verify_muts_at_pos(&strs, &positions);
        let expected = true;
        assert_eq!(result, expected);
    }
    
    #[test]
    fn check_base_drop(){
        let result = drop_bases("GATTACA".to_string(), &vec![1,3]);
        let expected = "GTACA".to_string();
        assert_eq!(result, expected);
    }

    #[test]    
    fn check_vec_base_drop(){
        let strs = vec!["GATTACA".to_string(), "GCTCACA".to_string()];
        let positions = vec![1,3];
        let result = drop_bases_all_kmers(&strs, &positions);
        let expected = vec!["GTACA".to_string(), "GTACA".to_string()];
        assert_eq!(result, expected);
    }

    //#[test]
    //fn check_hetmer_mismatch(){
    //    let a = vec!["GATTACA,GATCACA", "GATTACA,GATAACA"];
    //    let positions = vec![3];
    //    let result = verify_hetmers_mismatch(a, positions);
    //    let expected = vec![true, true];
    //    assert_eq!(result, expected);
    //}

    // helper function to test equality of floating point numbers
    fn round_to_decimals(x: f64, decimals: u32) -> f64 {
        let factor = 10f64.powi(decimals as i32);
        (x * factor).round() / factor
    }

    #[test]
    fn small_truncation(){
        let result = round_to_decimals(truncation_constant(5, 10.0), 7);
        let expected = 0.9707473;
        assert_eq!(result, expected);
    }

    #[test]
    fn seqs_from_hashmap() {
        let mut input = HashMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec!["TCGTC".to_string(), "AATAA".to_string(), "AAGAA".to_string(), "GATGA".to_string(), "TCATC".to_string(), "GAAGA".to_string()];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.0.into_iter().collect();
        let expected: HashSet<_> = vec!["TCGTC,TCATC".to_string(), "GATGA,GAAGA".to_string(), "AATAA,AAGAA".to_string()].into_iter().collect();

        //let expected = (extracted_seqs, vec!["1,6".to_string(), "2,100".to_string(), "10,9".to_string()], vec![9875, 2222, 1111]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn counts_from_hashmap(){
        let mut input = HashMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec!["TCGTC".to_string(), "AATAA".to_string(), "AAGAA".to_string(), "GATGA".to_string(), "TCATC".to_string(), "GAAGA".to_string()];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.1.into_iter().collect();
        let expected: HashSet<_> = vec!["1,6".to_string(), "2,100".to_string(), "10,9".to_string()].into_iter().collect();

        assert_eq!(actual, expected);

    }

    #[test]
    fn hashes_from_hashmap(){
        let mut input = HashMap::new();
        input.insert(9875, vec![0, 4]);
        input.insert(1111, vec![1, 2]);
        input.insert(2222, vec![3, 5]);

        let seqs = vec!["TCGTC".to_string(), "AATAA".to_string(), "AAGAA".to_string(), "GATGA".to_string(), "TCATC".to_string(), "GAAGA".to_string()];
        let counts = vec![1, 10, 9, 2, 6, 100];

        let result = extract_hetmers(input, seqs, counts);

        let actual: HashSet<_> = result.2.into_iter().collect();
        let expected: HashSet<_> = vec![9875, 2222, 1111].into_iter().collect();

        assert_eq!(actual, expected);

    }
    
    #[test]
    fn two_alleles() {
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]); // should be filtered out
        input.insert(3, vec![2, 3, 4, 5]); // should be filtered out
	input.insert(4, vec![2]); // should be filtered out
	input.insert(5, vec![5, 6]);

        let alleles = 2;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(1, vec![0, 1]);
        expected.insert(5, vec![5, 6]);

        assert_eq!(result, expected);
    }

    #[test]
    fn three_alleles() {
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 3;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(2, vec![2, 3, 4]);

        assert_eq!(result, expected);
    }

    #[test]
    fn four_alleles() {
        let mut input = HashMap::new();
        input.insert(1, vec![0, 1]);
        input.insert(2, vec![2, 3, 4]);
        input.insert(3, vec![2, 3, 4, 5]);
        input.insert(4, vec![2]);
        input.insert(5, vec![5, 6]);

        let alleles = 4;
        let result = filter_groups(input, alleles);

        let mut expected = HashMap::new();
        expected.insert(3, vec![2, 3, 4, 5]);

        assert_eq!(result, expected);
    }

    #[test]
    fn filter_groups_no_matches() {
        let mut input = HashMap::new();
        input.insert(1, vec![0]);
        input.insert(2, vec![1, 2, 3]);

        let result = filter_groups(input, 2);
        let expected: HashMap<u64, Vec<usize>> = HashMap::new();

        assert_eq!(result, expected);
    }
}