use clap::Parser; // argument parser
use itertools::izip;

mod one_hetmers;

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
    //#[arg(short = 'y', long, num_args = 1.., value_delimiter = ' ', required = true, value_parser = clap::builder::PossibleValuesParser::new(["hetmers", "emp_freq", "bayes_freq", "fst", "dxy", "fit"]))]
    //analyses: Vec<String>,

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

    /// shape parameter for prior distribution
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    alphas: Vec<f64>,

    /// shape parameter for prior distribution
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    betas: Vec<f64>,

    /// thresholds for determining if k-mer has abnormal copy number
    #[arg(short, long, num_args = 1.., value_delimiter = ' ', required = true)]
    sigmas: Vec<f64>,
}

// bring it all together! :D
fn main() {

    // arguments from cli
    let args = Args::parse(); 

    // loop over individual populations
    //if args.analyses.contains(&"hetmers".to_string()) { 
        for (input, output, minimum, coverage, pool, alpha, beta, sigma) in izip!(&args.inputs, &args.outputs, &args.minimums, args.coverages, args.pools, args.alphas, args.betas, args.sigmas) {
            one_hetmers::kmers_to_hetmers(input, output, *minimum, args.alleles, pool, coverage, alpha, beta, sigma);
        }
    //}

   println!("{}", "Done!");
}
