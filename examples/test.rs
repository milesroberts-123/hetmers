use std::f64;

fn posterior(x: i32, z: i32, n: i32, c: f64, alpha: f64, beta: f64) -> f64 {
    let mut likelihood_times_prior = Vec::new();
    let mut total_probability = 0.0;

    for i in 1..n {
        let p = i as f64 / n as f64;
        let lambdax = (i as f64) * c;
        let lambday = ((n - i) as f64) * c;
        
        let likelihood = (p.powf(x as f64 + alpha - 1.0)) * ((1.0 - p).powf((z - x) as f64 + beta - 1.0));
        let prior = 1.0 / ((lambdax.exp() - 1.0) * (lambday.exp() - 1.0));
        let value = likelihood * prior;
        
        likelihood_times_prior.push(value);
        total_probability += value;
    }
    
    let post: Vec<f64> = likelihood_times_prior.iter().map(|&v| v / total_probability).collect();
    
    if let Some((max_index, _)) = post.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()) {
        return (max_index + 1) as f64 / n as f64;
    }
    
    0.0
}

fn main() {
    let result = posterior(3, 7, 10, 1.5, 2.0, 3.0);
    println!("Result: {}", result);
}
