// functions to calculate fst or otherwise comapre two sets of k-mers/het-mers
//mod fst {
//
//    pub fn shared_hetmers(map1: Vec<u64>, map2: Vec<u64>) -> Vec<u64>{
//        println!("Find hetmers shared between two maps...");
//        let output = map1.into_iter()
//            .filter(|k| map2.contains(k))
//            .collect();
//        return output;
//    }
//    
//    pub fn population_specific_hetmers(map1: Vec<u64>, map2: Vec<u64>) -> Vec<u64>{
//        println!("Find hetmers specific to one population...");
//        let output = map1.into_iter()
//            .filter(|k| !map2.contains(k))
//            .collect();
//        return output;
//    }
//    
//    pub fn convert_u64_to_string(vec: &Vec<u64>) -> Vec<String> {
//        vec.iter().map(|num| num.to_string()).collect()
//    }
//    
//}

//fn verify_hetmers_mismatch(hetmers: &Vec<str>, positions: &[usize]) -> Vec<bool> {
//    println!("Checking candidate hetmers...");
//    let hetmer_checks: Vec<_> = hetmers.iter()
//        .filter_map(|s| {
//            let parts: Vec<&str> = s.split(',').collect();
//            if parts.len() == 2 {
//                let (str1, str2) = (parts[0], parts[1]);
//                return verify_muts_as_pos(str1, str2, positions);
//            }
//        None
//        })
//        .collect();

    // reformatting frequency list
//    return hetmer_checks;
//}