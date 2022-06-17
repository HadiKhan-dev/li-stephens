use ndarray::{Array1,Array2,Axis};
use std::cmp;


const BASES : [char;10] = ['a','c','g','t','n','A','C','G','T','N'];

/// Quick helper function for testing whether or not a character is a valid base pair
pub fn is_base(char_val : char) -> bool {
    if BASES.contains(&char_val) {
        return true;
    }
    return false;
}

/// Calculates the liklihood of a given test sequence being generated from the same population that generated the given haplotypes
/// and a rate for mutations and recombinations. The haplotypes and test sequence must be input as a matrix of characters. The lenght
/// of all the haplotypes and the test sequence must be the same and all characters in uppercase.
pub fn li_stephens_probs_matrix(haplotypes: &ndarray::Array2::<char>, test_haplotype: &ndarray::Array1::<char>, mutation_prob: f64, recombination_prob: f64) -> f64 {
    let mut max_val = haplotypes.len_of(Axis(1));

    max_val = cmp::max(max_val,test_haplotype.len());


    let new_test_array = test_haplotype.as_slice().unwrap();

    for i in 0..new_test_array.len()-1 {
        if !is_base(new_test_array[i]) {
            panic!("Not all characters in the test haplotype are a base");
        }
    }

    for i in 0..haplotypes.len_of(Axis(0)) {
        for j in 0..haplotypes.len_of(Axis(1)) {
            if !is_base(haplotypes[[i,j]]) {
                panic!("Not all characters in haplotype {} are a base",i);
            }
        }
    }   


    let mut prob_array = Array2::<f64>::zeros((haplotypes.len_of(Axis(0)),max_val));

    let char_array = haplotypes;

    let test_array = test_haplotype;

    for idx in 0..max_val {
        for inner_idx in 0..haplotypes.len_of(Axis(0)) {
            let mut transition: f64 = 1.0;
            if idx != 0 {
                transition = (1.0-recombination_prob)*prob_array[[inner_idx,idx-1]];
                for inner_idx2 in 0..haplotypes.len_of(Axis(0)) {
                    transition += (recombination_prob/(haplotypes.len() as f64))*prob_array[[inner_idx2,idx-1]];
                }
            }
            let emission = if test_array[idx] == 'N' || char_array[[inner_idx,idx]] == 'N' {
                1.0
            } else if test_array[idx] == char_array[[inner_idx,idx]] {
                1.0-mutation_prob
            } else {
                mutation_prob/3.0
            };
            prob_array[[inner_idx,idx]] = transition*emission;
        }
    }

    let mut ret_val : f64 = 0.0;

    for final_index in 0..haplotypes.len_of(Axis(0)) {
            if prob_array[[final_index,max_val-1]] > ret_val {
                ret_val = prob_array[[final_index,max_val-1]];
            }
    }

    return ret_val;



    
}



/// Calculates the liklihood of a given test sequence being generated from the same population that generated the given haplotypes
/// and a rate for mutations and recombinations. The haplotypes and test sequence must be input as a vector of strings and a string
/// respectively.
pub fn li_stephens_probs(haplotypes: &Vec<String>, test_haplotype: &String, mutation_prob: f64, recombination_prob: f64) -> f64 {

    let mut max_val = 0;

    for x in haplotypes {
        if x.len() > max_val {
            max_val = x.len();
        }
    }


    max_val = cmp::max(max_val,test_haplotype.len());

    let mut char_array = Array2::<char>::from_elem((haplotypes.len(),max_val),'N');

    let mut cur_index = 0;

    for x in haplotypes {
        for idx in 0..(x.len()) {
            char_array[[cur_index,idx]] = (x.as_bytes()[idx] as char).to_ascii_uppercase();
        }
        cur_index += 1;
    }

    let mut test_array = Array1::<char>::from_elem(max_val,'N');

    for idx in 0..(test_array.len()) {
        test_array[[idx]] = test_haplotype.as_bytes()[idx] as char;
    }



    return li_stephens_probs_matrix(&char_array,&test_array,mutation_prob,recombination_prob);

}
