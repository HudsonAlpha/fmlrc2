
#[derive(Clone,Debug,Eq,PartialEq)]
pub struct DynamicWFA<'a> {
    primary_seq: &'a [u8],
    other_seq: Vec<u8>,
    edit_distance: usize,
    distance_vector: Vec<usize>,
    max_distance: usize
}

impl<'a> Ord for DynamicWFA<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.max_distance as f64 / (self.edit_distance+1) as f64).partial_cmp(
            &(other.max_distance as f64 / (other.edit_distance+1) as f64)
        ).unwrap()
    }
}
impl<'a> PartialOrd for DynamicWFA<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> DynamicWFA<'a> {
    #[inline]
    pub fn new(primary_seq: &'a [u8]) -> DynamicWFA<'a> {
        DynamicWFA {
            primary_seq,
            other_seq: Vec::<u8>::new(),
            edit_distance: 0,
            distance_vector: vec![0],
            max_distance: 0
        }
    }

    #[inline]
    pub fn get_other_seq(&self) -> &Vec<u8> {
        &self.other_seq
    }

    #[inline]
    pub fn get_edit_distance(&self) -> usize {
        self.edit_distance
    }

    #[inline]
    pub fn get_max_distance(&self) -> usize {
        self.max_distance
    }

    #[inline]
    pub fn append(&mut self, symbol: u8) -> (usize, usize) {
        self.other_seq.push(symbol);
        let extension_made: bool = self.maximal_extend();

        if !extension_made {
            self.increase_edit_distance();
            self.maximal_extend();
        }

        self.max_distance = *self.distance_vector.iter().max().unwrap();
        
        //return the maximum extension and the current edit distance
        (self.max_distance, self.edit_distance)
    }


    #[inline]
    pub fn finalize(&mut self) -> (usize, usize) {
        while self.max_distance < self.primary_seq.len() {
            self.increase_edit_distance();
            self.maximal_extend();
            self.max_distance = *self.distance_vector.iter().max().unwrap();
        }

        //return the maximum extension and the current edit distance
        (self.max_distance, self.edit_distance)
    }

    #[inline]
    fn maximal_extend(&mut self) -> bool {
        let mut extension_made = false;
        for (i, d) in self.distance_vector.iter_mut().enumerate() {
            while *d < self.primary_seq.len() &&
                *d+self.edit_distance-i < self.other_seq.len() &&
                self.primary_seq[*d] == self.other_seq[*d+self.edit_distance-i]
            {
                *d += 1;
                extension_made = true;
            }
        }
        extension_made
    }

    #[inline]
    fn increase_edit_distance(&mut self) {
        //actually increase the edit distance
        self.edit_distance += 1;

        //create the new distance vector and populate it before overwriting the old one
        let new_dist_len = self.distance_vector.len()+2;
        let mut new_distance_vector: Vec<usize> = vec![0; new_dist_len];

        //edge cases: 
        //far left can only be a "deletion" of the other_seq, so no advance
        new_distance_vector[0] = self.distance_vector[0];
        //far right can only be a "deletion" of the primary_seq, so we increment 
        new_distance_vector[new_dist_len-1] = self.distance_vector[new_dist_len-3]+1;

        if self.edit_distance == 1 {
            //this is a special case also, the center point can only come from a mismatch
            new_distance_vector[1] = self.distance_vector[0]+1;
        } else {
            //two more special cases
            new_distance_vector[1] = std::cmp::max(
                self.distance_vector[0]+1, //mismatch
                self.distance_vector[1]     //deletion of the other_seq
                //deletion of the primary_seq is not available here
            );
            new_distance_vector[new_dist_len-2] = std::cmp::max(
                self.distance_vector[new_dist_len-3]+1, //mismatch
                //deletion of the other_seq is not available here
                self.distance_vector[new_dist_len-4]+1
            );

            //normal cases
            for i in 2..new_dist_len-2 {
                new_distance_vector[i] = std::cmp::max(
                    self.distance_vector[i-1]+1, //mismatch
                    std::cmp::max(
                        self.distance_vector[i], //deletion of the other_seq, doesn't progress us in the primary_seq
                        self.distance_vector[i-2]+1 //deletion of the primary_seq, does progress us
                    )
                )
            }
        }
        
        //save it off and send it back
        self.distance_vector = new_distance_vector;

        //println!("after edit extension \n{:?}", self);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::string_util;

    #[test]
    fn test_simple() {
        let test_seq = string_util::convert_stoi(&"AC");
        let modified_seq = string_util::convert_stoi(&"AACTG");

        //the third symbol still allows for a ED=1 (deletion)
        //after that, each one increases the "distance" (which is past the primary_seq now) and the edit_distance
        let expected: Vec<(usize, usize)> = vec![(1, 0), (2, 1), (2, 1), (3, 2), (4, 3)];

        //just test index allocation at the corners
        let mut dwfa = DynamicWFA::new(&test_seq);

        for (i, v) in modified_seq.iter().enumerate() {
            let result = dwfa.append(*v);
            assert_eq!(result, expected[i]);
        }
    }

    #[test]
    fn test_exact() {
        //tests an exact matching string
        let test_seq = string_util::convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let mut dwfa = DynamicWFA::new(&test_seq);
        for (i, v) in test_seq.iter().enumerate() {
            let result = dwfa.append(*v);
            assert_eq!(result, (i+1, 0));
        }
    }

    #[test]
    fn test_multiple() {
        let test_seq     = string_util::convert_stoi(&"AACGGATCAAGCTTACCAGTATTTACGT");
        let modified_seq = string_util::convert_stoi(&"AACGGACAAAAGCTTACCTGTATTACGT");
        let min_ed = 5; //modified_seq has 2 separate deletions, 1 2bp insertion, and 1 mismatch
        
        let mut dwfa = DynamicWFA::new(&test_seq);
        let mut result: (usize, usize) = (0, 0);
        
        for v in modified_seq.iter() {
            result = dwfa.append(*v);
        }
        assert_eq!(result, (test_seq.len(), min_ed));
    }

    #[test]
    fn test_big_insertion() {
        let test_seq     = string_util::convert_stoi(&"AACGGATTTTACGT");
        let modified_seq = string_util::convert_stoi(&"AACGGATAAAAGCTTACCTGTTTTACGT");//there's a big insert in the middle
        let min_ed = modified_seq.len() - test_seq.len();
        
        let mut dwfa = DynamicWFA::new(&test_seq);
        let mut result: (usize, usize) = (0, 0);
        
        for v in modified_seq.iter() {
            result = dwfa.append(*v);
        }
        assert!(result.0 >= test_seq.len() && result.1 == min_ed);
    }

    #[test]
    fn test_big_deletion() {
        //interesting phenomenon, the deletion must be <= right side flanking sequence to finalize correctly
        //otherwise, it think that mismatching is better, and it just hasn't mismatched enough yet
        //TODO: is there a way to make it do the whole thing correctly?
        let test_seq =     string_util::convert_stoi(&"ATTTTTTTTTTAAAAAAAAAA");
        let modified_seq = string_util::convert_stoi(&"AAAAAAAAAAA");//10 Ts are deleted
        let min_ed = test_seq.len() - modified_seq.len();
        
        let mut dwfa = DynamicWFA::new(&test_seq);
        let mut result: (usize, usize) = (0, 0);
        
        for v in modified_seq.iter() {
            result = dwfa.append(*v);
        }
        assert_eq!(result, (test_seq.len(), min_ed));
    }

    #[test]
    fn test_finalize() {
        //interesting phenomenon, the deletion must be <= right side flanking sequence to finalize correctly
        //otherwise, it think that mismatching is better, and it just hasn't mismatched enough yet
        //TODO: is there a way to make it do the whole thing correctly?
        let test_seq =     string_util::convert_stoi(&"ATTTTTTTTTTA");
        let modified_seq = string_util::convert_stoi(&"AA");//10 Ts are deleted
        let min_ed = test_seq.len() - modified_seq.len();
        
        let mut dwfa = DynamicWFA::new(&test_seq);
        for v in modified_seq.iter() {
            dwfa.append(*v);
        }
        let result = dwfa.finalize();
        assert_eq!(result, (test_seq.len(), min_ed));
    }

    #[test]
    fn test_cloning() {
        let test_seq =     string_util::convert_stoi(&"AAAAAAA");
        let modified_seq = string_util::convert_stoi(&"AAACAAA");
        
        let mut dwfa = DynamicWFA::new(&test_seq);
        let mut dwfa2: DynamicWFA = dwfa.clone();
        let mut restart_point: usize = 0;
        for (i, v) in modified_seq.iter().enumerate() {
            if test_seq[i] == *v {
                dwfa.append(*v);
            } else {
                restart_point = i;
                dwfa2 = dwfa.clone();
                break;
            }
        }

        dwfa.append(modified_seq[restart_point]);
        dwfa2.append(test_seq[restart_point]);
        println!("");
        println!("{:?}", dwfa);
        println!("{:?}", dwfa2);
        assert_ne!(dwfa, dwfa2);
        assert!(dwfa2 > dwfa);
    }
}