#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Kmer(u64);

fn encode_base(b: u8) -> Option<u64> {
    // only [atgcATGC], no [Nn]
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

impl Kmer {
    /// k-mer hashing, mm
    pub fn from_bytes(seq: &[u8]) -> Option<Kmer> {
        if seq.is_empty() || seq.len() > 31 {
            None
        } else {
            let mut code = 0;
            for &b in seq {
                // bitwise OR - left shift one base (2 bits), and then add the new encoded base with the OR
                code = (code << 2) | encode_base(b)?;
            }
            Some(Kmer(code))
        }
    }

    /// Reverse complement a Kmer
    pub fn revcomp(self, k: usize) -> Kmer {
        let mut code = 0;
        let mut x = self.0;
        for _ in 0..k {
            // hurts my brain
            let base = x & 3;
            code = (code << 2) | (3 - base);
            x = x >> 2;
        }
        Kmer(code)
    }

    pub fn canonical(self, k: usize) -> Kmer {
        let rc = self.revcomp(k);
        if rc > self { self } else { rc }
    }

    /// binary back to nucleotide
    pub fn decode(self, k: usize) -> String {
        let mut kmer_string = String::new();
        for i in 0..k {
            let current_base = match self.0 >> (2 * (k - 1 - i)) & 3 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => unreachable!(),
            };
            kmer_string.push(current_base);
        }
        kmer_string
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bitwise_practice() {
        let x: u64 = 0b00_01_10; // "ACG"
        println!("{x}");
        println!("{:06b}", x); // 000110
        println!("{:06b}", x >> 4); // 000000 -> A (0)
        println!("{:06b}", (x >> 2) & 3); // 000001 -> C (1)
        println!("{:06b}", x & 3); // 000010 -> G (2)
    }

    #[test]
    fn round_trip() {
        let km = Kmer::from_bytes(b"ACGTTGCA").unwrap();
        assert_eq!(km.decode(8), "ACGTTGCA");
    }

    #[test]
    fn revcomp_basic() {
        let km = Kmer::from_bytes(b"AACG").unwrap();
        assert_eq!(km.revcomp(4).decode(4), "CGTT");
    }

    #[test]
    fn revcomp_twice() {
        let km = Kmer::from_bytes(b"GATTACA").unwrap();
        assert_eq!(km.revcomp(7).revcomp(7), km);
    }

    #[test]
    fn canonical_same_both_strands() {
        let km = Kmer::from_bytes(b"GATTACA").unwrap();
        assert_eq!(km.canonical(7), km.revcomp(7).canonical(7));
    }

    #[test]
    fn canonical_picks_smaller() {
        let km = Kmer::from_bytes(b"TTTT").unwrap();
        assert_eq!(km.canonical(4).decode(4), "AAAA");
    }

    #[test]
    fn n_is_none() {
        assert_eq!(Kmer::from_bytes(b"ACNT"), None);
    }
}
