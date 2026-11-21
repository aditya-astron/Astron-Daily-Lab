// Astron Built this AP
// Project 051 | Telemetry framing / Error correction | CCSDS telecommand BCH(63,56) codec and CLTU link
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST

//! CCSDS 231.0-B telecommand synchronization and channel coding.
//!
//! Bit convention: a 63-bit codeword lives in the low 63 bits of a u64, bit i
//! holding the coefficient of x^i. The first bit on the wire is x^62, which is
//! also the MSB of the first byte of the codeblock.

// ---------------------------------------------------------------------------
// Polynomial arithmetic over GF(2) and generic cyclic codes
// ---------------------------------------------------------------------------

/// Remainder of v(x) mod g(x), with deg g = r.
pub fn poly_mod(mut v: u64, g: u64, r: u32) -> u64 {
    while v >> r != 0 {
        let top = 63 - v.leading_zeros();
        v ^= g << (top - r);
    }
    v
}

pub fn poly_mul(a: u64, b: u64) -> u64 {
    let mut p = 0;
    for i in 0..64 {
        if b >> i & 1 == 1 {
            p ^= a << i;
        }
    }
    p
}

/// Binary cyclic (n, n-r) code with generator g of degree r.
#[derive(Clone, Copy, Debug)]
pub struct CyclicCode {
    pub n: u32,
    pub r: u32,
    pub g: u64,
}

impl CyclicCode {
    pub fn k(&self) -> u32 {
        self.n - self.r
    }

    /// Systematic encoding: c(x) = i(x) x^r + (i(x) x^r mod g).
    pub fn encode(&self, info: u64) -> u64 {
        let shifted = info << self.r;
        shifted | poly_mod(shifted, self.g, self.r)
    }

    pub fn syndrome(&self, word: u64) -> u64 {
        poly_mod(word, self.g, self.r)
    }

    /// Column i of the parity-check matrix: x^i mod g.
    pub fn h_column(&self, i: u32) -> u64 {
        poly_mod(1u64 << i, self.g, self.r)
    }

    /// The 2^r words of the dual code. Word a has bit i = <a, x^i mod g>,
    /// so <u_a, t> = <a, syndrome(t)>, which is what the coset formula uses.
    pub fn dual_words(&self) -> Vec<u64> {
        let cols: Vec<u64> = (0..self.n).map(|i| self.h_column(i)).collect();
        (0..1u64 << self.r)
            .map(|a| cols.iter().enumerate().fold(0u64, |u, (i, &c)| u | (((a & c).count_ones() as u64) & 1) << i))
            .collect()
    }
}

/// The CCSDS TC code: g(x) = x^7 + x^6 + x^2 + 1 = (x + 1)(x^6 + x + 1).
/// x^6 + x + 1 is primitive, so this is the (63,57) Hamming code expurgated to
/// its even-weight half: d = 4, single-error correcting, double detecting.
pub const TC_CODE: CyclicCode = CyclicCode { n: 63, r: 7, g: 0b1100_0101 };

// ---------------------------------------------------------------------------
// Weight enumerators (exact integer arithmetic)
// ---------------------------------------------------------------------------

pub fn binom(n: i64, k: i64) -> i128 {
    if k < 0 || k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut c: i128 = 1;
    for j in 0..k {
        c = c * (n - j) as i128 / (j + 1) as i128;
    }
    c
}

/// Krawtchouk polynomial K_w(j; n) = sum_s (-1)^s C(j, s) C(n - j, w - s).
pub fn krawtchouk(w: i64, j: i64, n: i64) -> i128 {
    (0..=w).map(|s| if s % 2 == 0 { 1 } else { -1 } * binom(j, s) * binom(n - j, w - s)).sum()
}

/// Weight distribution of the coset t + C, from the dual code alone:
///   B_w = (1/|C_perp|) sum_{u in C_perp} (-1)^{<u,t>} K_w(wt u).
/// t = 0 gives MacWilliams' identity for the code itself.
pub fn coset_weight_distribution(code: &CyclicCode, t: u64) -> Vec<i128> {
    let n = code.n as i64;
    let dual = code.dual_words();
    let mut b = vec![0i128; code.n as usize + 1];
    for (w, bw) in b.iter_mut().enumerate() {
        let mut acc: i128 = 0;
        for &u in &dual {
            let sign = if (u & t).count_ones().is_multiple_of(2) { 1 } else { -1 };
            acc += sign * krawtchouk(w as i64, u.count_ones() as i64, n);
        }
        assert!(acc % dual.len() as i128 == 0, "MacWilliams sum not divisible");
        *bw = acc / dual.len() as i128;
    }
    b
}

/// Closed form for the even-weight subcode of the length n = 2^m - 1 Hamming
/// code: even part of A(z) = [(1+z)^n + n (1-z)(1-z^2)^((n-1)/2)] / (n+1).
pub fn expurgated_hamming_distribution(n: i64) -> Vec<i128> {
    let h = (n - 1) / 2;
    let mut a = vec![0i128; n as usize + 1];
    for w in 0..=n {
        let mut num = binom(n, w);
        // n (1 - z)(1 - z^2)^h: coefficient of z^w.
        let even = |e: i64| if e % 2 == 0 { (if (e / 2) % 2 == 0 { 1 } else { -1 }) * binom(h, e / 2) } else { 0 };
        num += n as i128 * (even(w) - if w >= 1 { even(w - 1) } else { 0 });
        assert!(num % (n + 1) as i128 == 0);
        a[w as usize] = num / (n + 1) as i128;
    }
    for w in (1..=n as usize).step_by(2) {
        a[w] = 0;
    }
    a
}

/// Weight-4 codewords counted from parity-check columns: a weight-4 word is a
/// 4-set of columns XORing to zero, i.e. two disjoint pairs with equal XOR.
/// Each 4-set splits into pairs three ways.
pub fn count_weight4_from_columns(code: &CyclicCode) -> i128 {
    let cols: Vec<u64> = (0..code.n).map(|i| code.h_column(i)).collect();
    let mut hist = vec![0i128; 1 << code.r];
    for i in 0..cols.len() {
        for j in i + 1..cols.len() {
            hist[(cols[i] ^ cols[j]) as usize] += 1;
        }
    }
    hist.iter().map(|&c| c * (c - 1) / 2).sum::<i128>() / 3
}

// ---------------------------------------------------------------------------
// Codeblock encoder / decoder
// ---------------------------------------------------------------------------

const PARITY_MASK: u64 = 0x7F;

/// Bit-serial LFSR, as in the encoder figure of 231.0-B: shift each info bit
/// in MSB first, feedback = input XOR register MSB.
pub fn parity_lfsr(info: &[u8; 7]) -> u8 {
    let mut reg: u8 = 0;
    for &byte in info {
        for k in (0..8).rev() {
            let fb = (byte >> k & 1) ^ (reg >> 6 & 1);
            reg = (reg << 1) & 0x7F;
            if fb == 1 {
                reg ^= 0x45; // g(x) without the x^7 term
            }
        }
    }
    reg
}

/// Byte-at-a-time encoder: the 7-bit register is kept left-aligned in a byte,
/// which turns it into an ordinary MSB-first CRC-8 with polynomial 0x8A.
pub struct TableEncoder {
    table: [u8; 256],
}

impl Default for TableEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl TableEncoder {
    pub fn new() -> Self {
        let mut table = [0u8; 256];
        for (i, t) in table.iter_mut().enumerate() {
            let mut c = i as u8;
            for _ in 0..8 {
                c = if c & 0x80 != 0 { (c << 1) ^ 0x8A } else { c << 1 };
            }
            *t = c;
        }
        Self { table }
    }

    pub fn parity(&self, info: &[u8; 7]) -> u8 {
        info.iter().fold(0u8, |c, &b| self.table[(c ^ b) as usize]) >> 1
    }

    /// 8-octet codeblock: 7 info octets, then the complemented parity bits and
    /// a zero filler bit. Complementing keeps an all-zero (or idle) input from
    /// producing an all-zero codeblock, which helps bit-transition density.
    pub fn codeblock(&self, info: &[u8; 7]) -> [u8; 8] {
        let mut cb = [0u8; 8];
        cb[..7].copy_from_slice(info);
        cb[7] = (!self.parity(info) & 0x7F) << 1;
        cb
    }
}

/// Undo the parity complement and drop the filler bit: the 63-bit word.
pub fn block_to_word(cb: &[u8; 8]) -> u64 {
    let hi = cb[..7].iter().fold(0u64, |w, &b| w << 8 | b as u64);
    hi << 7 | ((!(cb[7] >> 1) as u64) & PARITY_MASK)
}

pub fn word_to_block(w: u64) -> [u8; 8] {
    let mut cb = [0u8; 8];
    for (i, byte) in cb[..7].iter_mut().enumerate() {
        *byte = (w >> (7 + 8 * (6 - i))) as u8;
    }
    cb[7] = (!(w as u8) & 0x7F) << 1;
    cb
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecodeMode {
    /// Triple-error detection: accept only a zero syndrome.
    Detect,
    /// Single-error correction (also detects all double errors).
    Correct,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BlockResult {
    Clean([u8; 7]),
    Corrected([u8; 7], u32),
    Rejected,
}

pub struct Decoder {
    /// Syndrome -> error position + 1 for the 63 single-bit patterns, else 0.
    single: [u8; 128],
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    pub fn new() -> Self {
        let mut single = [0u8; 128];
        for i in 0..TC_CODE.n {
            let s = TC_CODE.h_column(i) as usize;
            assert_eq!(single[s], 0, "x^i mod g must be distinct");
            single[s] = (i + 1) as u8;
        }
        Self { single }
    }

    pub fn decode(&self, cb: &[u8; 8], mode: DecodeMode) -> BlockResult {
        let w = block_to_word(cb);
        let s = TC_CODE.syndrome(w) as usize;
        let info = |w: u64| {
            let mut out = [0u8; 7];
            out.copy_from_slice(&word_to_block(w)[..7]);
            out
        };
        if s == 0 {
            return BlockResult::Clean(info(w));
        }
        if mode == DecodeMode::Correct && self.single[s] != 0 {
            let pos = self.single[s] as u32 - 1;
            return BlockResult::Corrected(info(w ^ 1u64 << pos), pos);
        }
        BlockResult::Rejected
    }
}

// ---------------------------------------------------------------------------
// Exact codeblock statistics on a binary symmetric channel
// ---------------------------------------------------------------------------

/// Probability that the error pattern lands within Hamming distance `radius`
/// (0 or 1) of some word of a set with weight distribution `b`. Spheres of
/// radius 1 around words at mutual distance >= 4 never overlap, so the terms
/// just add: a word of weight w is reached by w-1, w or w+1 errors.
pub fn sphere_hit_probability(b: &[i128], p: f64, radius: u32) -> f64 {
    let n = (b.len() - 1) as i32;
    let q = 1.0 - p;
    let pr = |e: i32| if e < 0 || e > n { 0.0 } else { p.powi(e) * q.powi(n - e) };
    b.iter()
        .enumerate()
        .map(|(w, &bw)| {
            let w = w as i32;
            let mut t = pr(w);
            if radius == 1 {
                t += w as f64 * pr(w - 1) + (n - w) as f64 * pr(w + 1);
            }
            bw as f64 * t
        })
        .sum()
}

#[derive(Clone, Copy, Debug)]
pub struct BlockStats {
    pub correct: f64,
    pub undetected: f64,
    pub rejected: f64,
}

/// Exact outcome probabilities for one codeblock with i.i.d. bit errors p.
pub fn block_stats(a: &[i128], p: f64, mode: DecodeMode) -> BlockStats {
    let mut nonzero = a.to_vec();
    nonzero[0] = 0;
    let radius = if mode == DecodeMode::Correct { 1 } else { 0 };
    let n = (a.len() - 1) as i32;
    let q = 1.0 - p;
    let correct = q.powi(n) + if radius == 1 { n as f64 * p * q.powi(n - 1) } else { 0.0 };
    let undetected = sphere_hit_probability(&nonzero, p, radius);
    BlockStats { correct, undetected, rejected: 1.0 - correct - undetected }
}

/// Codewords whose leading `zero_info_bits` information bits are zero form
/// the shortened cyclic code of length n - zero_info_bits with the same g. Its
/// dual is the TC dual punctured to the low positions, so MacWilliams applies
/// unchanged. Returned padded to length n + 1 for use with n-bit error
/// patterns.
pub fn shortened_distribution(code: &CyclicCode, zero_info_bits: u32) -> Vec<i128> {
    let short = CyclicCode { n: code.n - zero_info_bits, ..*code };
    let mut a = coset_weight_distribution(&short, 0);
    a.resize(code.n as usize + 1, 0);
    a
}

/// Probability that a codeblock decodes (SEC) to information that agrees with
/// what was sent in its first `significant_bits` bits. In the last block of a
/// CLTU the rest is fill, so a miscorrection confined to the fill is harmless.
pub fn last_block_success(code: &CyclicCode, p: f64, significant_bits: u32, mode: DecodeMode) -> f64 {
    let a = shortened_distribution(code, significant_bits);
    sphere_hit_probability(&a, p, if mode == DecodeMode::Correct { 1 } else { 0 })
}

// ---------------------------------------------------------------------------
// TC transfer frame, randomizer, CLTU
// ---------------------------------------------------------------------------

/// CRC-16-CCITT, poly 0x1021, preset 0xFFFF, no reflection (the TC FECF).
pub fn crc16_ccitt(data: &[u8]) -> u16 {
    let mut crc: u16 = 0xFFFF;
    for &b in data {
        crc ^= (b as u16) << 8;
        for _ in 0..8 {
            crc = if crc & 0x8000 != 0 { (crc << 1) ^ 0x1021 } else { crc << 1 };
        }
    }
    crc
}

/// TC pseudo-randomizer: h(x) = x^8 + x^6 + x^4 + x^3 + x^2 + x + 1, all-ones
/// preset, period 255. Applied over the whole frame before BCH encoding.
pub fn randomizer_sequence(len: usize) -> Vec<u8> {
    let mut reg: u8 = 0xFF;
    (0..len)
        .map(|_| {
            let mut byte = 0u8;
            for _ in 0..8 {
                byte = byte << 1 | (reg & 1);
                let fb = (reg ^ reg >> 1 ^ reg >> 2 ^ reg >> 3 ^ reg >> 4 ^ reg >> 6) & 1;
                reg = reg >> 1 | fb << 7;
            }
            byte
        })
        .collect()
}

pub fn randomize(data: &mut [u8]) {
    let seq = randomizer_sequence(data.len());
    for (d, r) in data.iter_mut().zip(seq) {
        *d ^= r;
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TcFrame {
    pub bypass: bool,
    pub scid: u16,
    pub vcid: u8,
    pub seq: u8,
    pub data: Vec<u8>,
}

impl TcFrame {
    /// 5-octet primary header, data field, 2-octet FECF.
    pub fn to_bytes(&self) -> Vec<u8> {
        let total = 5 + self.data.len() + 2;
        assert!(total <= 1024 && self.scid < 1024 && self.vcid < 64);
        let len = (total - 1) as u16;
        let mut f = vec![
            (self.bypass as u8) << 5 | (self.scid >> 8) as u8 & 0x3,
            self.scid as u8,
            self.vcid << 2 | (len >> 8) as u8 & 0x3,
            len as u8,
            self.seq,
        ];
        f.extend_from_slice(&self.data);
        let crc = crc16_ccitt(&f);
        f.extend_from_slice(&crc.to_be_bytes());
        f
    }

    pub fn from_bytes(buf: &[u8]) -> Result<Self, FrameError> {
        if buf.len() < 7 {
            return Err(FrameError::Short);
        }
        let total = (((buf[2] & 0x3) as usize) << 8 | buf[3] as usize) + 1;
        if buf[0] >> 6 != 0 || total < 7 || total > buf.len() {
            return Err(FrameError::Short);
        }
        if crc16_ccitt(&buf[..total - 2]) != u16::from_be_bytes([buf[total - 2], buf[total - 1]]) {
            return Err(FrameError::Crc);
        }
        Ok(TcFrame {
            bypass: buf[0] >> 5 & 1 == 1,
            scid: ((buf[0] & 0x3) as u16) << 8 | buf[1] as u16,
            vcid: buf[2] >> 2,
            seq: buf[4],
            data: buf[5..total - 2].to_vec(),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameError {
    NoStart,
    Short,
    Crc,
}

pub const START_SEQ: u16 = 0xEB90;
pub const TAIL_SEQ: [u8; 8] = [0xC5, 0xC5, 0xC5, 0xC5, 0xC5, 0xC5, 0xC5, 0x79];
pub const FILL: u8 = 0x55;

/// CLTU = start sequence, codeblocks (last one padded with 0x55), tail.
pub fn build_cltu(frame: &[u8], enc: &TableEncoder) -> Vec<u8> {
    let mut data = frame.to_vec();
    randomize(&mut data);
    while !data.len().is_multiple_of(7) {
        data.push(FILL);
    }
    let mut out = START_SEQ.to_be_bytes().to_vec();
    for chunk in data.chunks(7) {
        out.extend_from_slice(&enc.codeblock(chunk.try_into().unwrap()));
    }
    out.extend_from_slice(&TAIL_SEQ);
    out
}

pub fn bytes_to_bits(b: &[u8]) -> Vec<u8> {
    b.iter().flat_map(|&x| (0..8).rev().map(move |k| x >> k & 1)).collect()
}

pub fn bits_to_bytes(bits: &[u8]) -> Vec<u8> {
    bits.chunks(8).map(|c| c.iter().fold(0u8, |a, &b| a << 1 | b)).collect()
}

#[derive(Clone, Debug)]
pub struct CltuReception {
    pub lock_bit: Option<usize>,
    pub blocks_accepted: usize,
    pub corrections: usize,
    pub frame: Result<TcFrame, FrameError>,
}

/// Receiver: bit-wise search for the start sequence with up to `start_tol`
/// errors, then decode 64-bit codeblocks until one is rejected (normally the
/// tail). If a lock yields no accepted codeblock it was false; search resumes.
pub fn receive_cltu(bits: &[u8], dec: &Decoder, mode: DecodeMode, start_tol: u32) -> CltuReception {
    let mut pos = 0;
    while pos + 16 <= bits.len() {
        let win = bits[pos..pos + 16].iter().fold(0u16, |a, &b| a << 1 | b as u16);
        if (win ^ START_SEQ).count_ones() > start_tol {
            pos += 1;
            continue;
        }
        let mut data = Vec::new();
        let mut at = pos + 16;
        let mut corrections = 0;
        while at + 64 <= bits.len() {
            let cb: [u8; 8] = bits_to_bytes(&bits[at..at + 64]).try_into().unwrap();
            match dec.decode(&cb, mode) {
                BlockResult::Clean(d) => data.extend_from_slice(&d),
                BlockResult::Corrected(d, _) => {
                    data.extend_from_slice(&d);
                    corrections += 1;
                }
                BlockResult::Rejected => break,
            }
            at += 64;
        }
        if data.is_empty() {
            pos += 1;
            continue;
        }
        let blocks = data.len() / 7;
        randomize(&mut data);
        return CltuReception { lock_bit: Some(pos), blocks_accepted: blocks, corrections, frame: TcFrame::from_bytes(&data) };
    }
    CltuReception { lock_bit: None, blocks_accepted: 0, corrections: 0, frame: Err(FrameError::NoStart) }
}

// ---------------------------------------------------------------------------
// Small deterministic RNG (xoshiro256**) for the Monte Carlo runs
// ---------------------------------------------------------------------------

pub struct Rng([u64; 4]);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ x >> 30).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ x >> 27).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ x >> 31
        };
        Rng([next(), next(), next(), next()])
    }

    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let out = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        out
    }

    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// n-bit i.i.d. Bernoulli(p) pattern via geometric gaps between errors.
    pub fn error_pattern(&mut self, n: u32, p: f64) -> u64 {
        let mut e = 0u64;
        let lq = (1.0 - p).ln();
        let mut i: f64 = -1.0;
        loop {
            i += 1.0 + ((1.0 - self.uniform()).ln() / lq).floor();
            if i >= n as f64 {
                return e;
            }
            e |= 1u64 << i as u32;
        }
    }
}

// ---------------------------------------------------------------------------
// Monte Carlo of the codeblock decoder
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default)]
pub struct BlockCounts {
    pub trials: u64,
    pub correct: u64,
    pub undetected: u64,
    pub rejected: u64,
}

pub fn simulate_blocks(p: f64, mode: DecodeMode, trials: u64, seed: u64) -> BlockCounts {
    let enc = TableEncoder::new();
    let dec = Decoder::new();
    let mut rng = Rng::new(seed);
    let mut c = BlockCounts { trials, ..Default::default() };
    for _ in 0..trials {
        let info: [u8; 7] = rng.next_u64().to_be_bytes()[..7].try_into().unwrap();
        let mut cb = enc.codeblock(&info);
        // Errors on the 63 coded bits; the filler bit is ignored by the decoder.
        let e = rng.error_pattern(63, p);
        let w = block_to_word(&cb) ^ e;
        cb = word_to_block(w);
        match dec.decode(&cb, mode) {
            BlockResult::Clean(d) | BlockResult::Corrected(d, _) if d == info => c.correct += 1,
            BlockResult::Rejected => c.rejected += 1,
            _ => c.undetected += 1,
        }
    }
    c
}
