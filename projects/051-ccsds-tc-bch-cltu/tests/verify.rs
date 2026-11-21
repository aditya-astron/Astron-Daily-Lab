// Astron Built this AP
// Project 051 | Telemetry framing / Error correction | CCSDS telecommand BCH(63,56) codec and CLTU link
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST

use tc_cltu::*;

const P6: u64 = 0b100_0011; // x^6 + x + 1

/// Small sibling of the TC code for brute-force checks: (x+1)(x^4+x+1), the
/// expurgated (15,10) Hamming code.
const SMALL: CyclicCode = CyclicCode { n: 15, r: 5, g: 0b11_0101 };

fn info_bytes(v: u64) -> [u8; 7] {
    v.to_be_bytes()[1..].try_into().unwrap()
}

#[test]
fn generator_factors_and_primitive_part() {
    assert_eq!(poly_mul(0b11, P6), TC_CODE.g);
    assert_eq!(poly_mul(0b11, 0b1_0011), SMALL.g);
    // x has order exactly 63 modulo x^6 + x + 1.
    let pow = |e: u32| (0..e).fold(1u64, |a, _| poly_mod(a << 1, P6, 6));
    assert_eq!(pow(63), 1);
    assert!(pow(21) != 1 && pow(9) != 1);
    // Hence x^n - 1 is divisible by g only at n = 63: the code is cyclic.
    assert_eq!(poly_mod(1u64 << 63 | 1, TC_CODE.g, 7), 0);
}

#[test]
fn three_encoders_agree() {
    let enc = TableEncoder::new();
    let mut rng = Rng::new(1);
    for t in 0..20_000 {
        let v = if t == 0 { 0 } else { rng.next_u64() >> 8 };
        let info = info_bytes(v);
        let p = enc.parity(&info);
        assert_eq!(p, parity_lfsr(&info));
        assert_eq!(TC_CODE.encode(v), v << 7 | p as u64);
        let w = block_to_word(&enc.codeblock(&info));
        assert_eq!(TC_CODE.syndrome(w), 0);
        assert_eq!(word_to_block(w), enc.codeblock(&info));
    }
    // Complemented parity: all-zero data still gives transitions.
    assert_eq!(enc.codeblock(&[0; 7]), [0, 0, 0, 0, 0, 0, 0, 0xFE]);
}

#[test]
fn weight_distribution_three_ways() {
    let a = coset_weight_distribution(&TC_CODE, 0);
    assert_eq!(a, expurgated_hamming_distribution(63));
    assert_eq!(a.iter().sum::<i128>(), 1i128 << 56);
    assert_eq!(a[1] + a[2] + a[3], 0);
    // A_4 of the Hamming code, n(n-1)(n-3)/24, kept by the even subcode.
    assert_eq!(a[4], 63 * 62 * 60 / 24);
    assert_eq!(count_weight4_from_columns(&TC_CODE), a[4]);
    assert!(a.iter().skip(1).step_by(2).all(|&v| v == 0));
    // Dual: simplex code (weight 32) plus its complements.
    let mut wts: Vec<u32> = TC_CODE.dual_words().iter().map(|u| u.count_ones()).collect();
    wts.sort();
    assert_eq!(wts[0], 0);
    assert!(wts[1..64].iter().all(|&w| w == 31) && wts[64..127].iter().all(|&w| w == 32));
    assert_eq!(wts[127], 63);
}

#[test]
fn macwilliams_and_coset_formula_match_enumeration() {
    let k = SMALL.k();
    let words: Vec<u64> = (0..1u64 << k).map(|i| SMALL.encode(i)).collect();
    let hist = |f: &dyn Fn(u64) -> u64| {
        let mut h = vec![0i128; 16];
        for &c in &words {
            h[f(c).count_ones() as usize] += 1;
        }
        h
    };
    assert_eq!(coset_weight_distribution(&SMALL, 0), hist(&|c| c));
    assert_eq!(coset_weight_distribution(&SMALL, 0), expurgated_hamming_distribution(15));
    let mut rng = Rng::new(7);
    for _ in 0..50 {
        let t = rng.next_u64() & 0x7FFF;
        assert_eq!(coset_weight_distribution(&SMALL, t), hist(&|c| c ^ t), "coset of {t:#x}");
    }
    // Shortened: leading 3 info bits forced to zero.
    let mut h = vec![0i128; 16];
    for &c in words.iter().filter(|&&c| c >> 12 == 0) {
        h[c.count_ones() as usize] += 1;
    }
    assert_eq!(shortened_distribution(&SMALL, 3), h);
}

/// Brute force over all 2^15 error patterns with a syndrome decoder for the
/// small code: exact outcome probabilities must equal the sphere formulas.
#[test]
fn outcome_probabilities_match_exhaustive_small_code() {
    let mut single = vec![None; 32];
    for i in 0..15 {
        single[SMALL.h_column(i) as usize] = Some(i);
    }
    let a = coset_weight_distribution(&SMALL, 0);
    for &p in &[0.002f64, 0.03, 0.1] {
        for mode in [DecodeMode::Detect, DecodeMode::Correct] {
            let (mut ok, mut undet, mut rej, mut last_ok) = (0.0, 0.0, 0.0, 0.0);
            for e in 0u64..1 << 15 {
                let w = e.count_ones() as i32;
                let pr = p.powi(w) * (1.0 - p).powi(15 - w);
                let s = SMALL.syndrome(e);
                let decoded = match (s, mode, single[s as usize]) {
                    (0, _, _) => Some(e),
                    (_, DecodeMode::Correct, Some(i)) => Some(e ^ 1 << i),
                    _ => None,
                };
                match decoded {
                    Some(0) => ok += pr,
                    Some(_) => undet += pr,
                    None => rej += pr,
                }
                if decoded.is_some_and(|c| c >> 12 == 0) {
                    last_ok += pr;
                }
            }
            let st = block_stats(&a, p, mode);
            assert!((st.correct - ok).abs() < 1e-12 && (st.undetected - undet).abs() < 1e-12 && (st.rejected - rej).abs() < 1e-12);
            assert!((last_block_success(&SMALL, p, 3, mode) - last_ok).abs() < 1e-12);
        }
    }
}

#[test]
fn decoder_exhaustive_up_to_three_errors() {
    let enc = TableEncoder::new();
    let dec = Decoder::new();
    let info = [0x1A, 0xCF, 0xFC, 0x1D, 0x00, 0x7E, 0x99];
    let w0 = block_to_word(&enc.codeblock(&info));
    let run = |e: u64, mode| dec.decode(&word_to_block(w0 ^ e), mode);
    for i in 0..63 {
        assert_eq!(run(1 << i, DecodeMode::Correct), BlockResult::Corrected(info, i));
        assert_eq!(run(1 << i, DecodeMode::Detect), BlockResult::Rejected);
    }
    let mut doubles = 0;
    let (mut triples, mut miscorrected) = (0, 0);
    for i in 0..63 {
        for j in i + 1..63 {
            let e2 = 1u64 << i | 1 << j;
            assert_eq!(run(e2, DecodeMode::Correct), BlockResult::Rejected);
            doubles += 1;
            for k in j + 1..63 {
                let e3 = e2 | 1 << k;
                assert_eq!(run(e3, DecodeMode::Detect), BlockResult::Rejected);
                triples += 1;
                if run(e3, DecodeMode::Correct) != BlockResult::Rejected {
                    miscorrected += 1;
                }
            }
        }
    }
    assert_eq!((doubles, triples), (1953, 39711));
    // A weight-3 error is miscorrected unless its syndrome is x^6 + x + 1, the
    // one odd syndrome that is not a single-bit column. Each weight-4
    // codeword has 4 weight-3 subsets, so the count must be 4 * A_4.
    assert_eq!(miscorrected, 4 * 9765);
    let b = coset_weight_distribution(&TC_CODE, P6);
    assert_eq!(triples - miscorrected, b[3] as i32);
}

#[test]
fn tail_sequence_is_never_a_codeblock() {
    let dec = Decoder::new();
    let t = block_to_word(&TAIL_SEQ);
    // The tail lands exactly on syndrome p(x), so it needs 3 flips to reach
    // the code and is rejected in either mode.
    assert_eq!(TC_CODE.syndrome(t), P6);
    for mode in [DecodeMode::Detect, DecodeMode::Correct] {
        assert_eq!(dec.decode(&TAIL_SEQ, mode), BlockResult::Rejected);
    }
    let b = coset_weight_distribution(&TC_CODE, t);
    assert_eq!((b[0], b[1], b[2], b[3]), (0, 0, 0, 651));
    // Every single error leaves it rejected; every double error makes it
    // correctable (p(x) + even syndrome is an odd non-p syndrome).
    let mut accepted2 = 0;
    for i in 0..63 {
        assert_eq!(dec.decode(&word_to_block(t ^ 1 << i), DecodeMode::Correct), BlockResult::Rejected);
        for j in i + 1..63 {
            if dec.decode(&word_to_block(t ^ 1 << i ^ 1 << j), DecodeMode::Correct) != BlockResult::Rejected {
                accepted2 += 1;
            }
        }
    }
    assert_eq!(accepted2, 1953);
    let p = 1e-4;
    let pr = sphere_hit_probability(&b, p, 1);
    assert!((pr / (1953.0 * p * p) - 1.0).abs() < 0.02);
}

#[test]
fn crc_and_randomizer_reference_values() {
    assert_eq!(crc16_ccitt(b"123456789"), 0x29B1);
    // First 40 bits of the TC randomizer sequence, CCSDS 231.0-B.
    let seq = randomizer_sequence(512);
    assert_eq!(&seq[..5], &[0xFF, 0x39, 0x9E, 0x5A, 0x68]);
    let bits = bytes_to_bits(&seq);
    assert_eq!(&bits[255..510], &bits[..255], "period 255");
    assert!((1..255).all(|s| bits[s..s + 255] != bits[..255]), "maximal length");
    assert_eq!(bits[..255].iter().filter(|&&b| b == 1).count(), 128);
}

#[test]
fn clean_and_singly_corrupted_cltus_round_trip() {
    let enc = TableEncoder::new();
    let dec = Decoder::new();
    let mut rng = Rng::new(3);
    for len in [0usize, 1, 2, 6, 7, 30, 249, 1017] {
        let frame = TcFrame { bypass: len % 2 == 0, scid: 0x3FF, vcid: 63, seq: len as u8, data: (0..len).map(|_| rng.next_u64() as u8).collect() };
        let bytes = frame.to_bytes();
        assert_eq!(TcFrame::from_bytes(&bytes), Ok(frame.clone()));
        let mut tx = vec![FILL; 4];
        tx.extend(build_cltu(&bytes, &enc));
        tx.extend([FILL; 4]);
        let mut bits = bytes_to_bits(&tx);
        let nb = bytes.len().div_ceil(7);
        let rx = receive_cltu(&bits, &dec, DecodeMode::Detect, 0);
        assert_eq!((rx.lock_bit, rx.blocks_accepted), (Some(32), nb));
        assert_eq!(rx.frame, Ok(frame.clone()));
        // One flip per codeblock, any position: SEC repairs all of them.
        for blk in 0..nb {
            bits[48 + 64 * blk + (rng.next_u64() % 63) as usize] ^= 1;
        }
        let rx = receive_cltu(&bits, &dec, DecodeMode::Correct, 0);
        assert_eq!((rx.corrections, rx.frame), (nb, Ok(frame.clone())));
        assert_eq!(receive_cltu(&bits, &dec, DecodeMode::Detect, 0).frame.ok(), None);
    }
}

#[test]
fn monte_carlo_agrees_with_exact() {
    let a = coset_weight_distribution(&TC_CODE, 0);
    let n = 400_000;
    for (mode, p) in [(DecodeMode::Correct, 0.02), (DecodeMode::Detect, 0.05)] {
        let ex = block_stats(&a, p, mode);
        let mc = simulate_blocks(p, mode, n, 99);
        for (e, k) in [(ex.correct, mc.correct), (ex.undetected, mc.undetected), (ex.rejected, mc.rejected)] {
            let z = (k as f64 / n as f64 - e) / (e * (1.0 - e) / n as f64).sqrt();
            assert!(z.abs() < 4.0, "{mode:?} p={p}: z = {z}");
        }
    }
}

#[test]
fn small_ber_limits() {
    // Leading terms: detect-only undetected ~ A_4 p^4; SEC reject ~ C(63,2) p^2;
    // SEC miscorrection ~ 4 A_4 p^3 (weight-3 patterns next to weight-4 words).
    let a = coset_weight_distribution(&TC_CODE, 0);
    let p = 1e-6;
    let d = block_stats(&a, p, DecodeMode::Detect);
    let s = block_stats(&a, p, DecodeMode::Correct);
    assert!((d.undetected / (9765.0 * p.powi(4)) - 1.0).abs() < 1e-3);
    assert!((s.rejected / (1953.0 * p * p) - 1.0).abs() < 1e-3);
    assert!((s.undetected / (4.0 * 9765.0 * p.powi(3)) - 1.0).abs() < 1e-3);
}
