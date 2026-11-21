// Astron Built this AP
// Project 051 | Telemetry framing / Error correction | CCSDS telecommand BCH(63,56) codec and CLTU link
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST

use tc_cltu::*;

const N_BLOCKS_FRAME_DATA: usize = 64;

fn sigma(p: f64, n: u64) -> f64 {
    (p * (1.0 - p) / n as f64).sqrt()
}

fn frame_for_test(rng: &mut Rng, seq: u8) -> TcFrame {
    let data = (0..N_BLOCKS_FRAME_DATA).map(|_| rng.next_u64() as u8).collect();
    TcFrame { bypass: true, scid: 0x2A5, vcid: 3, seq, data }
}

fn main() {
    let code = TC_CODE;
    let a = coset_weight_distribution(&code, 0);
    println!("CCSDS TC code: g(x) = x^7 + x^6 + x^2 + 1 = (x + 1)(x^6 + x + 1), (n, k) = ({}, {})", code.n, code.k());
    let dmin = a.iter().enumerate().skip(1).find(|(_, &v)| v != 0).unwrap().0;
    println!("d_min = {dmin}, codewords = 2^{}", code.k());
    println!("weight distribution (MacWilliams from the 128-word dual, = closed form):");
    for w in [4, 6, 8, 30, 32] {
        println!("   A_{w:<2} = {}", a[w]);
    }
    let dual_wts: Vec<u32> = code.dual_words().iter().map(|u| u.count_ones()).collect();
    let cnt = |w| dual_wts.iter().filter(|&&x| x == w).count();
    println!("   dual weights: 0 x{}, 31 x{}, 32 x{}, 63 x{}", cnt(0), cnt(31), cnt(32), cnt(63));

    // 1) Codeblock outcome probabilities, exact.
    println!("\n1) Codeblock outcomes on a BSC (exact, from the weight distribution)");
    println!("    BER   |  detect-only: reject      undetected  |  SEC: correct+fixed  reject      miscorrect");
    for &p in &[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2] {
        let d = block_stats(&a, p, DecodeMode::Detect);
        let s = block_stats(&a, p, DecodeMode::Correct);
        println!(
            "   {p:7.0e} |   {:10.3e}      {:10.3e}   |   {:.9}      {:10.3e}  {:10.3e}",
            d.rejected, d.undetected, s.correct, s.rejected, s.undetected
        );
    }

    // 2) Monte Carlo against the exact numbers.
    println!("\n2) Monte Carlo decoder vs exact (2e6 random codeblocks each)");
    println!("   mode     BER    | outcome      exact        simulated    z");
    for &(mode, p) in &[(DecodeMode::Detect, 0.03), (DecodeMode::Correct, 0.01), (DecodeMode::Correct, 0.03)] {
        let n = 2_000_000;
        let ex = block_stats(&a, p, mode);
        let mc = simulate_blocks(p, mode, n, 0x5151 + (p * 1e4) as u64);
        for (name, e, k) in [("correct", ex.correct, mc.correct), ("reject", ex.rejected, mc.rejected), ("undetected", ex.undetected, mc.undetected)] {
            let f = k as f64 / n as f64;
            println!("   {:<8} {p:<6} | {name:<10}  {e:.6e}  {f:.6e}  {:+5.2}", format!("{mode:?}"), (f - e) / sigma(e, n));
        }
    }

    // 3) The tail sequence as a received codeblock.
    let tail = block_to_word(&TAIL_SEQ);
    let b = coset_weight_distribution(&code, tail);
    let leader = b.iter().position(|&v| v != 0).unwrap();
    println!("\n3) Tail sequence C5C5C5C5C5C5C579");
    println!("   syndrome = 0x{:02X}, filler bit = {}, coset leader weight = {leader}", code.syndrome(tail), TAIL_SEQ[7] & 1);
    println!("   coset weights: B_{leader} = {}, B_{} = {}", b[leader], leader + 1, b[leader + 1]);
    println!("    BER   |  P(tail accepted, SEC)   P(tail accepted, detect)");
    for &p in &[1e-5, 1e-3, 1e-2] {
        println!("   {p:7.0e} |     {:10.3e}             {:10.3e}", sphere_hit_probability(&b, p, 1), sphere_hit_probability(&b, p, 0));
    }

    // 4) End-to-end CLTU link.
    let enc = TableEncoder::new();
    let dec = Decoder::new();
    let mut rng = Rng::new(20_261_121);
    let len = frame_for_test(&mut rng, 0).to_bytes().len();
    let nb = len.div_ceil(7);
    println!("\n4) CLTU link: {len}-octet TC frame -> {nb} codeblocks, 8 octets of 0x55 idle either side");
    let pre = bytes_to_bits(&[FILL, FILL, 0xEB, 0x90]);
    let side = (0..16)
        .map(|k| pre[k..k + 16].iter().fold(0u16, |a, &b| a << 1 | b as u16))
        .map(|w| (w ^ START_SEQ).count_ones())
        .min()
        .unwrap();
    let last_bits = 8 * (len - 7 * (nb - 1));
    println!("   start sequence vs every idle/start straddling window: min distance {side}");
    println!("   last block carries {} frame bits + {} fill bits", last_bits, 56 - last_bits);
    println!("   config            BER   | delivered (exact)  delivered (sim)    z   | CRC-caught  undetected  false locks");
    let trials = 200_000u64;
    for &(mode, tol, p) in &[
        (DecodeMode::Correct, 0, 1e-3),
        (DecodeMode::Correct, 0, 5e-3),
        (DecodeMode::Correct, 1, 5e-3),
        (DecodeMode::Detect, 0, 5e-3),
        (DecodeMode::Correct, 1, 2e-2),
    ] {
        let st = block_stats(&a, p, mode);
        let q = 1.0 - p;
        let p_start = q.powi(16) + if tol == 1 { 16.0 * p * q.powi(15) } else { 0.0 };
        // Every block must decode to the sent data, except that only the frame
        // octets of the last block matter; the rest of it is fill.
        let last = last_block_success(&code, p, 8 * (len - 7 * (nb - 1)) as u32, mode);
        let exact = p_start * st.correct.powi(nb as i32 - 1) * last;
        let (mut ok, mut crc, mut bad, mut false_lock) = (0u64, 0u64, 0u64, 0u64);
        let idle = [FILL; 8];
        for t in 0..trials {
            let frame = frame_for_test(&mut rng, t as u8);
            let mut tx = idle.to_vec();
            tx.extend(build_cltu(&frame.to_bytes(), &enc));
            tx.extend_from_slice(&idle);
            let mut bits = bytes_to_bits(&tx);
            for bit in bits.iter_mut() {
                if rng.uniform() < p {
                    *bit ^= 1;
                }
            }
            let rx = receive_cltu(&bits, &dec, mode, tol);
            if rx.lock_bit.is_some_and(|l| l != 64) {
                false_lock += 1;
            }
            match rx.frame {
                Ok(f) if f == frame => ok += 1,
                Ok(_) => bad += 1,
                Err(FrameError::Crc) => crc += 1,
                Err(_) => {}
            }
        }
        let f = ok as f64 / trials as f64;
        println!(
            "   {:<7} start<={tol} {p:<6.0e}|    {exact:.5}          {f:.5}       {:+5.2}  | {crc:>8}    {bad:>8}    {false_lock:>8}",
            format!("{mode:?}"),
            (f - exact) / sigma(exact, trials)
        );
    }
    println!("   ({trials} CLTUs per row)");
}
