# Astron Built this AP
# Project 025 | Signal processing | Polyphase DFT channelizer for a multi-satellite ground receiver
# Author: Aditya (aditya-astron)
# Date: 2026-11-09 IST
"""
Critically sampled polyphase DFT analysis filter bank.

A ground-station SDR digitises an fs = 8 Msps complex band that holds up to
M = 16 CubeSat downlinks on a 500 kHz raster. The channelizer splits it into
M baseband streams at fs/M each, at the cost of one prototype FIR plus one
M-point FFT per M input samples.

Channel k is centred at f_k = k fs/M (k > M/2 are the negative frequencies).
The reference definition is mix-down, low-pass, keep every M-th sample:

    y_k[m] = sum_l h[l] x[mM - l] exp(-j 2 pi k (mM - l)/M)
           = sum_l h[l] x[mM - l] exp(+j 2 pi k l / M)        (e^{-j2pi k m} = 1)

Split l = rM + p (p = 0..M-1):

    y_k[m] = sum_p exp(+j 2 pi k p / M) * sum_r h_p[r] x_p[m - r]
    h_p[r] = h[rM + p]         (polyphase branch p of the prototype)
    x_p[m] = x[mM - p]         (input commutator, branch p)

so each output instant is M short FIRs running at the low rate followed by an
inverse DFT across the branches:  y[:, m] = M * IFFT_p(u_p[m]).
"""
import numpy as np
from scipy import signal


def design_prototype(fs, M, pass_frac=0.4, stop_frac=0.6, atten_db=80.0):
    """Kaiser-window low-pass for a channel spacing of fs/M.

    Passband edge pass_frac*fs/M, stopband edge stop_frac*fs/M, cutoff at the
    midpoint (half the spacing). Content between pass_frac and 0.5 of the
    spacing would alias into the outer edge of the band after decimation, so
    signals must occupy at most +/- pass_frac of the spacing; anything beyond
    stop_frac is attenuated by atten_db before it can fold back.

    The tap count is forced odd so the group delay (N-1)/2 is an integer, and
    a multiple of M in input samples so it is also an integer number of
    output samples. The filter is then zero-padded to a multiple of M taps.
    """
    spacing = fs / M
    width = (stop_frac - pass_frac) * spacing / (fs / 2)      # Nyquist-normalised
    n_est, beta = signal.kaiserord(atten_db, width)
    # smallest odd N >= n_est with (N-1)/2 divisible by M
    half = int(np.ceil((n_est - 1) / 2 / M)) * M
    n_taps = 2 * half + 1
    cutoff = 0.5 * spacing
    h = signal.firwin(n_taps, cutoff, window=("kaiser", beta), fs=fs)
    pad = (-n_taps) % M
    h = np.concatenate([h, np.zeros(pad)])
    return h, half // M                       # taps, group delay in output samples


def channelize(x, h, M):
    """Polyphase DFT analysis bank. Returns array (M, len(x)//M)."""
    if len(h) % M:
        raise ValueError("prototype length must be a multiple of M")
    n_out = len(x) // M
    x = np.asarray(x, dtype=complex)[: n_out * M]
    # x_p[m] = x[mM - p]: prepend M-1 zeros so negative indices read as zero
    xpad = np.concatenate([np.zeros(M - 1, dtype=complex), x])
    idx = (np.arange(n_out)[None, :] * M) - np.arange(M)[:, None] + (M - 1)
    xp = xpad[idx]                                   # (M, n_out)
    hp = h.reshape(-1, M).T                          # (M, taps per branch)
    u = np.empty_like(xp)
    for p in range(M):
        u[p] = signal.lfilter(hp[p], 1.0, xp[p])
    return M * np.fft.ifft(u, axis=0)


def channelize_direct(x, h, M, channels=None):
    """Reference: per-channel mix, full-rate FIR, decimate by M."""
    x = np.asarray(x, dtype=complex)
    n_out = len(x) // M
    n = np.arange(len(x))
    ks = range(M) if channels is None else channels
    out = np.empty((len(ks), n_out), dtype=complex)
    for i, k in enumerate(ks):
        bb = x * np.exp(-2j * np.pi * k * n / M)
        out[i] = signal.lfilter(h, 1.0, bb)[: n_out * M : M]
    return out


def channel_freq(k, fs, M):
    """Centre frequency of channel k in the complex band [-fs/2, fs/2)."""
    f = k * fs / M
    return f - fs if f >= fs / 2 else f


def mults_per_input_sample(n_taps, M):
    """Real multiplies per complex input sample (complex x real tap = 2,
    complex x complex = 4; radix-2 FFT ~ (M/2) log2 M complex multiplies)."""
    direct = M * (4 + 2 * n_taps)                    # every channel: mixer + FIR
    poly = (2 * n_taps + 4 * (M / 2) * np.log2(M)) / M
    return direct, poly


def rrc(beta, sps, span):
    """Root-raised-cosine taps, unit energy, span symbols long."""
    t = np.arange(-span * sps / 2, span * sps / 2 + 1) / sps
    h = np.empty_like(t)
    for i, ti in enumerate(t):
        if abs(ti) < 1e-12:
            h[i] = 1 - beta + 4 * beta / np.pi
        elif abs(abs(ti) - 1 / (4 * beta)) < 1e-12:
            h[i] = beta / np.sqrt(2) * ((1 + 2 / np.pi) * np.sin(np.pi / (4 * beta))
                                        + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta)))
        else:
            h[i] = (np.sin(np.pi * ti * (1 - beta)) + 4 * beta * ti * np.cos(np.pi * ti * (1 + beta))) \
                   / (np.pi * ti * (1 - (4 * beta * ti) ** 2))
    return h / np.sqrt(np.sum(h ** 2))


def qpsk_burst(n_sym, sps, beta, span, rng):
    sym = ((2 * rng.integers(0, 2, n_sym) - 1) + 1j * (2 * rng.integers(0, 2, n_sym) - 1)) / np.sqrt(2)
    up = np.zeros(n_sym * sps, dtype=complex)
    up[::sps] = sym
    return sym, signal.lfilter(rrc(beta, sps, span), 1.0, up)


def evm_db(rx, ref):
    """EVM after a least-squares complex gain fit."""
    g = np.vdot(rx, ref) / np.vdot(rx, rx)
    err = g * rx - ref
    return 10 * np.log10(np.mean(np.abs(err) ** 2) / np.mean(np.abs(ref) ** 2))


def scenario(seed=7):
    """Three QPSK downlinks plus two strong neighbours, through the bank."""
    fs, M = 8e6, 16
    h, d_out = design_prototype(fs, M)
    # The transmit RRC is long so that its own truncation sidelobes (-81 dB
    # beyond 300 kHz at 48 symbols, -53 dB at 12) do not mask the channelizer.
    rsym, beta, tx_span, rx_span = 250e3, 0.35, 48, 16
    sps_in, sps_out = int(fs / rsym), int(fs / M / rsym)          # 32, 2
    rng = np.random.default_rng(seed)
    n_sym = 3000
    n = np.arange(n_sym * sps_in)

    # (channel, relative level dB). Channels 4 and 6 are 30 dB hotter than
    # their neighbour 5, the weak link of interest.
    plan = [(2, 0.0), (5, 0.0), (11, -6.0), (4, 30.0), (6, 30.0)]
    x = np.zeros(len(n), dtype=complex)
    syms = {}
    for k, lvl in plan:
        s, b = qpsk_burst(n_sym, sps_in, beta, tx_span, rng)
        syms[k] = s
        x += 10 ** (lvl / 20) * b * np.exp(2j * np.pi * k * n / M)

    y = channelize(x, h, M)
    g_rx = rrc(beta, sps_out, rx_span)
    # tx RRC delay + channelizer + rx RRC delay, in output samples
    delay = tx_span * sps_in // 2 // M + d_out + rx_span * sps_out // 2
    results = []
    for k, lvl in plan:
        z = signal.lfilter(g_rx, 1.0, y[k])
        rx = z[delay::sps_out][: n_sym]
        keep = slice(2 * tx_span, len(rx) - 2 * tx_span)        # drop filter edges
        results.append((k, lvl, channel_freq(k, fs, M), evm_db(rx[keep], syms[k][keep])))
    empty = [k for k in range(M) if k not in dict(plan)]
    floor = max(10 * np.log10(np.mean(np.abs(y[k, 200:]) ** 2)) for k in empty)
    ref = 10 * np.log10(np.mean(np.abs(y[4, 200:]) ** 2))
    return dict(fs=fs, M=M, h=h, d_out=d_out, results=results, floor_dbc=floor - ref)


def main():
    fs, M = 8e6, 16
    h, d_out = design_prototype(fs, M)
    w, H = signal.freqz(h, worN=1 << 16, fs=fs)
    mag = 20 * np.log10(np.abs(H) + 1e-300)
    sp = fs / M
    ripple = np.ptp(mag[w <= 0.4 * sp])
    stop = -mag[w >= 0.6 * sp].max()
    direct, poly = mults_per_input_sample(len(h), M)
    print(f"fs = {fs/1e6:.0f} Msps, M = {M}, spacing {sp/1e3:.0f} kHz, output {fs/M/1e3:.0f} ksps per channel")
    print(f"prototype: {len(h)} taps ({len(h)//M} per branch), group delay {d_out} output samples")
    print(f"  passband ripple {ripple:.5f} dB (0-200 kHz), stopband {stop:.1f} dB (>= 300 kHz)")
    print(f"real multiplies per input sample: direct {direct:.0f}, polyphase {poly:.1f} "
          f"({direct/poly:.0f}x fewer)")

    sc = scenario()
    print("\nQPSK 250 ksym/s, RRC 0.35, 3000 symbols per link, no noise")
    print(" ch   centre kHz   level dB   EVM dB")
    for k, lvl, f, e in sc["results"]:
        print(f"{k:3d}   {f/1e3:+9.0f}   {lvl:+7.1f}   {e:7.1f}")
    print(f"worst empty-channel power: {sc['floor_dbc']:.1f} dB relative to channel 4")


if __name__ == "__main__":
    main()
