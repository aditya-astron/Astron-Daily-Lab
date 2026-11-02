/*
 * Astron Built this AP
 * Project 008 | Error correction | CCSDS K=7 rate-1/2 convolutional code with soft-decision Viterbi decoding
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-02 IST
 */
/* Monte Carlo BER of the CCSDS K=7 R=1/2 code over BPSK/AWGN with hard,
 * 3-bit and unquantized soft decisions, against uncoded BPSK and the
 * soft-decision union bound  Pb <= sum_d B_d Q(sqrt(2 d R Eb/N0)). */
#include "channel.h"
#include "viterbi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK   4096
#define DMAX    24

static double union_bound(const double *b, double ebn0_db)
{
    double ebn0 = pow(10.0, ebn0_db / 10.0), pb = 0.0;
    for (int d = 0; d <= DMAX; d++)
        pb += b[d] * q_func(sqrt(2.0 * d * 0.5 * ebn0));
    return pb;
}

int main(int argc, char **argv)
{
    long max_bits = argc > 1 ? atol(argv[1]) : 8000000L;
    /* Viterbi errors arrive in bursts, so count many of them and never stop
     * on fewer than ~50 blocks. */
    const long min_errors = 500, min_bits = 50L * BLOCK;

    double a[DMAX + 1], b[DMAX + 1];
    int dfree = cc_distance_spectrum(DMAX, a, b);
    printf("CCSDS convolutional code K=7 R=1/2 (171, 133 octal)\n");
    printf("d_free = %d, asymptotic soft coding gain = 10log10(R*dfree) = %.2f dB\n\n",
           dfree, 10.0 * log10(0.5 * dfree));
    printf(" d    a_d        B_d\n");
    for (int d = dfree; d <= 20; d += 2)
        printf("%2d %6.0f %10.0f\n", d, a[d], b[d]);

    static uint8_t info[BLOCK], code[2 * (BLOCK + CC_MEM)], dec[BLOCK];
    static double rx[2 * (BLOCK + CC_MEM)];
    const enum rx_mode modes[3] = { RX_HARD, RX_Q3, RX_SOFT };

    printf("\nEb/N0  uncoded    hard       3-bit      soft       union bound\n");
    for (double snr = 1.0; snr <= 5.51; snr += 0.5) {
        double ber[3];
        for (int k = 0; k < 3; k++) {
            rng_t r;
            rng_seed(&r, 2026u + (uint64_t)(snr * 10) * 7 + (uint64_t)k);
            long errs = 0, bits = 0;
            while (bits < max_bits && (errs < min_errors || bits < min_bits)) {
                for (int i = 0; i < BLOCK; i++)
                    info[i] = (uint8_t)(rng_next(&r) >> 63);
                cc_encode(info, BLOCK, code);
                bpsk_awgn(code, 2 * (BLOCK + CC_MEM), snr, 0.5, modes[k], &r, rx);
                if (cc_viterbi_decode(rx, BLOCK, dec) != 0) {
                    fprintf(stderr, "out of memory\n");
                    return 1;
                }
                for (int i = 0; i < BLOCK; i++)
                    errs += info[i] != dec[i];
                bits += BLOCK;
            }
            ber[k] = errs ? (double)errs / bits : 0.0;
        }
        double ub = union_bound(b, snr);
        printf("%4.1f   %.3e  %.3e  %.3e  %.3e  %.3e\n", snr,
               q_func(sqrt(2.0 * pow(10.0, snr / 10.0))),
               ber[0], ber[1], ber[2], ub < 1.0 ? ub : 1.0);
    }
    printf("\n(0 means no errors in %ld bits)\n", max_bits);
    return 0;
}
