/*
 * Astron Built this AP
 * Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-16 IST
 *
 * File front end used by the libaec interoperability test.
 *   rice_tool enc n J rsi in.u16 out.rice     (in: little-endian uint16)
 *   rice_tool dec n J rsi N in.rice out.u16   (N: sample count)
 */
#include "rice.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint8_t *slurp(const char *path, size_t *len)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t *buf = malloc(sz > 0 ? (size_t)sz : 1);
    if (buf && fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf);
        buf = NULL;
    }
    fclose(f);
    *len = (size_t)sz;
    return buf;
}

static int spit(const char *path, const void *buf, size_t len)
{
    FILE *f = fopen(path, "wb");
    if (!f)
        return -1;
    size_t w = fwrite(buf, 1, len, f);
    fclose(f);
    return w == len ? 0 : -1;
}

int main(int argc, char **argv)
{
    if (argc < 7) {
        fprintf(stderr, "usage: %s enc n J rsi in out | dec n J rsi N in out\n", argv[0]);
        return 2;
    }
    rice_params p = {atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), -1};
    size_t len;

    if (!strcmp(argv[1], "enc")) {
        uint8_t *raw = slurp(argv[5], &len);
        if (!raw)
            return 1;
        size_t N = len / 2;
        uint32_t *x = malloc(N * sizeof *x);
        for (size_t i = 0; i < N; i++)
            x[i] = (uint32_t)raw[2 * i] | ((uint32_t)raw[2 * i + 1] << 8);
        uint8_t *out;
        size_t nb = rice_encode(x, N, &p, &out, NULL);
        if (!nb || spit(argv[6], out, nb))
            return 1;
        free(out);
        free(x);
        free(raw);
        return 0;
    }
    if (!strcmp(argv[1], "dec") && argc >= 8) {
        size_t N = strtoul(argv[5], NULL, 10);
        uint8_t *in = slurp(argv[6], &len);
        if (!in)
            return 1;
        uint32_t *x = malloc(N * sizeof *x);
        uint8_t *raw = malloc(2 * N);
        if (rice_decode(in, len, &p, x, N))
            return 1;
        for (size_t i = 0; i < N; i++) {
            raw[2 * i] = (uint8_t)x[i];
            raw[2 * i + 1] = (uint8_t)(x[i] >> 8);
        }
        if (spit(argv[7], raw, 2 * N))
            return 1;
        free(raw);
        free(x);
        free(in);
        return 0;
    }
    return 2;
}
