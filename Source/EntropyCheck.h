/*
    This file is part of EntropyCheck.

    Copyright (C) 2021 ReimuNotMoe <reimu@sudomaker.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the Apache License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    This file incorporates public domain works.
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

/* Bytes used as Monte Carlo co-ordinates.
 * This should be no more bits than the mantissa
 * of your "double" floating point type. */
#define __ENTCHK_MONTEN			6

typedef struct {
	int binary;			// Treat input as a bitstream
	uint64_t ccount[256];		// Bins to count occurrences of values
	uint64_t totalc;		// Total bytes counted
	double prob[256];		// Probabilities per bin for entropy

	int mp, sccfirst;
	unsigned int monte[__ENTCHK_MONTEN];
	long inmont, mcount;
	double cexp, incirc, montex, montey, montepi;
	double scc, sccun, sccu0, scclast, scct1, scct2, scct3;
	double ent, chisq, datasum;
} EntropyCheckContext;

typedef struct {
	uint64_t totalc;	// Total bytes counted
	double ent;		// Entropy bits per sample
	double chisq;		// Chi Square distribution
	double pochisq;		// Probability of measured Chi Square value
	double mean;		// Arithmetic mean value
	double montepi;		// Monte Carlo value for Pi
	double scc;		// Serial correlation coefficient
} EntropyCheckResult;


extern void EntropyCheck_Init(EntropyCheckContext *ctx, bool binary_mode);
extern void EntropyCheck_Feed(EntropyCheckContext *ctx, void *buf, size_t buf_len);
extern void EntropyCheck_Finalize(EntropyCheckContext *ctx, EntropyCheckResult *res);

#ifdef __cplusplus
}
#endif
