/*
    This file is part of EntropyCheck.

    Copyright (C) 2021 ReimuNotMoe <reimu@sudomaker.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the Apache License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include <EntropyCheck.hpp>

#include <cmath>
#include <cassert>
#include <cstdio>
#include <cinttypes>

#include <unistd.h>
#include <fcntl.h>


int main() {
	EntropyChecker ec;

	int fd = open("/dev/urandom", O_RDONLY);
	assert(fd > 0);

	uint8_t buf[32];

	for (size_t i=0; i<1024; i++) {
		ssize_t rc = read(fd, buf, sizeof(buf));
		assert(rc == sizeof(buf));

		ec.feed(buf, sizeof(buf));
	}

	auto results = ec.finalize();

	{
		printf("Entropy = %f bits per byte.\n", results.ent);

		printf("\nOptimum compression would reduce the size\n"
		       "of this %" PRIu64 " bytes file by %d percent.\n\n",
		       results.totalc, (short) ((100 * (8 - results.ent) / 8.0)));

		printf("Chi square distribution for %" PRIu64 " samples is %1.2f, and randomly\n",
		       results.totalc, results.chisq);
		if (results.pochisq < 0.0001) {
			printf("would exceed this value less than 0.01 percent of the times.\n\n");
		} else if (results.pochisq > 0.9999) {
			printf("would exceed this value more than than 99.99 percent of the times.\n\n");
		} else {
			printf("would exceed this value %1.2f percent of the times.\n\n",
			       results.pochisq * 100);
		}
		printf("Arithmetic mean value of data is %1.4f (%.1f = random).\n",
			results.mean, 127.5);

		printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
		       results.montepi, 100.0 * (fabs(M_PI - results.montepi) / M_PI));
		printf("Serial correlation coefficient is ");
		if (results.scc >= -99999) {
			printf("%1.6f (totally uncorrelated = 0.0).\n", results.scc);
		} else {
			printf("undefined (all values equal!).\n");
		}
	}

	return 0;
}