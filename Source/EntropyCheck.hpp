/*
    This file is part of EntropyCheck.

    Copyright (C) 2021 ReimuNotMoe <reimu@sudomaker.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the Apache License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#pragma once

#include "EntropyCheck.h"

class EntropyChecker {
private:
	EntropyCheckContext ctx;
public:
	EntropyChecker(bool binary_mode = false) {
		EntropyCheck_Init(&ctx, binary_mode);
	}

	void reset() {
		EntropyCheck_Init(&ctx, ctx.binary);
	}

	void feed(void *buf, size_t len) {
		EntropyCheck_Feed(&ctx, buf, len);
	}

	void finalize(EntropyCheckResult *res) {
		EntropyCheck_Finalize(&ctx, res);
	}

	EntropyCheckResult finalize() {
		EntropyCheckResult res;
		EntropyCheck_Finalize(&ctx, &res);
		return res;
	}
};