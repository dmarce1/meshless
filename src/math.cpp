/*
 * math.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#include "math.hpp"

#include <cmath>

static real bspline(real r);
static real d_bspline_dr(real r);
static real W_norm(real h);
static real dW_norm_dh(real h);

real rand_unit_box() {
	const auto x = (rand() + 0.5) / (real(RAND_MAX) + 1.0);
	return 2.0 * (x - 0.5);
}

real W(real r, real h) {
	auto norm = W_norm(h);
	return bspline(2 * r / h) / norm;
}

real dW_dh(real r, real h) {
	auto norm = W_norm(h);
	return -2 * r * d_bspline_dr(2 * r / h) / (h * h) / norm - dW_norm_dh(h) * W(r, h) / norm;
}

real distance(const vect &a, const vect &b) {
	real d = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		d += std::pow(a[dim] - b[dim], 2);
	}
	return std::sqrt(d);
}

real max_span(const range &r) {
	real m = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		m = std::max(m, r.second[dim] - r.first[dim]);
	}
	return m;
}

bool in_range(const vect &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.first[dim] || x[dim] > r.second[dim]) {
			return false;
		}
	}
	return true;
}

bool ranges_intersect(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (a.second[dim] < b.first[dim] || b.second[dim] < a.first[dim]) {
			return false;
		}
	}
	return true;
}

static real bspline(real r) {
	if (r < 1) {
		const auto r2 = r * r;
		return 1.0 - 1.5 * r2 + 0.75 * r2 * r;
	} else if (r <= 2) {
		return 0.25 * std::pow(2 - r, 3);
	} else {
		return 0.0;
	}
}

static real d_bspline_dr(real r) {
	if (r < 1) {
		return -3.0 * r + 2.25 * r * r;
	} else if (r <= 2) {
		return -0.75 * std::pow(2 - r, 2);
	} else {
		return 0.0;
	}
}

static real W_norm(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 * h / 8.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 7 * h * h / 80.0;
	} else {
		return 2 * 2.0 * M_PI * h * h * h / 32.0;
	}
}

static real dW_norm_dh(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 / 8.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 15 * h / 80.0;
	} else {
		return 2 * 6.0 * M_PI * h * h / 32.0;
	}
}
