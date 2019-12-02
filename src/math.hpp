/*
 * math.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include "vect.hpp"

using range = std::pair<vect,vect>;

real rand_unit_box();
real W(real, real);
real dW_dh(real, real);
real distance(const vect&, const vect&);
bool in_range(const vect&, const range&);
bool ranges_intersect(const range&, const range&);
real bspline(real r);
real d_bspline_dr(real r);
real W_norm(real h);
real dW_norm_dh(real h);
inline real box_volume(const range& R);

inline bool in_range(const vect &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.first[dim] || x[dim] > r.second[dim]) {
			return false;
		}
	}
	return true;
}

inline bool ranges_intersect(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (a.second[dim] < b.first[dim] || b.second[dim] < a.first[dim]) {
			return false;
		}
	}
	return true;
}

inline real W(real r, real h) {
	auto norm = W_norm(h);
	return bspline(2 * r / h) / norm;
}

inline real dW_dh(real r, real h) {
	auto norm = W_norm(h);
	return -2 * r * d_bspline_dr(2 * r / h) / (h * h) / norm - dW_norm_dh(h) * W(r, h) / norm;
}

inline real bspline(real r) {
	if (r < 1) {
		const auto r2 = r * r;
		return 1.0 - 1.5 * r2 + 0.75 * r2 * r;
	} else if (r <= 2) {
		return 0.25 * std::pow(2 - r, 3);
	} else {
		return 0.0;
	}
}

inline real d_bspline_dr(real r) {
	if (r < 1) {
		return -3.0 * r + 2.25 * r * r;
	} else if (r <= 2) {
		return -0.75 * std::pow(2 - r, 2);
	} else {
		return 0.0;
	}
}

inline real W_norm(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 * h / 4.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 7 * h * h / 80.0;
	} else {
		return 2 * 2.0 * M_PI * h * h * h / 32.0;
	}
}

inline real dW_norm_dh(real h) {
	if constexpr (NDIM == 1) {
		return 3.0 / 4.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 15 * h / 80.0;
	} else {
		return 2 * 6.0 * M_PI * h * h / 32.0;
	}
}

inline real box_volume(const range& R) {
	real vol = 1.0;
	for( int dim = 0; dim <NDIM; dim++) {
		vol *= R.second[dim] - R.first[dim];
	}
	return vol;
}

inline real rand_unit_box() {
	const auto x = (rand() + 0.5) / (real(RAND_MAX) + 1.0);
	return x - 0.5;
}

inline std::array<vect,NDIM> matrix_inverse(const std::array<vect,NDIM>& A) {
	std::array<vect,NDIM> Ainv;
#if(NDIM==1)
	Ainv[0][0] = 1.0 / A[0][0];
#else
#if(NDIM==2)
	real det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
	const real detinv = 1.0 / detinv;
	Ainv[0][0] = +A[1][1] * detinv;
	Ainv[0][1] = -A[1][0] * detinv;
	Ainv[1][0] = -A[0][1] * detinv;
	Ainv[1][1] = +A[0][0] * detinv;
#else
	real det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
	det /**/-= A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
	det /**/+= A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	const real detinv = 1.0 / det;
	Ainv[0][0] = (-A[1][2] * A[2][1] + A[1][1] * A[2][2]) * detinv;
	Ainv[0][1] = (+A[0][2] * A[2][1] - A[0][1] * A[2][2]) * detinv;
	Ainv[0][2] = (-A[0][2] * A[1][1] + A[0][1] * A[1][2]) * detinv;
	Ainv[1][0] = (+A[1][2] * A[2][0] - A[1][0] * A[2][2]) * detinv;
	Ainv[1][1] = (-A[0][2] * A[2][0] + A[0][0] * A[2][2]) * detinv;
	Ainv[1][2] = (+A[0][2] * A[1][0] - A[0][0] * A[1][2]) * detinv;
	Ainv[2][0] = (-A[1][1] * A[2][0] + A[1][0] * A[2][1]) * detinv;
	Ainv[2][1] = (+A[0][1] * A[2][0] - A[0][0] * A[2][1]) * detinv;
	Ainv[2][2] = (-A[0][1] * A[1][0] + A[0][0] * A[1][1]) * detinv;
#endif
#endif
	return Ainv;
}


inline real distance(const vect &a, const vect &b) {
	real d = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		d += std::pow(a[dim] - b[dim], 2);
	}
	return std::sqrt(d);
}

#endif /* MATH_HPP_ */

