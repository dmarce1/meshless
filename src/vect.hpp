/*
 * general_vect.hpp
 *
 *  Created on: Nov 30, 2019
 *      Author: dmarce1
 */

#ifndef VECT_HPP_
#define VECT_HPP_

#include "dim.hpp"
#include "real.hpp"

#include <array>
#include <atomic>
#include <cmath>

template<class T, int N>
class general_vect {
	std::array<T, N> v;
public:
	general_vect() = default;
	general_vect(std::array<real, N> a) :
			v(a) {
	}
	T& operator[](int i);
	T operator[](int i) const;
	general_vect operator-() const;
	general_vect operator-(const general_vect &other) const;
	general_vect operator+(const general_vect &other) const;
	general_vect operator*(T r) const;
	general_vect operator/(T r) const;
	T dot(const general_vect &other) const;
	bool operator!=(const general_vect &other) const {
		for (int dim = 0; dim < NDIM; dim++) {
			if (other[dim] != (*this)[dim]) {
				return true;
			}
		}
		return false;
	}
};

template<class T, int N>
inline T& general_vect<T, N>::operator[](int i) {
	return v[i];
}

template<class T, int N>
inline T general_vect<T, N>::operator[](int i) const {
	return v[i];
}

template<class T, int N>
inline general_vect<T, N> general_vect<T, N>::operator-() const {
	general_vect<T, N> result;
	for (int dim = 0; dim < N; dim++) {
		result[dim] = -v[dim];
	}
	return result;
}

template<class T, int N>
inline general_vect<T, N> general_vect<T, N>::operator-(const general_vect<T, N> &other) const {
	general_vect<T, N> result;
	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] - other[dim];
	}
	return result;
}

template<class T, int N>
inline general_vect<T, N> general_vect<T, N>::operator+(const general_vect<T, N> &other) const {
	general_vect<T, N> result;
	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] + other[dim];
	}
	return result;
}

template<class T, int N>
inline general_vect<T, N> general_vect<T, N>::operator*(T r) const {
	general_vect<T, N> result;
	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] * r;
	}
	return result;
}

template<class T, int N>
inline general_vect<T, N> general_vect<T, N>::operator/(T r) const {
	general_vect<T, N> result;
	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] / r;
	}
	return result;
}

template<class T, int N>
inline T general_vect<T, N>::dot(const general_vect<T, N> &other) const {
	T result = 0.0;
	for (int dim = 0; dim < N; dim++) {
		result += v[dim] * other[dim];
	}
	return result;
}

template<class T, int N>
inline T abs(const general_vect<T, N> &v) {
	return std::sqrt(v.dot(v));
}

template<class T, int N>
inline general_vect<T, N> max(const general_vect<T, N> &a, const general_vect<T, N> &b) {
	general_vect<T, N> c;
	for (int i = 0; i < N; i++) {
		c[i] = std::max(a[i], b[i]);
	}
	return c;
}

template<class T, int N>
inline general_vect<T, N> min(const general_vect<T, N> &a, const general_vect<T, N> &b) {
	general_vect<T, N> c;
	for (int i = 0; i < N; i++) {
		c[i] = std::min(a[i], b[i]);
	}
	return c;
}

using vect = general_vect<real, NDIM>;
using atomic_vect = general_vect<std::atomic<real>, NDIM>;

#endif /* VECT_HPP_ */
