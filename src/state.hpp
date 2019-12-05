/*
 * state.hpp
 *
 *  Created on: Dec 2, 2019
 *      Author: dmarce1
 */

#ifndef STATE_HPP_
#define STATE_HPP_

#include <cassert>

static int constexpr STATE_SIZE = 2 + NDIM;
static double constexpr fgamma = 7.0 / 5.0;

struct state: public general_vect<real, STATE_SIZE> {
	state(const general_vect<real, STATE_SIZE> &other) :
			general_vect<real, STATE_SIZE>(other) {
	}
	state() = default;
	state& operator=(const general_vect<real, STATE_SIZE> &other) {
		for (int f = 0; f < STATE_SIZE; f++) {
			(*this)[f] = other[f];
		}
		return *this;
	}
	real mass() const {
		assert( (*this)[0] > 0.0);
		return (*this)[0];
	}
	real& mass() {
		return (*this)[0];
	}
	real energy() const {
		return (*this)[1];
	}
	real& energy() {
		return (*this)[1];
	}
	vect momentum() const {
		vect m;
		for (int dim = 0; dim < NDIM; dim++) {
			m[dim] = (*this)[2 + dim];
		}
		return m;
	}
	real momentum(int i) const {
		return (*this)[2 + i];
	}
	real& momentum(int i) {
		return (*this)[2 + i];
	}
	std::pair<state, real> flux(const vect &vf, const vect &norm) const {
		std::pair<state, real> p;
		state &f = p.first;
		const auto v0 = momentum() / mass();
		const auto v = v0 - vf;
		const auto vnorm = v.dot(norm);
		const auto ein = std::max(energy() - momentum().dot(momentum()) * (0.5 / mass()), 0.0);
		const auto pre = (fgamma - 1.0) * ein;
		f.mass() = mass() * vnorm;
		for (int dim = 0; dim < NDIM; dim++) {
			f.momentum(dim) = momentum(dim) * vnorm + norm[dim] * pre;
		}
		f.energy() = energy() * vnorm + v0.dot(norm) * pre;
		const auto cs = std::sqrt(ein / mass() * fgamma * (fgamma - 1));
		p.second = cs + std::abs(vnorm);
		return p;
	}
};

inline std::pair<state, real> flux(const state &L, const state &R, const vect &vf, const vect &norm) {
	std::pair<state, real> F;
	const auto tmpL = L.flux(vf, norm);
	const auto tmpR = R.flux(vf, norm);
	const state &fL = tmpL.first;
	const state &fR = tmpR.first;
	const real &cL = tmpL.second;
	const real &cR = tmpR.second;
	const real c = std::max(cL, cR);
	F.first = ((fL + fR) - (R - L) * c) * 0.5;
	F.second = cR + cL;
	return F;
}

#endif /* STATE_HPP_ */
