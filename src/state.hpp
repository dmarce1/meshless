/*
 * state.hpp
 *
 *  Created on: Dec 2, 2019
 *      Author: dmarce1
 */

#ifndef STATE_HPP_
#define STATE_HPP_


#include <cassert>

#include "vect.hpp"

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
	std::pair<state, real> flux(const vect &vf, const vect &norm) const;
	state to_prim() const;
	state to_con() const;
};


std::pair<state, real> flux(const state &L, const state &R, const vect &vf, const vect &norm);

#endif /* STATE_HPP_ */
