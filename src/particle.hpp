/*
 * particle.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include "math.hpp"
#include "state.hpp"
#include <vector>

struct particle;

struct neighbor {
	particle *ptr;
	neighbor *ret;
	atomic_vect area;
	state flux;
	inline neighbor(particle *_ptr) :
			ptr(_ptr) {
	}
	inline bool operator<(const neighbor &other) const {
		return ptr < other.ptr;
	}
};

struct particle {
	vect x;
	real h;
	real V;
	state st;
	real rho() const {
		return st.mass() / V;
	}
	real E() const {
		return st.energy() / V;
	}
	vect v() const {
		return st.momentum() / st.mass();
	}
	std::vector<neighbor> neighbors;
	particle();
	bool h_not_set() const;
	real omega() const;
};

std::vector<particle> random_particle_set(int);
std::vector<particle> cartesian_particle_set(int);

inline bool particle::h_not_set() const {
	return h <= 0;
}

#endif /* PARTICLE_HPP_ */
