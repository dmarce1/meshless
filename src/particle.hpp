/*
 * particle.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include "math.hpp"
#include <set>
#include <vector>

struct particle;

struct neighbor {
	particle *ptr;
	vect psi_a;
	inline neighbor(particle *_ptr) :
			ptr(_ptr) {
	}
	inline bool operator<(const neighbor& other) const {
		return ptr < other.ptr;
	}
};

struct particle {
	vect x;
	real h;
	real V;
	std::set<neighbor> neighbors;
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
