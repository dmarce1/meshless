/*
 * particle.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include "math.hpp"
#include <vector>

struct particle {
	vect x;
	real h;
	std::vector<const particle*> neighbors;
	particle();
	bool h_not_set() const;
};

#endif /* PARTICLE_HPP_ */
