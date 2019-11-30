/*
 * particle.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#include "particle.hpp"

particle::particle() {
	h = -1.0;
}

bool particle::h_not_set() const {

	return h <= 0;
}
