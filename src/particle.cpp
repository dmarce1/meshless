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

std::vector<particle> random_particle_set(int N) {
	std::vector<particle> parts(N);
	for (auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] = rand_unit_box();
		}
	}
	return std::move(parts);
}

std::vector<particle> cartesian_particle_set(int N) {
	std::vector<particle> parts;
	parts.reserve(std::pow(N, NDIM));
	particle part;
#if(NDIM==3)
	for (int l = 0; l < N; l++) {
		part.x[2] = (l + 0.5) / N - 0.5;
#endif
#if(NDIM>=2)
		for (int j = 0; j < N; j++) {
			part.x[1] = (j + 0.5) / N - 0.5;
#endif
			for (int i = 0; i < N; i++) {
				part.x[0] = (i + 0.5) / N - 0.5;
				parts.push_back(part);
			}
#if(NDIM>=2)
		}
#endif
#if(NDIM==3)
	}
#endif
		return std::move(parts);
	}
