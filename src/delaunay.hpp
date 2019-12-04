/*
 * delauny.hpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#ifndef DELAUNAY_HPP_
#define DELAUNAY_HPP_

#include "particle.hpp"
#include "vect.hpp"

using delauny_region = general_vect<int,NDIM+1>;

std::vector<delauny_region> compute_delaunay_regions(const std::vector<particle>&);

#endif /* DELAUNAY_HPP_ */
