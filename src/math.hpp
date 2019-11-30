/*
 * math.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef MATH_HPP_
#define MATH_HPP_

#include "dim.hpp"
#include "real.hpp"

#include <array>

using vect = std::array<real,NDIM>;
using range = std::pair<vect,vect>;

real rand_unit_box();
real W(real, real);
real dW_dh(real, real);
real distance(const vect&, const vect&);
bool in_range(const vect&, const range&);
real max_span(const range&);
bool ranges_intersect(const range&, const range&);

#endif /* MATH_HPP_ */
