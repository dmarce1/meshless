/*
 * state.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include "state.hpp"

state state::to_prim() const {
	state p;
	p.mass() = mass();
	p.energy() = energy() - momentum().dot(momentum()) / mass() / 2.0;
	for (int dim = 0; dim < NDIM; dim++) {
		p.momentum(dim) = momentum(dim) / mass();
	}
	return p;
}

state state::to_con() const {
	state c;
	c.energy() = energy() + momentum().dot(momentum()) * mass() / 2.0;
	for (int dim = 0; dim < NDIM; dim++) {
		c.momentum(dim) = momentum(dim) * mass();
	}
	c.mass() = mass();
	return c;
}


std::pair<state, real> state::flux(const vect &vf, const vect &norm) const {
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


std::pair<state, real> flux(const state &L, const state &R, const vect &vf, const vect &norm) {
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
