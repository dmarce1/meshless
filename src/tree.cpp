/*
 * tree.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#include "tree.hpp"
#include <cmath>

#if(NDIM == 1 )
#define C 1
#define NGB 8
#else
#if( NDIM == 2 )
#define C M_PI
#define NGB 16
#else
#define C (4.0 * M_PI / 3.0)
#define NGB 32
#endif
#endif

tree::tree(std::vector<particle> &&parts) {
	for (int dim = 0; dim < NDIM; dim++) {
		box.first[dim] = +std::numeric_limits<real>::max();
		box.second[dim] = -std::numeric_limits<real>::max();
	}
	for (auto part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			const auto &x = part.x[dim];
			auto &min = box.first[dim];
			auto &max = box.second[dim];
			min = std::min(min, x);
			max = std::max(max, x);
		}
	}
	make_tree(std::move(parts));
}
tree::tree(std::vector<particle> &&parts, const range &_box, tree_ptr _parent) :
		parent(_parent), box(_box) {
	make_tree(std::move(parts));
}
void tree::make_tree(std::vector<particle> &&parts) {
	if (parts.size() <= NMAX) {
		refined = false;
	} else {
		std::vector<particle*> ptrs(parts.size());
		const int sz = parts.size();
		for (int i = 0; i < sz; i++) {
			ptrs[i] = &parts[i];
		}
		vect mid;
		for (int dim = 0; dim < NDIM; dim++) {
			std::sort(ptrs.begin(), ptrs.end(), [dim](particle *a, particle *b) {
				return a->x[dim] < b->x[dim];
			});
			mid[dim] = (ptrs[sz / 2]->x[dim] + ptrs[sz / 2 + 1]->x[dim]) * 0.5;
		}
		range this_box;
		for (int n = 0; n < NCHILD; n++) {
			int m = n;
			for (int dim = 0; dim < NDIM; dim++) {
				if (m & 1) {
					this_box.first[dim] = mid[dim];
					this_box.second[dim] = box.second[dim];
				} else {
					this_box.first[dim] = box.first[dim];
					this_box.second[dim] = mid[dim];
				}
				m >>= 1;
			}
			int sz = parts.size();
			std::vector<particle> cparts;
			for (int i = 0; i < sz; i++) {
				auto &part = parts[i];
				if (in_range(part.x, this_box)) {
					cparts.push_back(std::move(part));
					part = std::move(parts[sz - 1]);
					sz--;
				}
			}
			children[n] = new_tree(std::move(cparts), this_box, self);
		}
	}

}

void tree::compute_smoothing_lengths() {
	compute_smoothing_lengths(self);
}

void tree::compute_smoothing_lengths(tree_ptr root) {
	const real h0 = max_span(box) / 2.0;
	for (auto &part : parts) {
		if (part.h_not_set()) {
			part.h = h0;
		}
	}

	/**********/
	bool done = false;
	std::vector<const particle*> others;
	do {
		range R = box;
		for (auto &part : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				R.first[dim] = std::min(R.first[dim], part.x[dim] - part.h);
				R.second[dim] = std::min(R.second[dim], part.x[dim] + part.h);
			}
		}
		done = true;
		others = particles_in_range(R);
		for (auto &part : parts) {
			double N = 0.0;
			double dNdh = 0.0;
			auto &h = part.h;
			for (auto &other : others) {
				if (other != &part) {
					const auto d = distance(part.x, other->x);
					N += C * std::pow(h, NDIM) * W(d, h);
					dNdh += C * NDIM * std::pow(h, NDIM - 1) * W(d, h);
					dNdh += C * std::pow(h, NDIM) * dW_dh(d, h);
				}
			}
			const auto dh = 0.5 * -(N - NGB) / dNdh;
			h += dh;
			if (std::abs(NGB - N) > 0.5) {
				done = false;
			}
		}
	} while (!done);
	for( auto& part : parts) {
		part.neighbors.clear();
		for( auto& other : others ) {
			if( distance(other->x, part.x) < part.h) {
				part.neighbors.push_back(other);
			}
		}
	}
 }

void tree::find_siblings() {
	std::array<tree_ptr, NSIBLING> null_sibs;
	find_siblings(std::move(null_sibs));
}

void tree::find_siblings(std::array<tree_ptr, NSIBLING> &&my_sibs) {
	siblings = std::move(my_sibs);
	if (refined) {
		std::array<tree_ptr, NSIBLING> child_sibs;
		for (int ci = 0; ci < NCHILD; ci++) {
			for (int dim = 0; dim < NDIM; dim++) {
				const int mask = 1 << dim;
				const int right = 2 * dim + 1;
				const int left = 2 * dim;
				if (mask & ci) {
					if (siblings[right] != nullptr) {
						child_sibs[right] = siblings[right]->siblings[left];
					} else {
						child_sibs[right] = nullptr;
					}
					child_sibs[left] = children[ci ^ mask];
				} else {
					if (siblings[left] != nullptr) {
						child_sibs[left] = siblings[left]->siblings[right];
					} else {
						child_sibs[left] = nullptr;
					}
					child_sibs[right] = children[ci ^ mask];
				}
			}
			children[ci]->find_siblings(std::move(child_sibs));
		}
	}
}

std::vector<const particle*> tree::particles_in_range(const range &R) const {
	std::vector<const particle*> these_parts;
	if (ranges_intersect(R, box)) {
		if (refined) {
			for (int ci = 0; ci < NCHILD; ci++) {
				const auto tmp = children[ci]->particles_in_range(R);
				these_parts.insert(these_parts.end(), tmp.cbegin(), tmp.cend());
			}
		} else {
			for (const auto &part : parts) {
				if (in_range(part.x, R)) {
					these_parts.push_back(&part);
				}
			}
		}
	}
	return std::move(these_parts);
}

void tree::destroy() {
	parent = null_tree_ptr;
	self = null_tree_ptr;
	for (int i = 0; i < NSIBLING; i++) {
		siblings[i] = null_tree_ptr;
	}
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->destroy();
			children[i] = null_tree_ptr;
		}
	}
}
