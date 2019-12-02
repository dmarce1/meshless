/*
 * tree.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#include "tree.hpp"
#include <cmath>

#include <hpx/include/async.hpp>
#include <hpx/include/future.hpp>

#if(NDIM == 1 )
#define C 1
#define NGB 4
#else
#if( NDIM == 2 )
#define C M_PI
#define NGB 16
#else
#define C (4.0 * M_PI / 3.0)
#define NGB 32
#endif
#endif

tree::tree(std::vector<particle> &&parts) :
		level(0) {
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

void tree::compute_interactions() {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->compute_interactions();
		}
	} else {
		for (auto &part : parts) {
			std::array<vect, NDIM> E;
			for (int n = 0; n < NDIM; n++) {
				for (int m = 0; m < NDIM; m++) {
					E[n][m] = 0.0;
				}
			}
			const auto &pi = part;
			for (auto &other : part.neighbors) {
				const auto &pj = *(other.ptr);
				const auto psi_j = W(abs(pi.x - pj.x), pi.h) / pi.V;
				for (int n = 0; n < NDIM; n++) {
					for (int m = 0; m < NDIM; m++) {
						E[n][m] += (pj.x[n] - pi.x[n]) * (pj.x[m] - pi.x[m]) * psi_j;
					}
				}
			}
			const auto B = matrix_inverse(E);
			for (auto &other : part.neighbors) {
				const auto &pj = *(other.ptr);
				const auto psi_j = W(abs(pi.x - pj.x), pi.h) / pi.V;
				vect psi_a_j;
				for (int n = 0; n < NDIM; n++) {
					psi_a_j[n] = 0.0;
					for (int m = 0; m < NDIM; m++) {
						psi_a_j[n] += B[n][m] * (pj.x[m] - pi.x[m]) * psi_j;
					}
				}
			}
		}
	}
}

void tree::find_neighbors() {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->find_neighbors();
		}
	} else {
		for (auto &part : parts) {
			for (auto this_n : part.neighbors) {
				auto &nset = this_n.ptr->neighbors;
				const auto ptr = &part;
				if (nset.find(ptr) == nset.end()) {
					nset.insert(ptr);
				}
			}
		}
	}
}

tree_stats tree::compute_tree_statistics() const {
	tree_stats st;
	if (refined) {
		st.n_nodes = 1;
		st.n_leaves = 0;
		st.max_level = 0;
		st.n_part = 0;
		st.n_neighbor = 0;
		for (int ci = 0; ci < NCHILD; ci++) {
			auto tmp = children[ci]->compute_tree_statistics();
			st.n_nodes += tmp.n_nodes;
			st.n_leaves += tmp.n_leaves;
			st.max_level = std::max(st.max_level, tmp.max_level);
			st.n_part += tmp.n_part;
			st.n_neighbor += tmp.n_neighbor;
		}
	} else {
		st.n_nodes = 1;
		st.n_leaves = 1;
		st.max_level = level;
		st.n_part = parts.size();
		st.n_neighbor = 0;
		for (const auto &part : parts) {
			st.n_neighbor += part.neighbors.size();
		}
	}
	return st;
}

tree::tree(std::vector<particle> &&parts, const range &_box) :
		box(_box) {
	make_tree(std::move(parts));
}

void tree::make_tree(std::vector<particle> &&_parts) {
	if (_parts.size() <= NMAX) {
		refined = false;
		parts = std::move(_parts);
	} else {
		refined = true;
		std::vector<particle*> ptrs(_parts.size());
		const int sz = _parts.size();
		for (int i = 0; i < sz; i++) {
			ptrs[i] = &_parts[i];
		}
		vect mid;
		for (int dim = 0; dim < NDIM; dim++) {
			std::sort(ptrs.begin(), ptrs.end(), [dim](particle *a, particle *b) {
				return a->x[dim] < b->x[dim];
			});
			mid[dim] = ptrs[sz / 2]->x[dim];
//			mid[dim] = (box.first[dim] + box.second[dim]) / 2.0;
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
			int sz = _parts.size();
			std::vector<particle> cparts;
			for (int i = 0; i < sz; i++) {
				auto &part = _parts[i];
				if (in_range(part.x, this_box)) {
					cparts.push_back(std::move(part));
					part = std::move(_parts[sz - 1]);
					sz--;
					i--;
					_parts.resize(sz);
				}
			}
			children[n] = new_tree(std::move(cparts), this_box);
		}
	}

}

void tree::form_tree() {
	form_tree(nullptr, 0);
}

void tree::form_tree(tree_ptr _parent, int _level) {
	level = _level;
	parent = _parent;
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->form_tree(self, level + 1);
		}
	}
}

void tree::compute_smoothing_lengths() {
	compute_smoothing_lengths(self);
}

void tree::compute_smoothing_lengths(tree_ptr root) {
	if (refined) {
		std::vector<hpx::future<void>> futs;
		for (int ci = 0; ci < NCHILD; ci++) {
			futs.push_back(hpx::async([this, ci, root]() {
				children[ci]->compute_smoothing_lengths(root);
			}));
		}
		hpx::wait_all(futs);
	} else {
		std::vector<const particle*> others;
		const real vol = box_volume(box);
		const real h0 = std::pow(NGB * vol / parts.size(), 1.0 / NDIM);
		for (auto &part : parts) {
			if (part.h_not_set()) {
				part.h = h0;
			}
		}
		const int sz = parts.size();
		for (int i = 0; i < sz; i++) {
			auto &part = parts[i];
			bool done = false;
			double N;
			int iters = 0;
			//		printf("\n");
			double max_dh = std::numeric_limits<double>::max();
			do {
				root->particles_in_sphere(others, part.x, part.h);
				N = 0.0;
				auto &h = part.h;
				double dNdh = 0.0;
				for (auto &other : others) {
					if (other != &part) {
						const auto d = distance(part.x, other->x);
						N += C * std::pow(h, NDIM) * W(d, h);
						dNdh += C * NDIM * std::pow(h, NDIM - 1) * W(d, h);
						dNdh += C * std::pow(h, NDIM) * dW_dh(d, h);
					}
				}
				//		printf("%e %e %e\n", h, N, max_dh);
				if (dNdh == 0.0) {
					h *= 2.0;
				} else {
					auto dh = -(N - NGB) / dNdh;
					dh = std::min(h, std::max(-h / 2.0, dh));
					max_dh = std::min(0.99 * max_dh, std::abs(dh));
					dh = std::copysign(std::min(max_dh, std::abs(dh)), dh);
					h += 0.99 * dh;
				}
				if (iters > 1000) {
					printf("Failed to converge\n");
					abort();
				}
				iters++;
			} while (std::abs(NGB - N) > 0.5);
			part.neighbors.clear();
			for (auto &other : others) {
				part.neighbors.insert(const_cast<particle*>(other));
			}
		}
	}
}

real tree::compute_volumes() {
	real total = 0.0;
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			total += children[ci]->compute_volumes();
		}
	} else {
		for (auto &part : parts) {
			double omega = 0.0;
			for (const auto &neighbor : part.neighbors) {
				omega += W(abs(part.x - neighbor.ptr->x), part.h);
			}
			//	printf( "%e %e\n", 1.0 / omega, part.h);
			part.V = 1.0 / omega;
			total += part.V;
		}
	}
	return total;
}

void tree::particles_in_range(std::vector<const particle*> &these_parts, const range &R) const {
	if (ranges_intersect(R, box)) {
		if (refined) {
			for (int ci = 0; ci < NCHILD; ci++) {
				children[ci]->particles_in_range(these_parts, R);
			}
		} else {
			for (const auto &part : parts) {
				if (in_range(part.x, R)) {
					these_parts.push_back(&part);
				}
			}
		}
	}
}

void tree::particles_in_sphere(std::vector<const particle*> &these_parts, const vect &c, real r) const {
	range this_box;
	these_parts.clear();
	for (int dim = 0; dim < NDIM; dim++) {
		this_box.first[dim] = c[dim] - r;
		this_box.second[dim] = c[dim] + r;
	}
	particles_in_range(these_parts, this_box);
	int sz = these_parts.size();
	for (int i = 0; i < sz; i++) {
		auto &part = these_parts[i];
		if (abs(c - part->x) > r) {
			part = these_parts[sz - 1];
			i--;
			sz--;
		}
	}
	these_parts.resize(sz);
}

void tree::destroy() {
	parent = null_tree_ptr;
	self = null_tree_ptr;
	if (refined) {
		for (int i = 0; i < NCHILD; i++) {
			children[i]->destroy();
			children[i] = null_tree_ptr;
		}
	}
}
