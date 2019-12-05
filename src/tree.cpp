/*
 * tree.cpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#include "tree.hpp"
#include <cmath>
#include <cstdlib>

#define Ncond_crit (100)

#if(NDIM == 1 )
#define C 1
#define NGB 2
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
	for (int dim = 0; dim < NDIM; dim++) {
		box.first[dim] -= 1.0e-12 * std::abs(box.first[dim]);
		box.second[dim] += 1.0e-12 * std::abs(box.second[dim]);
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
			auto &pi = part;
			for (auto &other : part.neighbors) {
				const auto &pj = *(other->ptr);
				const auto psi_j = W(abs(pi.x - pj.x), pi.h) * pi.V;
				for (int n = 0; n < NDIM; n++) {
					for (int m = 0; m < NDIM; m++) {
						E[n][m] += (pj.x[n] - pi.x[n]) * (pj.x[m] - pi.x[m]) * psi_j;
					}
				}
			}
			decltype(E) B;
			pi.Ncond = condition_number(E, B);
			for (auto &other : part.neighbors) {
				const auto &pj = *(other->ptr);
				const auto psi_j = W(abs(pi.x - pj.x), pi.h) * pi.V;
				real psi_a_j;
				for (int n = 0; n < NDIM; n++) {
					psi_a_j = 0.0;
					for (int m = 0; m < NDIM; m++) {
						psi_a_j += B[n][m] * (pj.x[m] - pi.x[m]) * psi_j;
					}
					other->psi_a[n] = psi_a_j;
					const real da = part.V * psi_a_j;
					//			printf( "%e %e\n", part.V, psi_j);
					other->area[n] += da;
					other->ret->area[n] -= da;
				}
			}
		}
	}
}

void tree::boundary_conditions() {
	boundary_conditions(box);
}

void tree::boundary_conditions(const range &rng) {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->boundary_conditions(rng);
		}
	} else {
		int sz = parts.size();
		for (int i = 0; i < sz; i++) {
			if (!in_range(parts[i].x, rng)) {
				sz--;
				parts[i] = std::move(parts[sz]);
				i--;
			}
		}
		parts.resize(sz);
	}
}
real tree::compute_fluxes() {
	real tmin = std::numeric_limits<real>::max();
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto this_tmin = children[ci]->compute_fluxes();
			tmin = std::min(tmin, this_tmin);
		}
	} else {
		for (auto &part : parts) {
			auto &pi = part;
			for (auto &other : part.neighbors) {
				auto &pj = *(other->ptr);
				if (&pj > &pi) {
					vect A;
					for (int dim = 0; dim < NDIM; dim++) {
						A[dim] = other->area[dim];
					}
					auto dx = pj.x - pi.x;
					auto xij = pi.x + dx * (pi.h) / (pi.h + pj.h);
					auto vij = pi.v() + (pj.v() - pi.v()) * (xij - pi.x).dot(dx) / (dx.dot(dx));
					if (abs(A) != 0.0) {
						auto norm = A / abs(A);
						const auto dxL = other->xij - pi.x;
						const auto dxR = other->xij - pj.x;
						state L = other->U;
						state R = other->ret->U;
						L = L.to_con();
						R = R.to_con();
						const auto tmp = flux(L, R, vij, norm);
						other->flux = tmp.first;
						other->ret->flux = -tmp.first;
						const real vsig = tmp.second - std::min(0.0, (pi.v() - pj.v()).dot(pi.x - pj.x) / abs(pi.x - pj.x));
						tmin = std::min(tmin, pi.h / vsig);
					} else {
						for (int i = 0; i < STATE_SIZE; i++) {
							other->flux[i] = 0.0;
							other->ret->flux[i] = 0.0;
						}
					}
				}
			}
		}
	}
	return tmin;
}

void tree::output(const char *filename) {
	if (level == 0) {
		FILE *fp = fopen(filename, "wt");
		if (fp == NULL) {
			printf("System error on line %i of %s\n", __LINE__, __FILE__);
			abort();
		}
		fclose(fp);
	}
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->output(filename);
		}
	} else {
		FILE *fp = fopen("parts.txt", "at");
		for (auto &part : parts) {
			for (int dim = 0; dim < NDIM; dim++) {
				fprintf(fp, "%e ", part.x[dim]);
			}
			fprintf(fp, "%e ", part.V);
			fprintf(fp, "%e ", part.h);
			fprintf(fp, "%e ", part.rho());
			fprintf(fp, "%e ", part.E());
			for (int dim = 0; dim < NDIM; dim++) {
				fprintf(fp, "%e ", part.v()[dim]);
			}
			fprintf(fp, "%i\n", int(part.neighbors.size()));
		}
		fclose(fp);
	}
}

void tree::compute_next_step(real dt) {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->compute_next_step(dt);
		}
	} else {
		for (auto &part : parts) {
			auto &pi = part;
			pi.x = pi.x + pi.v() * dt;
			for (auto &other : part.neighbors) {
				const auto area = abs(other->area);
				pi.st = pi.st - other->flux * dt * area;
			}
		}
	}
}

void tree::initialize(const std::function<state(vect)> &f) {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->initialize(f);
		}
	} else {
		for (auto &part : parts) {
			part.st = f(part.x) * part.V;
		}
	}
}
std::vector<particle> tree::gather_particles() const {
	std::vector<particle> return_parts;
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			const auto tmp = children[ci]->gather_particles();
			return_parts.insert(return_parts.end(), tmp.begin(), tmp.end());
		}
	} else {
		return_parts = parts;
	}
	return return_parts;
}

void tree::find_neighbors() {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->find_neighbors();
		}
	} else {
		for (auto &part : parts) {
			for (auto &this_n : part.neighbors) {
				auto &others = this_n->ptr->neighbors;
				const auto ptr = &part;
				bool found = false;
				for (auto &other : others) {
					if (ptr == other->ptr) {
						found = true;
						this_n->ret = other;
						other->ret = this_n;
						break;
					}
				}
				if (!found) {
					others.push_back(std::make_shared<neighbor>(ptr));
					auto &other = others[others.size() - 1];
					this_n->ret = other;
					other->ret = this_n;
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
			if (sz % 2 == 0) {
				mid[dim] = (ptrs[sz / 2]->x[dim] + ptrs[sz / 2 - 1]->x[dim]) * 0.5;
			} else {
				mid[dim] = ptrs[sz / 2]->x[dim];
			}
		}
		range this_box;
		for (int ci = 0; ci < NCHILD; ci++) {
			int m = ci;
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
			children[ci] = new_tree(std::move(cparts), this_box);
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
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->compute_smoothing_lengths(root);
		}
	} else if (parts.size()) {
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
			double N;
			int iters = 0;
			//		printf("\n");
			double max_dh = std::numeric_limits<double>::max();
			bool done = false;
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
				if (std::abs(NGB - N) < 1e-6) {
					done = true;
					root->particles_in_sphere(others, part.x, h);
				}
			} while (!done);
			part.neighbors.resize(0);
			for (auto &other : others) {
				if (other != &part) {
					part.neighbors.push_back(std::make_shared<neighbor>(const_cast<particle*>(other)));
				}
			}
		}
	}
}

void tree::compute_gradients() {
	if (refined) {
		for (int ci = 0; ci < NCHILD; ci++) {
			children[ci]->compute_gradients();
		}
	} else if (parts.size()) {

		std::array<state, NDIM> gradient;

		for (auto &part : parts) {
			auto &pi = part;
			for (int dim = 0; dim < NDIM; dim++) {
				for (int i = 0; i < STATE_SIZE; i++) {
					gradient[dim][i] = 0.0;
				}
			}
			state max_ngb;
			state min_ngb;
			const state st_i = state(pi.st / pi.V).to_prim();
			for (int i = 0; i < STATE_SIZE; i++) {
				max_ngb[i] = min_ngb[i] = st_i[i];
			}
			for (auto &other : part.neighbors) {
				const auto &pj = *(other->ptr);
				const state st_j = state(pj.st / pj.V).to_prim();
				for (int dim = 0; dim < NDIM; dim++) {
					gradient[dim] = gradient[dim] + (st_j - st_i) * other->psi_a[dim];
				}
			}
			real max_dx = 0.0;
			for (auto &other : part.neighbors) {
				const auto &pj = *(other->ptr);
				auto dx = pj.x - pi.x;
				auto &xij = other->xij;
				xij = pi.x + dx * (pi.h) / (pi.h + pj.h);
				const state st_j = state(pj.st / pj.V).to_prim();
				const auto mid_dx = other->xij - pi.x;
				max_dx = std::max(max_dx, std::sqrt(mid_dx.dot(mid_dx)));
				max_ngb = max(max_ngb, st_j);
				min_ngb = min(min_ngb, st_j);
			}
			const auto beta = std::max(1.0, std::min(2.0, Ncond_crit / pi.Ncond));
			real alpha;
			for (int i = 0; i < STATE_SIZE; i++) {
				const auto dmax_ngb = max_ngb[i] - st_i[i];
				const auto dmin_ngb = st_i[i] - min_ngb[i];
				real grad_abs = 0.0;
				for (int dim = 0; dim < NDIM; dim++) {
					grad_abs += gradient[dim][i] * gradient[dim][i];
				}
				grad_abs = std::sqrt(grad_abs);
				const real den = grad_abs * max_dx;
				if (den == 0.0) {
					alpha = 1.0;
				} else {
					alpha = std::min(1.0, beta * std::min(dmax_ngb / den, dmin_ngb / den));
				}
				for (auto &other : part.neighbors) {
					const auto &pj = *(other->ptr);
					for (int dim = 0; dim < NDIM; dim++) {
						gradient[dim][i] = gradient[dim][i] * alpha;
					}
				}

			}
			constexpr auto psi1 = 0.5;
			constexpr auto psi2 = 0.25;
			for (auto &other : part.neighbors) {
				const auto &pj = *(other->ptr);
				state &phi_mid = other->U;
				state phi_i = state(pi.st / pi.V).to_prim();
				state phi_j = state(pj.st / pj.V).to_prim();
				phi_mid = phi_i;
				for (int dim = 0; dim < NDIM; dim++) {
					const auto dx = other->xij - pi.x;
					phi_mid = phi_mid + gradient[dim] * dx[dim];
				}
				const auto dphi_abs = abs(phi_j, phi_i);
				const auto delta_1 = dphi_abs * psi1;
				const auto delta_2 = dphi_abs * psi2;
				const auto phi_min = min(phi_i, phi_j);
				const auto phi_max = max(phi_i, phi_j);
				const auto phi_bar = phi_i + (phi_j - phi_i) * abs(other->xij - pi.x) / abs(pi.x - pj.x);
				for (int i = 0; i < STATE_SIZE; i++) {
					if (phi_i[i] == phi_j[i]) {
						phi_mid[i] = phi_i[i];
					} else if (phi_i[i] < phi_j[i]) {
						state phi_m;
						if ((phi_min[i] - delta_1[i]) * phi_min[i] > 0.0) {
							phi_m[i] = phi_min[i] - delta_1[i];
						} else {
							phi_m[i] = phi_min[i] * std::abs(phi_min[i]) / (std::abs(phi_min[i]) + delta_1[i]);
						}
						phi_mid[i] = std::max(phi_m[i], std::min(phi_bar[i] + delta_2[i], phi_mid[i]));
					} else {
						state phi_p;
						if ((phi_max[i] + delta_1[i]) * phi_max[i] > 0.0) {
							phi_p[i] = phi_max[i] + delta_1[i];
						} else {
							phi_p[i] = phi_max[i] * std::abs(phi_max[i]) / (std::abs(phi_max[i]) + delta_1[i]);
						}
						phi_mid[i] = std::max(phi_p[i], std::min(phi_bar[i] - delta_2[i], phi_mid[i]));
					}
				}
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
			for (auto &neighbor : part.neighbors) {
				omega += W(abs(part.x - neighbor->ptr->x), part.h);
				for (int dim = 0; dim < NDIM; dim++) {
					neighbor->area[dim] = 0.0;
				}
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
