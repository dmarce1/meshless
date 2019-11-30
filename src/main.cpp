/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */


//#include <array>
//#include <cmath>
//#include <limits>
//#include <memory>
//#include <vector>
//
//template<int NDIM>
//struct particle {
//	double h;
//	double m;
//	vect<NDIM> x;
//	vect<NDIM> v;
//	std::vector<particle*> neighbors;
//};
//
//template<int NDIM>
//class tree;
//
//template<int NDIM>
//using tree_ptr = std::shared_ptr<tree<NDIM>>;
//
//template<int NDIM, class ...Args>
//auto new_tree(Args &&... args) {
//	return std::make_shared<tree<NDIM>>(std::forward<Args>(args)...);
//}
//
//template<int NDIM>
//class tree {
//	static constexpr int NCHILD = 1 << NDIM;
//	static constexpr int MAXPARTS = 1000;
//	std::vector<particle<NDIM>> particles;
//	tree_ptr<NDIM> parent;
//	std::array<tree_ptr<NDIM>, NCHILD> children;
//	vect<NDIM> min;
//	vect<NDIM> max;
//	bool refined;
//
//	tree() :
//			refined(false) {
//	}
//
//public:
//
//	tree(std::vector<particle<NDIM>> &&parts) {
//		for (int dim = 0; dim < NDIM; dim++) {
//			min[dim] = +std::numeric_limits<double>::max();
//			max[dim] = -std::numeric_limits<double>::max();
//		}
//		for (auto part : parts) {
//			for (int dim = 0; dim < NDIM; dim++) {
//				min[dim] = std::min(min[dim], part.x[dim]);
//				max[dim] = std::max(max[dim], part.x[dim]);
//			}
//		}
//		if (parts.size() > MAXPARTS) {
//			refined = true;
//			for (int n = 0; n < NCHILD; n++) {
//				int m = n;
//				vect<NDIM> this_min;
//				vect<NDIM> this_max;
//				for (int dim = 0; dim < NDIM; dim++) {
//					const auto mid = 0.5 * (min[dim] + max[dim]);
//					if (m & 1) {
//						this_min[dim] = mid;
//						this_max[dim] = max[dim];
//					} else {
//						this_min[dim] = min[dim];
//						this_max[dim] = mid;
//					}
//					m >>= 1;
//				}
//				int sz = particles.size();
//				std::vector<particle<NDIM>> child_parts;
//				for (int i = 0; i < sz; i++) {
//					auto &part = particles[i];
//					if (in_range(part.x, this_min, this_max)) {
//						child_parts.push_back(std::move(part));
//						std::swap(part, particles[sz - 1]);
//						sz--;
//					}
//				}
//				children[n] = new_tree<NDIM>(std::move(child_parts));
//			}
//		} else {
//			particles = std::move(parts);
//			refined = false;
//		}
//	}
//
//};
//
//template<int NDIM>
//auto compute_smoothing_lengths(std::vector<particle<NDIM>> &parts) {
//	const int sz = parts.size();
//	constexpr auto C = NDIM == 1 ? 1.0 : (NDIM == 2 ? M_PI : 4.0 * M_PI / 3.0);
//	constexpr auto Nngb = NDIM == 1 ? 8 : (NDIM == 2 ? 16 : 32);
//	for (int i = 0; i < sz; i++) {
//		auto &part = parts[i];
//		auto &h = part.h;
//		double mind = std::numeric_limits<double>::max();
//		for (int j = 0; j < sz; j++) {
//			if (i != j) {
//				const auto d = distance<NDIM>(part.x, parts[j].x);
//				mind = std::min(mind, d);
//			}
//		}
//		part.neighbors.clear();
//		h = 2.0 * mind;
//		for (int j = 0; j < sz; j++) {
//			if (i != j) {
//				const auto d = distance<NDIM>(part.x, parts[j].x);
//				if (d < h) {
//					part.neighbors.push_back(&parts[j]);
//				}
//			}
//		}
//		double N = 0.0;
//		int iters = 0;
//		do {
//			N = 0.0;
//			double dNdh = 0.0;
//			for (const auto *neighbor : part.neighbors) {
//				const auto d = distance<NDIM>(part.x, neighbor->x);
//				N += C * std::pow(h, NDIM) * W<NDIM>(d, h);
//				dNdh += C * NDIM * std::pow(h, NDIM - 1) * W<NDIM>(d, h);
//				dNdh += C * std::pow(h, NDIM) * dW_dh<NDIM>(d, h);
//			}
//			const auto dh = 0.5 * -(N - Nngb) / dNdh;
//			h += dh;
//			part.neighbors.clear();
//			for (int j = 0; j < sz; j++) {
//				if (i != j) {
//					const auto d = distance<NDIM>(part.x, parts[j].x);
//					if (d < h) {
//						part.neighbors.push_back(&parts[j]);
//					}
//				}
//			}
//			iters++;
//			if (iters > 250) {
//				printf("Error - could not find smoothing length\n");
//				abort();
//			}
//		} while (std::abs(N - Nngb) > 0.5 || part.neighbors.size() == 0);
//	}
//}
//
//template<int NDIM>
//auto compute_avg_Nngb(const std::vector<particle<NDIM>> &parts) {
//	const int sz = parts.size();
//	double sum = 0.0;
//	for (int i = 0; i < sz; i++) {
//		const int count = parts[i].neighbors.size();
//		sum += count;
//		printf("%i\n", count);
//	}
//	sum /= sz;
//	printf("%f\n", sum);
//	return sum;
//}
//
//template<int NDIM>
//auto random_particles(int N) {
//	std::vector<particle<NDIM>> parts(N);
//	for (auto &part : parts) {
//		for (int dim = 0; dim < NDIM; dim++) {
//			part.x[dim] = rand_unit_box();
//		}
//	}
//	return std::move(parts);
//}

int main() {

}

