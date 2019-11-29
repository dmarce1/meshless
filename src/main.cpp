/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include <array>
#include <cmath>
#include <limits>
#include <vector>

auto bspline(double r) {
	if (r < 1) {
		const auto r2 = r * r;
		return 1.0 - 1.5 * r2 + 0.75 * r2 * r;
	} else if (r <= 2) {
		return 0.25 * std::pow(2 - r, 3);
	} else {
		return 0.0;
	}
}

auto d_bspline_dr(double r) {
	if (r < 1) {
		return -3.0 * r + 2.25 * r * r;
	} else if (r <= 2) {
		return -0.75 * std::pow(2 - r, 2);
	} else {
		return 0.0;
	}
}

template<int NDIM>
auto W_norm(double h) {
	if constexpr (NDIM == 1) {
		return 3.0 * h / 8.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 7 * h * h / 80.0;
	} else {
		return 2 * 2.0 * M_PI * h * h * h / 32.0;
	}
}

template<int NDIM>
auto dW_norm_dh(double h) {
	if constexpr (NDIM == 1) {
		return 3.0 / 8.0;
	} else if constexpr (NDIM == 2) {
		return 2 * M_PI * 15 * h / 80.0;
	} else {
		return 2 * 6.0 * M_PI * h * h / 32.0;
	}
}

template<int NDIM>
auto W(double r, double h) {
	auto norm = W_norm<NDIM>(h);
	return bspline(2 * r / h) / norm;
}

template<int NDIM>
auto dW_dh(double r, double h) {
	auto norm = W_norm<NDIM>(h);
	return -2 * r * d_bspline_dr(2 * r / h) / (h * h) / norm - dW_norm_dh<NDIM>(h) * W<NDIM>(r, h) / norm;
}

template<int NDIM>
using vect = std::array<double,NDIM>;

template<int NDIM>
struct particle {
	double h;
	double m;
	vect<NDIM> x;
	vect<NDIM> v;
	std::vector<particle*> neighbors;
};

template<int NDIM>
auto distance(const vect<NDIM> &a, const vect<NDIM> &b) {
	double d = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		d += std::pow(a[dim] - b[dim], 2);
	}
	return std::sqrt(d);
}

template<int NDIM>
auto compute_smoothing_lengths(std::vector<particle<NDIM>> &parts) {
	const int sz = parts.size();
	constexpr auto C = NDIM == 1 ? 1.0 : (NDIM == 2 ? M_PI : 4.0 * M_PI / 3.0);
	constexpr auto Nngb = NDIM == 1 ? 8 : (NDIM == 2 ? 16 : 32);
	for (int i = 0; i < sz; i++) {
		auto &part = parts[i];
		auto &h = part.h;
		double mind = std::numeric_limits<double>::max();
		for (int j = 0; j < sz; j++) {
			if (i != j) {
				const auto d = distance<NDIM>(part.x, parts[j].x);
				mind = std::min(mind, d);
			}
		}
		part.neighbors.clear();
		h = 2.0 * mind;
		for (int j = 0; j < sz; j++) {
			if (i != j) {
				const auto d = distance<NDIM>(part.x, parts[j].x);
				if (d < h) {
					part.neighbors.push_back(&parts[j]);
				}
			}
		}
		double N = 0.0;
		int iters = 0;
		do {
			N = 0.0;
			double dNdh = 0.0;
			for (const auto *neighbor : part.neighbors) {
				const auto d = distance<NDIM>(part.x, neighbor->x);
				N += C * std::pow(h, NDIM) * W<NDIM>(d, h);
				dNdh += C * NDIM * std::pow(h, NDIM - 1) * W<NDIM>(d, h);
				dNdh += C * std::pow(h, NDIM) * dW_dh<NDIM>(d, h);
			}
			const auto dh = 0.5 * -(N - Nngb) / dNdh;
			h += dh;
			part.neighbors.clear();
			for (int j = 0; j < sz; j++) {
				if (i != j) {
					const auto d = distance<NDIM>(part.x, parts[j].x);
					if (d < h) {
						part.neighbors.push_back(&parts[j]);
					}
				}
			}
			iters++;
			if (iters > 250) {
				printf("Error - could not find smoothing length\n");
				abort();
			}
		} while (std::abs(N - Nngb) > 0.5 || part.neighbors.size() == 0);
	}
}

template<int NDIM>
auto compute_avg_Nngb(const std::vector<particle<NDIM>> &parts) {
	const int sz = parts.size();
	double sum = 0.0;
	for (int i = 0; i < sz; i++) {
		const int count = parts[i].neighbors.size();
		sum += count;
		printf("%i\n", count);
	}
	sum /= sz;
	printf("%f\n", sum);
	return sum;
}

auto rand_unit_box() {
	const auto x = (rand() + 0.5) / (double(RAND_MAX) + 1.0);
	return 2.0 * (x - 0.5);
}

template<int NDIM>
auto random_particles(int N) {
	std::vector<particle<NDIM>> parts(N);
	for (auto &part : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			part.x[dim] = rand_unit_box();
		}
	}
	return std::move(parts);
}

int main() {
	constexpr auto ndim = 3;
	auto parts = random_particles<ndim>(10000);
	compute_smoothing_lengths<ndim>(parts);
	compute_avg_Nngb<ndim>(parts);
	return 0;

}

