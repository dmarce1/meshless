/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include "tree.hpp"
#include <fenv.h>
#include <hpx/hpx_init.hpp>

int hpx_main(int argc, char *argv[]) {
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);

	auto parts = cartesian_particle_set(64);
	auto t = tree::new_tree(std::move(parts));
	t->form_tree();
	t->compute_smoothing_lengths();
	auto stats = t->compute_tree_statistics();
	t->find_neighbors();
	printf("Effective volume = %e\n", t->compute_volumes());
	printf("nodes              = %i\n", stats.n_nodes);
	printf("leaves             = %i\n", stats.n_leaves);
	printf("n_parts            = %i\n", stats.n_part);
	printf("n_parts/leaf       = %e\n", stats.n_part / double(stats.n_leaves));
	printf("n_neighbor/n_parts = %e\n", stats.n_neighbor / double(stats.n_part));
	printf("max_level          = %i\n", stats.max_level);
	return hpx::finalize();

}

int main(int argc, char *argv[]) {
	hpx::init(argc, argv);
}
