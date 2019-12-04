/*
 * main.cpp
 *
 *  Created on: Nov 28, 2019
 *      Author: dmarce1
 */

#include "tree.hpp"
#include <fenv.h>
#include <hpx/hpx_init.hpp>

state sod(const vect &x) {
	state U;
	for (int i = 0; i < STATE_SIZE; i++) {
		U[i] = 0.0;
	}
	if (x[0] > 0.0) {
		U.mass() = 1.0;
		U.energy() = 2.5;
	} else {
		U.mass() = 0.125;
		U.energy() = 0.25;
	}
	return U;
}

int hpx_main(int argc, char *argv[]) {
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);

	auto parts = cartesian_particle_set(10000);
	auto t = tree::new_tree(std::move(parts));
	t->form_tree();
	t->compute_smoothing_lengths();
	t->compute_volumes();
	t->find_neighbors();
	t->initialize(sod);
	double tm = 0.0;
	real dt = 0.0;
	int iter = 0;
	while (tm < 0.1) {
		t->compute_interactions();
		real dt = 0.4 * t->compute_fluxes();
		t->compute_next_step(dt);
		parts = t->gather_particles();
		t = tree::new_tree(std::move(parts));
		t->form_tree();
		t->compute_smoothing_lengths();
		t->find_neighbors();
		t->compute_volumes();
		tm += dt;
		printf( "%e %e\n", tm, dt);
		iter++;
	//	break;
	}
	parts = t->gather_particles();
	t = tree::new_tree(std::move(parts));
	t->form_tree();
	t->compute_smoothing_lengths();
	t->compute_volumes();
	auto stats = t->compute_tree_statistics();
	t->output("parts.txt");
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
