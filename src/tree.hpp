/*
 * tree.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef TREE_HPP_
#define TREE_HPP_

#include "particle.hpp"

#include <algorithm>
#include <limits>
#include <functional>
#include <memory>

class tree;

using tree_ptr = std::shared_ptr<tree>;
static constexpr auto null_tree_ptr = nullptr;

struct tree_stats {
	int max_level;
	int n_nodes;
	int n_leaves;
	int n_part;
	int n_neighbor;
};

class tree {

	static constexpr int NMAX = 100;
	std::array<tree_ptr, NCHILD> children;
	tree_ptr parent;
	tree_ptr self;
	std::vector<particle> parts;
	bool refined;
	range box;
	int level;
	void make_tree(std::vector<particle> &&p);
	void compute_smoothing_lengths(tree_ptr root);
	void particles_in_range(std::vector<const particle*>&, const range&) const;
	void form_tree(tree_ptr, int);

public:
	void initialize(const std::function<state(vect)>&);
	void compute_interactions();
	real compute_fluxes();
	void compute_next_step(real dt);
	void find_neighbors();
	tree(std::vector<particle>&&, const range&);
	tree(std::vector<particle>&&);
	tree_stats compute_tree_statistics() const;
	void compute_smoothing_lengths();
	real compute_volumes();
	void destroy();
	void form_tree();
	void particles_in_sphere(std::vector<const particle*>&, const vect&, real) const;
	void output(const char* filename);
	std::vector<particle> gather_particles() const;

	template<class ...Args>
	static tree_ptr new_tree(Args &&...args) {
		const auto ptr = std::make_shared<tree>(std::forward<Args>(args)...);
		ptr->self = ptr;
		return ptr;
	}
};

#endif /* TREE_HPP_ */
