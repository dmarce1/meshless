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
#include <memory>

class tree;

using tree_ptr = std::shared_ptr<tree>;
static constexpr auto null_tree_ptr = nullptr;

class tree {

	static constexpr int NMAX = 1000;
	std::array<tree_ptr, NCHILD> children;
	std::array<tree_ptr, NSIBLING> siblings;
	tree_ptr parent;
	tree_ptr self;
	std::vector<particle> parts;
	bool refined;
	range box;

	void make_tree(std::vector<particle> &&p);
	void find_siblings(std::array<tree_ptr, NSIBLING>&&);
	void compute_smoothing_lengths(tree_ptr root);

public:
	tree(std::vector<particle>&&);
	tree(std::vector<particle>&&, const range&, tree_ptr parent);
	void compute_smoothing_lengths();
	void destroy();
	void find_siblings();
	std::vector<const particle*> particles_in_range(const range&) const;

	template<class ...Args>
	static tree_ptr new_tree(Args &&...args) {
		const auto ptr = std::make_shared<tree>(std::forward<Args>(args)...);
		ptr->self = ptr;
		return ptr;
	}

};
#endif /* TREE_HPP_ */
