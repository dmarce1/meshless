/*
 * delauny.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include "delaunay.hpp"
#include <cassert>

std::vector<delauny_region> compute_delaunay_regions(const std::vector<particle> &parts) {
	std::vector<delauny_region> regions;
	const int sz = parts.size();
#if(NDIM==1)
	regions.resize(sz - 1);
	std::vector<const particle*> ptrs(parts.size());
	for (int i = 0; i < sz; i++) {
		ptrs[i] = &parts[i];
	}
	std::sort(ptrs.begin(), ptrs.end(), [](const particle *a, const particle *b) {
		return a->x[0] < b->x[0];
	});
	for (int i = 0; i < sz - 1; i++) {
		delauny_region &r = regions[i];
		r[0] = ptrs[i]->x[0];
		r[1] = ptrs[i + 1]->x[0];
	}
#else
	FILE *fp = fopen("points.dat", "wt");
	fprintf(fp, "%i D%i\n", NDIM, NDIM);
	fprintf(fp, "%i\n", sz);
	for (const auto &pi : parts) {
		for (int dim = 0; dim < NDIM; dim++) {
			fprintf(fp, "%e ", pi.x[dim]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	if (system("qdelaunay i Qt < points.dat > delaunay.dat") != 0) {
		printf("qhull must be installed\n");
		abort();
	}
	int dsz;
	constexpr auto buffer_size = 1024;
	static char buffer[buffer_size + 1];
	fp = fopen("delaunay.dat", "rt");
	if (fgets(buffer, buffer_size, fp) == 0) {
		assert(false);
		goto SYSTEM_ERROR;
	}
	dsz = std::atoi(buffer);
	regions.reserve(dsz);
	for (int i = 0; i < dsz; i++) {
		if (fgets(buffer, buffer_size, fp) == 0) {
			assert(false);
			goto SYSTEM_ERROR;
		}
		delauny_region r;
		char *ptr = buffer;
		for (int n = 0; n < NDIM + 1; n++) {
			r[n] = std::atoi(ptr);
			while (std::isspace(*ptr)) {
				ptr++;
			}
			while (!std::isspace(*ptr)) {
				ptr++;
			}
		}
		regions.push_back(r);
	}
	fclose(fp);
	return regions;
#endif
	SYSTEM_ERROR: printf("Problem with Delauay\n");
	abort();
}
