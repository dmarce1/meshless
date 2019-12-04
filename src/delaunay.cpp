/*
 * delauny.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include "delaunay.hpp"

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
	system("qdelaunay i Qt < points.dat > delaunay.dat");
	system("rm points.dat");
	fp = fopen("delaunay.dat", "rt");
	constexpr auto buffer_size = 1024;
	static char buffer[buffer_size + 1];
	fgets(buffer, buffer_size, fp);
	const int dsz = std::atoi(buffer);
	regions.resize(dsz);
	for( int i = 0; i < dsz; i++) {
		fgets(buffer, buffer_size, fp);
		delauny_region r;
		char *ptr = buffer;
		for (int n = 0; n < NDIM + 1; n++) {
			r[n] = std::atoi(ptr);
			while (!std::isspace(*ptr)) {
				ptr++;
			}
		}
		regions.push_back(r);
	}
	fclose(fp);
#endif
	return regions;
}
