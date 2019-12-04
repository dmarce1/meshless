/*
 * silo.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include "delaunay.hpp"
#include <silo.h>
#include <string>

void output_silo(const std::vector<particle> &parts, const char *filename) {
	DBfile *db = DBCreateReal(filename, DB_CLOBBER, DB_LOCAL, "Meshless",
	DB_HDF5);
#if( NDIM == 1)
	const int shapetypes[1] = {[DB_ZONETYPE_BEAM};
#elif( NDIM ==2)
	const int shapetypes[1] = { DB_ZONETYPE_TRIANGLE };
#else
	const int shapetypes[1] = {DB_ZONETYPE_TET};
#endif
	const int shapesize[1] = { NDIM + 1 };
	const auto delaunay_regions = compute_delaunay_regions(parts);
	const int nzones = delaunay_regions.size();
	const int nnodes = parts.size();
	const int shapecount[1] = { nzones };
	std::vector<int> node_list;
	node_list.reserve(NDIM * delaunay_regions.size());
	for (const auto r : delaunay_regions) {
		for (int n = 0; n < NDIM + 1; n++) {
			node_list.push_back(r[n]);
		}
	}
	real* coords[NDIM];
	char* coordnames[NDIM];
	for (int dim = 0; dim < NDIM; dim++) {
		coords[dim] = new real[nnodes];
		coordnames[dim] = new char[2];
		coordnames[dim][0] = 'x' + dim;
		coordnames[dim][1] = '\0';
		for (int i = 0; i < nnodes; i++) {
			coords[dim][i] = parts[i].x[dim];
		}
	}
	DBPutZonelist2(db, "zonelist", nzones, NDIM, node_list.data(), node_list.size(), 0, 0, 0, shapetypes, shapesize, shapecount, 1,
	NULL);
	DBPutUcdmesh(db, "mesh", NDIM, coordnames, coords, nnodes, nzones, "zonelist", NULL, DB_DOUBLE, NULL);
	std::vector<real> rho;
	std::vector<real> energy;
	std::array<std::vector<real>, NDIM> mom;
	rho.reserve(nnodes);
	energy.reserve(nnodes);
	for (int dim = 0; dim < NDIM; dim++) {
		mom[dim].reserve(nnodes);
	}
	for (const auto& pi : parts) {
		rho.push_back(pi.st.mass() / pi.V);
		energy.push_back(pi.st.energy() / pi.V);
		for (int dim = 0; dim < NDIM; dim++) {
			mom[dim].push_back(pi.st.momentum(dim)/ pi.V) ;
		}
	}
	DBPutUcdvar1(db, "density", "mesh", rho.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	DBPutUcdvar1(db, "energy", "mesh", energy.data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	for (int dim = 0; dim < NDIM; dim++) {
		std::string mom_name = std::string() + char('x' + char(dim)) + "_momentum";
		DBPutUcdvar1(db, mom_name.c_str(), "mesh", mom[dim].data(), nnodes, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
	}
	for (int dim = 0; dim < NDIM; dim++) {
		delete[] coordnames[dim];
		delete[] coords[dim];
	}
	DBClose(db);
}
