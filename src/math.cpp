/*
 * math.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */




#include "math.hpp"


real condition_number(const std::array<vect, NDIM> &A, std::array<vect, NDIM> &Ainv) {
	static constexpr auto NMAX = std::numeric_limits<real>::max();
	if (matrix_inverse(A, Ainv)) {
		real Asum = 0.0;
		real Ainvsum = 0.0;
		for (int n = 0; n < NDIM; n++) {
			Asum += A[n].dot(A[n]);
			Ainvsum += Ainv[n].dot(Ainv[n]);
		}
		return std::sqrt(Asum * Ainvsum) / NDIM;
	} else {
		return NMAX;
	}
}


bool matrix_inverse(const std::array<vect, NDIM> &A, std::array<vect, NDIM> &Ainv) {
#if(NDIM==1)
	if( A[0][0] == 0.0 ) {
		return false;
	}
	Ainv[0][0] = 1.0 / A[0][0];
#else
#if(NDIM==2)
	real det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
	if (det == 0.0) {
		return false;
	}
	const real detinv = 1.0 / detinv;
	Ainv[0][0] = +A[1][1] * detinv;
	Ainv[0][1] = -A[0][1] * detinv;
	Ainv[1][0] = -A[1][0] * detinv;
	Ainv[1][1] = +A[0][0] * detinv;
#else
	real det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]);
	if( det == 0.0 ) {
		return false;
	}
	det /**/-= A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]);
	det /**/+= A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	const real detinv = 1.0 / det;
	Ainv[0][0] = (-A[1][2] * A[2][1] + A[1][1] * A[2][2]) * detinv;
	Ainv[0][1] = (+A[0][2] * A[2][1] - A[0][1] * A[2][2]) * detinv;
	Ainv[0][2] = (-A[0][2] * A[1][1] + A[0][1] * A[1][2]) * detinv;
	Ainv[1][0] = (+A[1][2] * A[2][0] - A[1][0] * A[2][2]) * detinv;
	Ainv[1][1] = (-A[0][2] * A[2][0] + A[0][0] * A[2][2]) * detinv;
	Ainv[1][2] = (+A[0][2] * A[1][0] - A[0][0] * A[1][2]) * detinv;
	Ainv[2][0] = (-A[1][1] * A[2][0] + A[1][0] * A[2][1]) * detinv;
	Ainv[2][1] = (+A[0][1] * A[2][0] - A[0][0] * A[2][1]) * detinv;
	Ainv[2][2] = (-A[0][1] * A[1][0] + A[0][0] * A[1][1]) * detinv;
#endif
#endif
	return true;
}

