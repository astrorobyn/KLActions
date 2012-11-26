/* vector.h: Inlined vector functions.

Copyright (C) 2007 Will M. Farr <farr@mit.edu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include<math.h>
#include<stddef.h>

static inline double dot(double x[3], double y[3]) {
	size_t i;
	double sum = 0.0;
	
	for(i = 0; i < 3; i++) {
		sum += x[i]*y[i];
	}
	
	return sum;
}

static inline double norm(double x[3]) {
	return sqrt(dot(x,x));
}

static inline void cross(double x[3], double y[3], double cp[3]){

	cp[0] = x[1] * y[2] - x[2] * y[1];
	cp[1] = x[2] * y[0] - x[0] * y[2];
	cp[2] = x[0] * y[1] - x[1] * y[0];
	
	
}

static inline double distance_squared(double x[3], double y[3]) {
	size_t i;
	double sum = 0.0;
	
	for (i = 0; i < 3; i++) {
		double temp = x[i] - y[i];
		sum += temp*temp;
	}
	
	return sum;
}

static inline double distance(double x[3], double y[3]) {
	return sqrt(distance_squared(x,y));
}

static inline double theta(double x[3]) {
	return atan2(sqrt(x[0] * x[0] + x[1] * x[1]), x[2]);
}

static inline double phi(double x[3]) {

	return atan2(x[1], x[0]);
	
}

static inline double vr(double x[3], double v[3]) {

	double th = theta(x);
	double ph = phi(x);
	
	return v[0] * sin(th) * cos(ph) + v[1] * sin(th) * sin(ph) + v[2] * cos(th);
	
}

static inline double vphi(double x[3], double v[3]) {

	double ph = phi(x);

	return -1.0 * v[0] * sin(ph) + v[1] * cos(ph);
}


static inline void rotateIntoZ(double v[3], double vp[3], double L[3]) {

	//L is the vector that should be rotated to point along the Z axis

	double R[3][3];

	double Lx = L[0];
	double Ly = L[1];
	double Lz = L[2];
	
	double lsq = Lx * Lx + Ly * Ly;
	double lmag = norm(L);
	
	int i,j;
	
	//compute rotation matrix
	
	R[0][0] = 1.0 + Lx * Lx / lsq * (Lz / lmag - 1.0);
	R[1][1] = 1.0 + Ly * Ly / lsq * (Lz / lmag - 1.0);

	R[0][1] = Lx * Ly / lsq * (Lz / lmag - 1.0);
	R[1][0] = R[0][1];
	
	R[2][0] = Lx / lmag;
	R[2][1] = Ly / lmag;
	
	R[0][2] = -R[2][0];
	R[1][2] = -R[2][1];
	
	R[2][2] = Lz / lmag;
	
	// initialize rotated vector
	
	for (i=0; i<3; i++) {
		vp[i] = 0;
	}
	
	// do dot product
	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			vp[i] += R[i][j] * v[j];
		}
	}
}

static inline void rotateFromZ(double v[3], double vp[3], double L[3]) {
	
	// undoes rotateIntoZ
	
	double R[3][3];
	
	double Lx = L[0];
	double Ly = L[1];
	double Lz = L[2];
	
	double lsq = Lx * Lx + Ly * Ly;
	double lmag = norm(L);
	
	int i,j;
	
	//compute rotation matrix
	
	R[0][0] = 1.0 + Lx * Lx / lsq * (Lz / lmag - 1.0);
	R[1][1] = 1.0 + Ly * Ly / lsq * (Lz / lmag - 1.0);
	
	R[0][1] = Lx * Ly / lsq * (Lz / lmag - 1.0);
	R[1][0] = R[0][1];
	
	R[2][0] = -Lx / lmag;
	R[2][1] = -Ly / lmag;
	
	R[0][2] = -R[2][0];
	R[1][2] = -R[2][1];
	
	R[2][2] = Lz / lmag;
	
	// initialize rotated vector
	
	for (i=0; i<3; i++) {
		vp[i] = 0;
	}
	
	// do dot product
	
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			vp[i] += R[i][j] * v[j];
		}
	}
}
#endif /* __VECTOR_H__ */
