/*
 *  bodyAux.c
 *  cloudMap
 *
 *  Created by Robyn Sanderson on 3/23/12.
 *  Auxiliary functions for taking info from body type, rotating it into xy plane, and converting to polar coords
 *
 */
#include "bodyAux.h"

void bodyRotate(body * b, body * bp, double L[3])
{
	
	//L is the vector that should be rotated to point along the Z axis
	
	rotateIntoZ(b->q, bp->q, L);
	rotateIntoZ(b->p, bp->p, L);
	
	//check for consv of L
	
	
	double oldL = norm(L);
	cross(bp->q, bp->p, L);
	double newL = norm(L);
	
       	assert(fabs((oldL-newL)/oldL)<1e-12);
}

void bodyInverseRotate(body * b, body * bp, double L[3])
{
	
	//Undoes bodyRotate
	
	rotateFromZ(b->q, bp->q, L);
	rotateFromZ(b->p, bp->p, L);
	
	//check for consv of L
	
	double oldL = norm(L);
	cross(bp->q, bp->p, L);
	double newL = norm(L);
	
	assert(fabs(oldL-newL)<1e-12);
}

void bodyStrip(body * b, double q0[2], double qdot[2]) {
	
	body * bp;
	
	bp = malloc(sizeof(body));
	
	double L[3];
	
	cross(b->q, b->p, L);
	
	bodyRotate(b, bp, L);
	
	q0[0] = norm(bp->q);
	q0[1] = phi(bp->q);
	
	double v[3];
	int i;
	
	for(i=0;i<3;i++) {
		v[i] = bp->p[i]/b->m;	
	}
	
	qdot[0] = vr(bp->q, v);
	qdot[1] = vphi(bp->q,v)/q0[0];
	
	free(bp);
	
}

void bodyReplace(body * b, double x[2], double v[2]) {
	
	body * bp;
	
	bp = malloc(sizeof(body));
	
	bp->q[0] = x[0] * cos(x[1]);
	bp->q[1] = x[0] * sin(x[1]);
	bp->q[2] = 0.0;
	
	bp->p[0] = b->m * (v[0] * cos(x[1]) - x[0] * v[1] * sin(x[1]));
	bp->p[1] = b->m * (v[0] * sin(x[1]) + x[0] * v[1] * cos(x[1]));
	bp->p[2] = 0.0;
	
	double L[3];
	
	cross(b->q, b->p, L);
	
	bodyInverseRotate(bp, b, L);
	
	free(bp);
}
