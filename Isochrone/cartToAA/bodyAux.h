/*
 *  bodyAux.h
 *  actionCalc
 *
 *  Created by Robyn Sanderson on 7/14/10.
 *  Auxiliary functions for taking info from body type, rotating it into xy plane, and converting to polar coords
 */
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "vector.h"

void bodyRotate(body * b, body * bp, double L[3]);
void bodyInverseRotate(body * b, body * bp, double L[3]);
void bodyStrip(body * b, double q0[2], double qdot[2]);
void bodyReplace(body * b, double x[2], double v[2]);
