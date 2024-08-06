#pragma once
#include "Vector.h"
#include "functionalMatrix.h"
#include "functionalVector.h"

// x0 is the inital guess, F contains the expressions for the functions, h determines the spacing to compute derivatives
// for the time being im choosing an uniform spacing for all derivatives altough I'm planning to change it with an adaptive method 
Vector newton(Vector x0, functionalVector F, float h, int maxIter, float tol);

// m is the order or the method, or the number of inner steps for each iteration
// for a more clear explanation see ISBN 978-3-0365-9215-2 ch. 7
Vector newtonMultiStep(int m, Vector x0, functionalVector F, float h, int maxIter, float tol);
