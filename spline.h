#pragma once
#include "polynomial.h"
#include "Vector.h"
#include "Matrix.h"


class spline
{
public:
	int nNodes;
	float* nodes; // this can be better
	polynomial* polynomials;
public:
	spline() {
		this->nNodes = 3;
		this->nodes = new float[nNodes] {0.0};
		this->polynomials = new polynomial[1];
	}
	//defaults to natural
	spline(Vector& nodesX, Vector& nodesY);

	// the flag determines the type of spline 
	//	1: natural spline 
	//	2: periodic spline 
	//	3: not-a-knot
	//spline(Vector& nodesX, Vector& nodesY, int flag);
	
};

