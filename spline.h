#pragma once
#include "polynomial.h"
#include "Vector.h"
#include "Matrix.h"




class spline
{
public:
	Vector nodes;
	polynomial* polynomials;
	

public:
	spline() {
		this->nodes = Vector();
		this->polynomials = {};
	}
	//defaults to natural
	spline(Vector& nodesX, Vector& nodesY);

	// the flag determines the type of spline 
	//	1: natural spline 
	//	2: periodic spline 
	//	3: not-a-knot
	spline(Vector& nodesX, Vector& nodesY, int flag);
	
	// this constructor is only for clamped splines for they require additional arguments 
	spline(Vector& nodesX, Vector& nodesY, float y10, float y1n);


	// identifies the correct interval for the point x0
	// this is a first implementation for the sake of testing the results, 
	// need to implement a more elegant solution 
	int findInterval(float x0);
	//evaluates the spline in x0
	float eval(float x0);

};

