#pragma once
#include "Vector.h"
#include <functional>
#include "Matrix.h"
#include "functionalVector.h"

float centeredDifference(float x0, float delta, std::function<float(float)> f);

float centeredDifference(std::function<float(float)> f, float xip, float xim);

//only for equally spaced grids 
float centeredDifference4(std::function<float(float)> f, float xip2, float xip, float xim, float xim2);

float backwardDifference(float x0, float delta, std::function<float(float)> f);

float backwardDifference(std::function<float(float)> f, float xi, float xim);

//only for equally spaced grids 
float bacwardDifference3(std::function<float(float)> f, float xip, float xi, float xim, float xim2);

float forwardDifference(float x0, float delta, std::function<float(float)> f);

float forwardDifference(std::function<float(float)> f, float xip, float xi);

//only for equally spaced grids 
float forwardDifference3(std::function<float(float)> f, float xip2, float xip, float xi, float xim);



// SECOND DERIVATIVE
float centered2Difference(float x0, float delta, std::function<float(float)> f);




// MULTIVARIABLE DIFFERENCES 

// x0 contains the values of the indipendent vaiables where we want to evaluate the derivative, 
// i is the index of the variable with respect to which we're differentiating, h is the spacing of the discretization 
float centeredDifferenceMulti(Vector x0, int i, float h, std::function<float(float*)> f);

// evaluates the jacobian of the field F in the point x0, h is the, uniform, spacing of the discreatization
Matrix Jacobian(functionalVector F, float h, Vector x0);