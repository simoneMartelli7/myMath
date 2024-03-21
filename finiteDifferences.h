#pragma once

#include <functional>

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

