#include "finiteDifferences.h"

float centeredDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0 + delta) - f(x0 - delta)) / delta;
}

float backwardDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0;
}

float forwardDifference(float x0, float delta, std::function<float(float)> f);

float centered2Difference(float x0, float delta, std::function<float(float)> f);

float padel(float x0, float delta, std::function<float(float)> f);

