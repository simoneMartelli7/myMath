#include "finiteDifferences.h"

float centeredDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0 + delta) - f(x0 - delta)) / delta;
}

float centeredDifference(std::function<float(float)> f, float xip, float xim)
{
	return (f(xip) - f(xim)) / (xip - xim);
}

float centeredDifference4(std::function<float(float)> f, float xip2, float xip, float xim, float xim2)
{
	return (-f(xip2) + 8 * f(xip) - 8 * f(xim) + f(xim2)) / (12 * (xip2 - xip));
}

float backwardDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0) - f(x0 - delta)) / delta;
}

float backwardDifference(std::function<float(float)> f, float xi, float xim)
{
	return (f(xi) - f(xim)) / (xi - xim);
}

float bacwardDifference3(std::function<float(float)> f, float xip, float xi, float xim, float xim2)
{
	return (2 * f(xip) + 3 * f(xi) - 6 * f(xim) + f(xim2)) / (6 * (xip - xi));
}

float forwardDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0 + delta) - f(x0)) / delta;
}

float forwardDifference(std::function<float(float)> f, float xi, float xip)
{
	return (f(xip) - f(xi)) / (xip - xi);
}

float forwardDifference3(std::function<float(float)> f, float xip2, float xip, float xi, float xim)
{
	return (-f(xip2) + 6 * f(xip) - 3 * f(xi) - 2 * f(xim)) / (6 * (xip - xi));
}

float centered2Difference(float x0, float delta, std::function<float(float)> f)
{
	return 0.25 * (f(x0 + delta) - 2 * f(x0) + f(x0 - delta)) / (delta * delta);
}



