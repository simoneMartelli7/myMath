#pragma once

#include <functional>

class ODE
{
private:
	float order;
	std::function<float(float)> f;
	float cauchy;
	float neumann;

};

