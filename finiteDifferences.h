#pragma once

#include <functional>

float centeredDifference(float x0, float delta, std::function<float(float)> f);

float backwardDifference(float x0, float delta, std::function<float(float)> f);

float forwardDifference(float x0, float delta, std::function<float(float)> f);

float centered2Difference(float x0, float delta, std::function<float(float)> f);

float padel(float x0, float delta, std::function<float(float)> f);
