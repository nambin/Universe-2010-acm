#pragma once

#include <cmath>
#include "../../../Basis/Basis/Constants.h"

namespace arithmetic {

template<typename T>
class Sin {
public : 
	Sin(const T _a = 1, const T _b = 0, const T _k = 1) : a(_a), b(_b), k(_k) {}

	const T operator()(const T x) const	{ return k * sin(a*x + b);		}
	const T dx(const T x) const			{ return a * k * cos(a*x + b);	}

private : 
	T a,b,k;
};

template<typename T>
class Cos {
public : 
	Cos(const T _a = 1, const T _b = 0, const T _k = 1) : a(_a), b(_b), k(_k) {}

	const T operator()(const T x) const	{ return k * cos(a*x + b);			}
	const T dx(const T x) const			{ return - a * k * sin(a*x + b);	}

private : 
	T a,b,k;
};

};			// arithmetic