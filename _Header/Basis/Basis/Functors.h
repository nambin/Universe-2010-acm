#ifndef _BASIS_FUNCTORS_H_
#define _BASIS_FUNCTORS_H_

namespace basis {

template<typename T>
class Identity {
public :
	inline T operator()(T x) const { return x; }
};

template<typename T>
class Negative {
public :
	inline T operator()(T x) const { return -x; }
};

template<typename T>
class Multiply {
public :
	Multiply(const T _mul) : mul(_mul) {}
	inline T operator()(T x) const { return mul * x; }

	T mul;
};

template<typename T>
class Absolute {
public :
	inline T operator()(T x) const { return abs(x); }
};

template<typename T>
class NegativeAbsolute {
public :
	inline T operator()(T x) const { return -abs(x); }
};

};		// basis

#endif