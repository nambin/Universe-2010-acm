#ifndef _BASIS_BASIC_LIBRARY_H_
#define _BASIS_BASIC_LIBRARY_H_

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include "Constants.h"
#include "../Basis.h"

namespace basis {

inline void throwError(const char* _errorMsgFormat, ...) {
	va_list ap;
	va_start(ap, _errorMsgFormat);
	vfprintf_s(stderr, _errorMsgFormat, ap);
	va_end(ap);

	std::cout << std::endl;

	system("pause");
}

template <typename T> inline void release(T *&x)	{ if(x != NULL) delete(x); x = NULL;	}

template <typename T> inline const T Rad(const T x)	{ return static_cast<T>(DegToRad) * x;	}
template <typename T> inline const T Deg(const T x)	{ return static_cast<T>(RadToDeg) * x;	}

template <typename T> inline const T MIN(const T a, const T b)				{ return (a < b) ? a : b;	}
template <typename T> inline const T MAX(const T a, const T b)				{ return (a > b) ? a : b;	}
template <typename T> inline const T MIN3(const T a, const T b, const T c)	{ return MIN(a, MIN(b,c));	}
template <typename T> inline const T MAX3(const T a, const T b, const T c)	{ return MAX(a, MAX(b,c));	}

template <typename T> inline const T ABS(const T a)		{ return ::abs(a);	}
template <typename T> inline const T SQ(const T a)		{ return a*a;	}
template <typename T> inline const T SQUARE(const T a)	{ return a*a;	}
template <typename T> inline const T POW2(const T a)	{ return a*a;	}
template <typename T> inline const T CB(const T a)		{ return a*a*a;	}
template <typename T> inline const T CUBE(const T a)	{ return a*a*a;	}
template <typename T> inline const T POW3(const T a)	{ return a*a*a;	}

template <typename T> 
class TypeChecker {
public : 
	TypeChecker() {}
	template<typename S> const bool operator()(const S* const ptr) const {
		const T* const ret = dynamic_cast<const T* const>(ptr);

		if(ret) return true;
		else	return false;
	}
};

template <typename T> 
inline const bool IsAlmostSame( const T x, const T y, const T eps = static_cast<T>(EPSILON) ) {
	if( abs( x - y ) <= eps ) return true;
	else return false;
}

template <typename T> 
inline const bool IsAlmostZero( const T x, const T eps = static_cast<T>(EPSILON) ) {
	if( abs( x ) <= eps ) return true;
	else return false;
}

template <typename T> 
inline const int SIGN( const T x ) {
	if( x > 0 )	  return  1;
	else		  return -1;
}

template <typename T> 
inline const int Sign( const T x, const T eps = static_cast<T>(EPSILON) ) {
	if( IsAlmostZero(x,eps) )	return  0;
	else if( x > 0 )			return  1;
	else						return -1;
}

template <typename T>
inline const T clamp(T val, T min, T max) {
	if(val < min)		return min;
	else if(val > max)	return max;
	return val;
}

BASISDLL void byteSwap( unsigned char* b, int n );
#define BYTE_SWAP(x) byteSwap((unsigned char *) &x,sizeof(x))

};		// namespace basis

#endif