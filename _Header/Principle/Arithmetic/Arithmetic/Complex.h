#ifndef _ARITHMETIC_COMPLEX_H_
#define _ARITHMETIC_COMPLEX_H_

namespace arithmetic {

template<typename T>
class Complex {
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	Complex();
	Complex(T _x, T _y = 0);

	//////////////////////////////////////////////////////////////////////////
	// Getter
	T	real() const		{ return x; }
	T	imaginary() const	{ return y; }
	T&	real()				{ return x; }
	T&	imaginary()			{ return y; }

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	Complex& zero()			{ x = y = 0;		}	
	void set(T _x, T _y)	{ x = _x; y = _y;	}	

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline Complex&	normalize();
	inline const T	mag() const;
	inline const T	sqrmag() const;

	inline Complex	getInverse() const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Complex& operator+=(T d);	
	inline Complex& operator-=(T d);	
	inline Complex& operator*=(T d);	
	inline Complex& operator/=(T d);	

	inline Complex& operator+=(const Complex&); 
	inline Complex& operator-=(const Complex&); 
	inline Complex& operator*=(const Complex&); 
	inline Complex& operator/=(const Complex&); 

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<( std::ostream&, const Complex<T>& );
	template <typename T> friend std::istream& operator>>( std::istream&, Complex<T>& );

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend const Complex<T> operator+(const Complex<T> &a);
	template <typename T> friend const Complex<T> operator-(const Complex<T> &a);
	template <typename T> friend const Complex<T> operator+(const Complex<T> &a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator-(const Complex<T> &a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator*(const Complex<T> &a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator/(const Complex<T> &a,const Complex<T> &b);

	template <typename T> friend const Complex<T> operator+(const Complex<T> &a,T b);
	template <typename T> friend const Complex<T> operator-(const Complex<T> &a,T b);
	template <typename T> friend const Complex<T> operator*(const Complex<T> &a,T b);
	template <typename T> friend const Complex<T> operator/(const Complex<T> &a,T b);

	template <typename T> friend const Complex<T> operator+(T a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator-(T a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator*(T a,const Complex<T> &b);
	template <typename T> friend const Complex<T> operator/(T a,const Complex<T> &b);

protected : 
	T	x, y;
};

template<typename T> Complex<T>::Complex() : x(0), y(0) {}
template<typename T> Complex<T>::Complex(T _x, T _y) : x(_x), y(_y) {}	

template<typename T> inline Complex<T>& Complex<T>::normalize()		{ return (*this) /= mag();			}
template<typename T> inline const T Complex<T>::mag() const			{ return sqrt(x*x + y*y);			}
template<typename T> inline const T Complex<T>::sqrmag() const			{ return (x*x + y*y);				}
template<typename T> inline Complex<T> Complex<T>::getInverse() const	{ return Complex<T>(x,-y) / mag();	}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template<typename T> inline Complex<T>& Complex<T>::operator+=(T d) { x += d;			}
template<typename T> inline Complex<T>& Complex<T>::operator-=(T d) { x -= d;			}
template<typename T> inline Complex<T>& Complex<T>::operator*=(T d) { x *= d; y *= d;	}	
template<typename T> inline Complex<T>& Complex<T>::operator/=(T d) { x /= d; y /= d;	}	

template<typename T> inline Complex<T>& Complex<T>::operator+=(const Complex& c) { x += c.x; y += c.y; return (*this); }
template<typename T> inline Complex<T>& Complex<T>::operator-=(const Complex& c) { x -= c.x; y -= c.y; return (*this); }
template<typename T> inline Complex<T>& Complex<T>::operator*=(const Complex& c) { return (*this) = (*this) * c; }
template<typename T> inline Complex<T>& Complex<T>::operator/=(const Complex& c) { return (*this) = (*this) / c; }

//////////////////////////////////////////////////////////////////////////
// IO
template<typename T> std::ostream& operator<<( std::ostream& os, const Complex<T>& c ) {
	os << "( " << c.x << " + " << c.y << " i )";
	return os;
}

template<typename T> std::istream& operator>>( std::istream& is, Complex<T>& c ) {
	static char	buf[64];	
	is >> buf >> c.x >> buf >> c.y >> buf;
	return is;
}

//////////////////////////////////////////////////////////////////////////
// Friends Operators
template <typename T> inline const Complex<T> operator+(const Complex<T> &a) {	return a;						}
template <typename T> inline const Complex<T> operator-(const Complex<T> &a) {	return Complex<T>(-a.x, -a.y);	}

template <typename T> inline const Complex<T> operator+(const Complex<T>& a,const Complex<T>& b) { return Complex<T>(a.x + b.x, a.y + b.y); }
template <typename T> inline const Complex<T> operator-(const Complex<T>& a,const Complex<T>& b) { return Complex<T>(a.x - b.x, a.y - b.y); }
template <typename T> inline const Complex<T> operator*(const Complex<T>& a,const Complex<T>& b) { return Complex<T>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
template <typename T> inline const Complex<T> operator/(const Complex<T>& a,const Complex<T>& b) { return a * b.getInverse(); }

template <typename T> inline const Complex<T> operator+(const Complex<T>& a,T b) { return Complex<T>(a.x + b, a.y); }
template <typename T> inline const Complex<T> operator-(const Complex<T>& a,T b) { return Complex<T>(a.x - b, a.y); }
template <typename T> inline const Complex<T> operator*(const Complex<T>& a,T b) { return Complex<T>(a.x * b, a.y * b); }
template <typename T> inline const Complex<T> operator/(const Complex<T>& a,T b) { return Complex<T>(a.x / b, a.y / b); }

template <typename T> inline const Complex<T> operator+(T a,const Complex<T> &b) { return b + a; }
template <typename T> inline const Complex<T> operator-(T a,const Complex<T> &b) { return Complex<T>(a-b.x, -b.y); }
template <typename T> inline const Complex<T> operator*(T a,const Complex<T> &b) { return b * a; }
template <typename T> inline const Complex<T> operator/(T a,const Complex<T> &b) { return Complex<T>(a) / b;		}

};		// arithmetic

#endif