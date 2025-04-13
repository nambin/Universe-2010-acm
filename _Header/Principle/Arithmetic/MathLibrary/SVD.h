#ifndef _ARITHMETIC_SVD_H_
#define _ARITHMETIC_SVD_H_

#include <vector>
#include <limits>
#include "../../Principle.h"
#include "MathLibrary.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template<typename T> 
class SVD {	
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	SVD(const MatrixN<T> &a);

	//////////////////////////////////////////////////////////////////////////
	// About Matrix Space
	const int rank(double thres = -1.0f) const;
	const int nullity(double thres = -1.0f) const;
	void getRange(MatrixN<T> range, double thres = -1.0f) const;
	void getNullspace(MatrixN<T> nullspace, double thres = -1.0f) const;

	//////////////////////////////////////////////////////////////////////////
	// Solvers
	void solve(VectorN<T> &x, const VectorN<T> &b, double thres = -1.0f) const;
	void solve(MatrixN<T> &x, const MatrixN<T> &b, double thres = -1.0f) const; 

	void getPseudoInverse(MatrixN<T> &inv, double thres = -1.0f) const;

	//////////////////////////////////////////////////////////////////////////
	// Condition Number
	const T inv_ConditionNumber() const;

private : 
	//////////////////////////////////////////////////////////////////////////
	// SVD Routines
	void decompose();
	void reorder();
	const T pythag(const T a, const T b);

public : 
	int m,n;
	MatrixN<T> u,v;
	VectorN<T> w;
	mutable T eps, tsh;
};

template<typename T>
SVD<T>::SVD(const MatrixN<T> &a) : m(static_cast<int>(a.nRows())), n(static_cast<int>(a.nCols())), u(a), v(n,n), w(n) {
	eps = std::numeric_limits<T>::epsilon();
	decompose();
	reorder();
	tsh = 0.5f * sqrt(m+n+1.0f) * w[0] * eps;
}

template<typename T>
const int SVD<T>::rank(double thres) const {
	int j,nr = 0;
	tsh = static_cast<T>(thres >= 0.0f ? thres : 0.5f * sqrt(m+n+1.0f) * w[0] * eps);
	for (j=0;j<n;++j) 
		if (w[j] > tsh) 
			nr++;
	return nr;
}

template<typename T>
const int SVD<T>::nullity(double thres) const {
	int j,nn = 0;
	tsh = static_cast<T>(thres >= 0.0f ? thres : 0.5f * sqrt(m+n+1.0f) * w[0] * eps);
	for (j=0;j<n;++j) 
		if (w[j] <= tsh) 
			nn++;
	return nn;
}

template<typename T>
void SVD<T>::getRange(MatrixN<T> range, double thres) const {
	int i,j,nr=0;
	range.resize(m, rank(thres));
	
	for (j=0;j<n;++j) {
		if (w[j] > tsh) {
			for (i=0;i<m;++i) 
				range(i,nr) = u(i,j);
			++nr;
		}
	}
}

template<typename T>
void SVD<T>::getNullspace(MatrixN<T> nullspace, double thres) const {
	int j,jj,nn=0;
	nullspace.resize(n,nullity(thres));

	for (j=0;j<n;++j) {
		if (w[j] <= tsh) {
			for (jj=0;jj<n;++jj) 
				nullspace(jj,nn) = v(jj,j);
			++nn;
		}
	}
}

template<typename T>
const T SVD<T>::inv_ConditionNumber() const {
	return (w[0] <= 0.0f || w[n-1] <= 0.0f) ? 0.0f : w[n-1]/w[0];
}

template<typename T>
void SVD<T>::getPseudoInverse(MatrixN<T> &inv, double thres = -1.0f) const {

	inv.resize(n, m);

	int i,j;
	tsh = static_cast<T>(thres >= 0.0f ? thres : 0.5f * sqrt(m+n+1.0f)*w[0]*eps);

	MatrixN<T> ut(n, m);
	for (j=0;j<n;++j) {
		if (w[j] > tsh) {
			T wInv = 1.0f / w[j];
			for (i=0;i<m;++i) 
				ut(j,i) = u(i,j) * wInv;
		}
	}

	inv = v * ut;
}

template<typename T>
void SVD<T>::solve(VectorN<T> &x, const VectorN<T> &b, double thresh) const {
	
	int i,j,jj;
	T s;

	if (b.size() != m || x.size() != n) 
		throwError("SVD::solve(VectorN) bad sizes");

	VectorN<T> tmp(n);
	tsh = static_cast<T>(thresh >= 0.0f ? thresh : 0.5f * sqrt(m+n+1.0f)*w[0]*eps);

	for (j=0;j<n;++j) {
		s = 0.0f;
		if (w[j] > tsh) {
			for (i=0;i<m;++i) 
				s += u(i,j) * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j=0;j<n;++j) {
		s = 0.0f;
		for (jj=0;jj<n;++jj) 
			s += v(j,jj) * tmp[jj];
		x[j] = s;
	}
}

template<typename T>
void SVD<T>::solve(MatrixN<T> &x, const MatrixN<T> &b, double thres) const {

	int i,j,m=b.nCols();

	if (b.nRows() != n || x.nRows() != n || b.nCols() != x.nCols())
		throwError("SVD::solve(MatrixN) bad sizes");

	VectorN<T> xx(n);
	for (j=0;j<m;++j) {
		for (i=0;i<n;i++) xx[i] = b(i,j);
		solve(xx,xx,thres);
		for (i=0;i<n;i++) x(i,j) = xx[i];
	}
}

template<typename T>
void SVD<T>::decompose() {
	bool flag;
	int i,its,j,jj,k,l,nm;
	T anorm,c,f,g,h,s,scale,x,y,z;
	VectorN<T> rv1(n);
	g = scale = anorm = 0.0f;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0f;
		if (i < m) {
			for (k=i;k<m;k++) scale += ::abs(u(k,i));
			if (scale != 0.0f) {
				for (k=i;k<m;k++) {
					u(k,i) /= scale;
					s += u(k,i) * u(k,i);
				}
				f=u(i,i);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u(i,i)=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0f,k=i;k<m;k++) s += u(k,i) * u(k,j);
					f=s/h;
					for (k=i;k<m;k++) u(k,j) += f*u(k,i);
				}
				for (k=i;k<m;k++) u(k,i) *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0f;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += ::abs(u(i,k));
			if (scale != 0.0f) {
				for (k=l-1;k<n;k++) {
					u(i,k) /= scale;
					s += u(i,k) * u(i,k);
				}
				f=u(i,l-1);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u(i,l-1)=f-g;
				for (k=l-1;k<n;k++) rv1[k]=u(i,k)/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0f,k=l-1;k<n;k++) s += u(j,k)*u(i,k);
					for (k=l-1;k<n;k++) u(j,k) += s*rv1[k];
				}
				for (k=l-1;k<n;k++) u(i,k) *= scale;
			}
		}
		anorm=MAX(anorm,(::abs(w[i])+::abs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0f) {
				for (j=l;j<n;j++)
					v(j,i) = (u(i,j)/u(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0f,k=l;k<n;k++) s += u(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j) = v(j,i) = 0.0f;
		}
		v(i,i)=1.0f;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) u(i,j)=0.0f;
		if (g != 0.0f) {
			g=1.0f/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += u(k,i)*u(k,j);
				f=(s/u(i,i))*g;
				for (k=i;k<m;k++) u(k,j) += f*u(k,i);
			}
			for (j=i;j<m;j++) u(j,i) *= g;
		} else for (j=i;j<m;j++) u(j,i)=0.0;
		++u(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (l == 0 || ::abs(rv1[l]) <= eps*anorm) {
					flag=false;
					break;
				}
				if (::abs(w[nm]) <= eps*anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (::abs(f) <= eps*anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0f/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u(j,nm);
						z=u(j,i);
						u(j,nm)=y*c+z*s;
						u(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0f) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == 29) throwError("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0f*h*y);
			g=pythag(f,1.0f);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0f;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0f/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=u(jj,j);
					z=u(jj,i);
					u(jj,j)=y*c+z*s;
					u(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0f;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

template<typename T>
void SVD<T>::reorder() {
	int i,j,k,s,inc=1;
	T sw;
	VectorN<T> su(m), sv(n);
	do { inc *= 3; inc++; } while (inc <= n);
	do {
		inc /= 3;
		for (i=inc;i<n;i++) {
			sw = w[i];
			for (k=0;k<m;k++) su[k] = u(k,i);
			for (k=0;k<n;k++) sv[k] = v(k,i);
			j = i;
			while (w[j-inc] < sw) {
				w[j] = w[j-inc];
				for (k=0;k<m;k++) u(k,j) = u(k,j-inc);
				for (k=0;k<n;k++) v(k,j) = v(k,j-inc);
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for (k=0;k<m;k++) u(k,j) = su[k];
			for (k=0;k<n;k++) v(k,j) = sv[k];

		}
	} while (inc > 1);
	for (k=0;k<n;k++) {
		s=0;
		for (i=0;i<m;i++) if (u(i,k) < 0.) s++;
		for (j=0;j<n;j++) if (v(j,k) < 0.) s++;
		if (s > (m+n)/2) {
			for (i=0;i<m;i++) u(i,k) = -u(i,k);
			for (j=0;j<n;j++) v(j,k) = -v(j,k);
		}
	}
}

template<typename T>
const T SVD<T>::pythag(const T a, const T b) {
	T absa=::abs(a), absb=::abs(b);
	return (absa > absb ? absa*sqrt(1.0f + SQ(absb/absa)) :
		   (absb == 0.0f ? 0.0f : absb*sqrt(1.0f + SQ(absa/absb))));
}

};		// arithmetic

#endif
