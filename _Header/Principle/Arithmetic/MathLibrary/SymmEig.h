#ifndef _ARITHMETIC_SYMMETRIC_EIGENVALUE_H_
#define _ARITHMETIC_SYMMETRIC_EIGENVALUE_H_

#include <limits>
#include "MathLibrary.h"
#include "../Vector/VectorN.h"
#include "../Matrix/MatrixN.h"

namespace arithmetic {

template<typename T> void eigsrt(VectorN<T> &d, MatrixN<T> *v=NULL);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute eigenvalues and eigenvectors for real symmetric matrix a[0..n-1][0..n-1]
// Input : a		_ Input Matrix
//         yesvec	_ if true, eigenvectors are computed
//
// Output : d[0..n-1]			_ the eigenvalues of a sorted into descending order
//          z[0..n-1][0..n-1]	_ a matrix whose columns contain the corresponding normalized eigenvectors
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> class SymmEig {
public : 
	SymmEig(const MatrixN<T> &a, bool yesvec = true);
	SymmEig(const VectorN<T> &dd, const VectorN<T> &ee, bool yesvec=true);

	void sort();
	void tred2();
	void tqli();
	const T pythag(const T a, const T b);

public : 
	int n;
	MatrixN<T> z;
	VectorN<T> d,e;
	bool yesvecs;
};

template<typename T> SymmEig<T>::SymmEig(const MatrixN<T> &a, bool yesvec) : n(static_cast<int>(a.dimRow())), z(a), d(n), e(n), yesvecs(yesvec) {
	tred2();
	tqli();
	sort();
}

template<typename T> SymmEig<T>::SymmEig(const VectorN<T> &dd, const VectorN<T> &ee, bool yesvec) : n((int)dd.size()), d(dd), e(ee), z(n,n,(T)0.0f), yesvecs(yesvec) {
	z.identity();
	tqli();
	sort();
}

template<typename T> void SymmEig<T>::sort() {
	if (yesvecs)
		eigsrt(d,&z);
	else
		eigsrt(d);
}

template<typename T> void SymmEig<T>::tred2()	{
	int l,k,j,i;
	T scale,hh,h,g,f;
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += ::abs(z(i,k));
			if (scale == 0.0)
				e[i]=z(i,l);
			else {
				for (k=0;k<i;k++) {
					z(i,k) /= scale;
					h += z(i,k)*z(i,k);
				}
				f=z(i,l);
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z(i,l)=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					if (yesvecs)
						z(j,i)=z(i,j)/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z(j,k)*z(i,k);
					for (k=j+1;k<i;k++)
						g += z(k,j)*z(i,k);
					e[j]=g/h;
					f += e[j]*z(i,j);
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z(i,j);
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z(j,k) -= (f*e[k]+g*z(i,k));
				}
			}
		} else
			e[i]=z(i,l);
		d[i]=h;
	}
	if (yesvecs) d[0]=0.0;
	e[0]=0.0;
	for (i=0;i<n;i++) {
		if (yesvecs) {
			if (d[i] != 0.0) {
				for (j=0;j<i;j++) {
					g=0.0;
					for (k=0;k<i;k++)
						g += z(i,k)*z(k,j);
					for (k=0;k<i;k++)
						z(k,j) -= g*z(k,i);
				}
			}
			d[i]=z(i,i);
			z(i,i)=1.0;
			for (j=0;j<i;j++) z(j,i)=z(i,j)=0.0;
		} else {
			d[i]=z(i,i);
		}
	}
}

template<typename T> void SymmEig<T>::tqli()	{

	int m,l,iter,i,k;
	T s,r,p,g,f,dd,c,b;
	const T EPS=numeric_limits<T>::epsilon();
	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=::abs(d[m])+::abs(d[m+1]);
				if (::abs(e[m]) <= EPS*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) throw("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					if (yesvecs) {
						for (k=0;k<n;k++) {
							f=z(k,i+1);
							z(k,i+1)=s*z(k,i)+c*f;
							z(k,i)=c*z(k,i)-s*f;
						}
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

template<typename T> const T SymmEig<T>::pythag(const T a, const T b) {

	T absa=::abs(a), absb=::abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQ(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQ(absa/absb))));

}

template<typename T> void eigsrt(VectorN<T> &d, MatrixN<T> *v)	{

	int k;
	int n = static_cast<int>(d.size());
	for (int i=0;i<n-1;i++) {
		T p=d[k=i];
		for (int j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			if (v != NULL)
				for (int j=0;j<n;j++) {
					p=(*v)(j,i);
					(*v)(j,i)=(*v)(j,k);
					(*v)(j,k)=p;
				}
		}
	}
}

};		// namepsace arithmetic

#endif