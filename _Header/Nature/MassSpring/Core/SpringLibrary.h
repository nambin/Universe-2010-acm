#ifndef _MASS_SPRING_SPACE_SPRING_LIBRARY_H_
#define _MASS_SPRING_SPACE_SPRING_LIBRARY_H_

#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../../../Principle/Arithmetic/Matrix/Matrix3.h"
#include "../Core/Spring.h"
#include "../../../Basis/Basis/Array.h"

namespace springspace {

using basis::Array;
using arithmetic::Vector3;
using arithmetic::Matrix3;

template<typename T> class SpringFunctor;
class SpringLibrary {
public : 
	template<typename T>
	static inline const Vector3<T> elasticForce(const Spring<T>& s, T k, const Vec3Array<T>& x) {
		const Vector3<T>& x0 = x[s.idx0];
		const Vector3<T>& x1 = x[s.idx1];
		return (k * (1 - s.restlen/length(x0, x1))) * (x1-x0);
	}
	template<typename T>
	static void DelasticForceDxi(Matrix3<T>& mat, const Spring<T>& s, T k, const Vec3Array<T>& x) {
		const Vector3<T> xji = x[s.idx0] - x[s.idx1];
		T len = arithmetic::mag(xji);

		xji.getTensorProduct(mat);
		mat *= (s.restlen / (len*len*len));

		T tmp = (len - s.restlen) / len;
		for(Matrix3<T>::index_type i=0;i<3;++i) 
			mat(i,i) += tmp;

		mat *= (-k);
	}
	template<typename T>
	static inline const Vector3<T> dampingForce(const Spring<T>& s, T d, const Vec3Array<T>& v) {
		const Vector3<T>& v0 = v[s.idx0];
		const Vector3<T>& v1 = v[s.idx1];
		return d * (v1 - v0);
	}
	template<typename T>
	static inline void DdampingForceDvi(Matrix3<T>& mat, const Spring<T>& s, T d, const Vec3Array<T>& v) {
		mat.zero();
		for(Matrix3<T>::index_type i=0;i<3;++i)
			mat(i,i) = -d;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Other Library
	template<typename T>
	static void addSpringForce(Vec3Array<T>& force, const Array<Spring<T>>& springs, 
							   const SpringFunctor<T>* func, const Vec3Array<T>& pos, const Vec3Array<T>& vel) {
		Vector3<T> sforce;
		for(size_t i=0;i<springs.size();++i) {
			const Spring<T>& s = springs[i];

			sforce = func->f(s, pos, vel);
			force[s.idx0] += sforce;
			force[s.idx1] -= sforce;
		}
	}
	template<typename T>
	static void addSpringElasticForce(Vec3Array<T>& force, const Array<Spring<T>>& springs, 
									  const SpringFunctor<T>* func, const Vec3Array<T>& pos) {
		Vector3<T> sforce;
		for(size_t i=0;i<springs.size();++i) {
			const Spring<T>& s = springs[i];

			sforce = SpringLibrary::elasticForce(s, func->k, pos);
			force[s.idx0] += sforce;
			force[s.idx1] -= sforce;
		}
	}

	template<typename T>
	static const T maxStrain(const Array<Spring<T>>& springs, const Vec3Array<T>& pos) {
		T _max = 0, len;
		Vector3<T> v2v;

		for(size_t i=0;i<springs.size();++i) {
			const Spring<real>& s = springs[i];

			v2v = pos[s.idx1] - pos[s.idx0];
			len = mag(v2v);

			_max = basis::MAX(_max, abs((len - s.restlen) / s.restlen));
		}

		return _max;
	}
};

template<typename T>
class SpringFunctor {
public :
	SpringFunctor(T _k) : k(_k) {}
	virtual const T elasticCoeff() const { return k; }
	virtual const T dampingCoeff() const { return 0; }

	virtual const Vector3<T> f(const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const = 0;
	virtual void DfDx(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const = 0;
	virtual void DfDv(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const = 0;

	T k;
};

template<typename T>
class SpringFunctor_Elastic : public SpringFunctor<T> {
public :
	SpringFunctor_Elastic(T k) : SpringFunctor(k) {}
	const Vector3<T> f(const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		return SpringLibrary::elasticForce(s, k, x);
	}

	void DfDx(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		SpringLibrary::DelasticForceDxi(mat, s, k, x);
	}

	void DfDv(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		mat.zero();
	}	
};

template<typename T>
class SpringFunctor_ElasticDamping : public SpringFunctor<T> {
public :
	SpringFunctor_ElasticDamping(T _k, T _d) : SpringFunctor(_k), d(_d) {}
	virtual const T dampingCoeff() const { return d; }

	const Vector3<T> f(const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		return SpringLibrary::elasticForce(s, k, x) + SpringLibrary::dampingForce(s, d, v);
	}

	void DfDx(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		SpringLibrary::DelasticForceDxi(mat, s, k, x);
	}

	void DfDv(Matrix3<T>& mat, const Spring<T>& s, const Vec3Array<T>& x, const Vec3Array<T>& v) const {
		SpringLibrary::DdampingForceDvi(mat, s, d, v);
	}

	T d;
};

};		// springspace

#endif
