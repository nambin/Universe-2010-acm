#pragma once

#include "VelocityFields2D.h"

namespace fluidspace {

template<typename T>
class TranslationVelocityField3D : public VelocityField3D<T> {
public : 
	TranslationVelocityField3D(const Vector3<T>& v) : vel(v) {}

	TranslationVelocityField3D<T>* clone() const {
		return new TranslationVelocityField3D<T>(vel);
	}

	const T getU(const Vector3<T> &pos) const { return vel[0];	}
	const T getV(const Vector3<T> &pos) const {	return vel[1];	}	
	const T getW(const Vector3<T> &pos) const {	return vel[2];	}	

protected :
	void calc_maxVel() const { maxVel = vel; }

	Vector3<T> vel;
};

template<typename T>
class RotationVelocityField3D : public VelocityField3D<T> {
public : 
	RotationVelocityField3D(T _step, T _L, const Vector3<T>& _pivot_center, const Vector3<T>& _pivot_direction)
	: stepPerOneRotation(_step), L(_L), pivot_center(_pivot_center), pivot_direction(unitVector(_pivot_direction)) {}

	RotationVelocityField3D<T>* clone() const {
		return new RotationVelocityField3D(stepPerOneRotation, L, pivot_center, pivot_direction);
	}

	const Vector3<T> getVel(const Vector3<T> &pos) const { 
		const Vector3<T> diff = pos - pivot_center;
		const Vector3<T> cross_vector = cross(pivot_direction, diff);
		const T len_cross_vector = arithmetic::mag(cross_vector);

		if(len_cross_vector < 1.0e-5)
			return Vector3<T>(0);

		const Vector3<T> vel_dir = cross_vector / len_cross_vector;		

		const T radius = arithmetic::mag(diff - pivot_direction * dot(diff, pivot_direction));

		return vel_dir * T(2 * PI * radius / T(stepPerOneRotation));	
	}

	const T getU(const Vector3<T> &pos) const { return getVel(pos)[0];	}
	const T getV(const Vector3<T> &pos) const {	return getVel(pos)[1];	}	
	const T getW(const Vector3<T> &pos) const {	return getVel(pos)[2];	}	

protected : 
	void calc_maxVel() const {
		size_t nX, nY, nZ; nX = nY = nZ = 200;

		Vector3<T> vel, pos;
		for(size_t i=0;i<nX;++i) { for(size_t j=0;j<nY;++j) { for(size_t k=0;k<nZ;++k) {

			pos[0] = L * T(i) / T(nX-1);
			pos[1] = L * T(j) / T(nY-1);
			pos[2] = L * T(k) / T(nZ-1);

			vel = getVel(pos);

			if(abs(maxVel[0]) + abs(maxVel[1]) + abs(maxVel[2]) < abs(vel[0]) + abs(vel[1]) + abs(vel[2]))
				maxVel = vel;
		} } }

		cout << "Maximum Velocity in RotationVelocityField3D : " << maxVel << endl;
	}

protected : 
	T			stepPerOneRotation, L;
	Vector3<T>	pivot_center;
	Vector3<T>	pivot_direction;
};

template<typename T>
class ZalesakVelocityField3D : public VelocityField3D<T> {
public : 
	ZalesakVelocityField3D(std::size_t _step = 1, T L = 1) : zalesak(_step, L) {}

	ZalesakVelocityField3D<T>* clone() const { return new ZalesakVelocityField3D<T>(zalesak.get_num_steps(), zalesak.get_Lx()); }

	const T getU(const Vector3<T> &pos) const { return zalesak.getU(Vector2<T>(pos[0], pos[1]));	}
	const T getV(const Vector3<T> &pos) const {	return zalesak.getV(Vector2<T>(pos[0], pos[1]));	}	
	const T getW(const Vector3<T> &pos) const {	return 0;											}	

protected : 
	void calc_maxVel() const { 
		maxVel = Vector3<T>(zalesak.get_maxVel()[0], zalesak.get_maxVel()[1], 0);
	}

	ZalesakVelocityField2D<T> zalesak;
};

template<typename T>
class LeVequeVelocityField3D : public VelocityField3D<T> {
public : 
	LeVequeVelocityField3D(T _L = 1) : L(_L) {}
	const T getU(const Vector3<T> &pos) const {	return 2.0f * SQ(sin(PI*pos[0]/L)) * sin(2.0f*PI*pos[1]/L) * sin(2.0f*PI*pos[2]/L) * L;	}
	const T getV(const Vector3<T> &pos) const {	return - sin(2.0f*PI*pos[0]/L) * SQ(sin(PI*pos[1]/L)) * sin(2.0f*PI*pos[2]/L) * L;		}
	const T getW(const Vector3<T> &pos) const {	return - sin(2.0f*PI*pos[0]/L) * sin(2.0f*PI*pos[1]/L) * SQ(sin(PI*pos[2]/L)) * L;		}

	//ret[0] = 2*sin(M_PI*x[0])*sin(M_PI*x[0])*sin(2*M_PI*x[1])*sin(2*M_PI*x[2]);
	//ret[1] = -sin(2*M_PI*x[0])*sin(M_PI*x[1])*sin(M_PI*x[1])*sin(2*M_PI*x[2]);
	//ret[2] = -sin(2*M_PI*x[0])*sin(2*M_PI*x[1])*sin(M_PI*x[2])*sin(M_PI*x[2]);

protected : 
	void calc_maxVel() const {
		size_t nX, nY, nZ; nX = nY = nZ = 200;

		Vector3<T> vel, pos;
		for(size_t i=0;i<nX;++i) { for(size_t j=0;j<nY;++j) { for(size_t k=0;k<nZ;++k) {

			pos[0] = L * T(i) / T(nX-1);
			pos[1] = L * T(j) / T(nY-1);
			pos[2] = L * T(k) / T(nZ-1);

			vel = getVel(pos);

			if(abs(maxVel[0]) + abs(maxVel[1]) + abs(maxVel[2]) < abs(vel[0]) + abs(vel[1]) + abs(vel[2]))
				maxVel = vel;
		} } }
	}

	T	L;
};

template<typename T>
class EnrightVelocityField3D : public LeVequeVelocityField3D<T> {
public : 
	EnrightVelocityField3D(T _L = 1) : LeVequeVelocityField3D<T>(_L) {}

	EnrightVelocityField3D<T>* clone() const { return new EnrightVelocityField3D<T>(L); }
	
	//const T getU(const Vector3<T> &pos) const {	return cos(PI*time/3.0) * LeVequeVelocityField3D<real>::getU(pos);	}
	//const T getV(const Vector3<T> &pos) const {	return cos(PI*time/3.0) * LeVequeVelocityField3D<real>::getV(pos);	}
	//const T getW(const Vector3<T> &pos) const {	return cos(PI*time/3.0) * LeVequeVelocityField3D<real>::getW(pos);	}

	const T getU(const Vector3<T> &pos) const {	if(time < 1.5) return LeVequeVelocityField3D<real>::getU(pos); else return - LeVequeVelocityField3D<real>::getU(pos); }
	const T getV(const Vector3<T> &pos) const {	if(time < 1.5) return LeVequeVelocityField3D<real>::getV(pos); else return - LeVequeVelocityField3D<real>::getV(pos); }
	const T getW(const Vector3<T> &pos) const {	if(time < 1.5) return LeVequeVelocityField3D<real>::getW(pos); else return - LeVequeVelocityField3D<real>::getW(pos); }
};

};			// fluidspace
