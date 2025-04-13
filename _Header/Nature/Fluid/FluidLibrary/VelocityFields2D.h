#pragma once

namespace fluidspace {

template<typename T>
class TranslationVelocityField2D : public VelocityField2D<T> {
public : 
	TranslationVelocityField2D(const Vector2<T>& v) : vel(v) {}

	const T getU(const Vector2<T> &pos) const { return vel[0];	}
	const T getV(const Vector2<T> &pos) const {	return vel[1];	}	

protected :
	void calc_maxVel() const { maxVel = vel; }

	Vector2<T> vel;
};

template<typename T>
class ZeroVelocityField2D : public VelocityField2D<T> {
public : 
	ZeroVelocityField2D() {}

	const T getU(const Vector2<T> &pos) const { return T(0.0);	}
	const T getV(const Vector2<T> &pos) const {	return T(0.0);	}	

protected :
	void calc_maxVel() const { maxVel = Vector2<T>(0,0); }
};

template<typename T>
class RotationVelocityField2D : public VelocityField2D<T> {
public : 
	RotationVelocityField2D(T _step = 1, T _Lx = 1, T _Ly = 1) : stepPerOneRotation(_step), Lx(_Lx), Ly(_Ly) {}

	const T getU(const Vector2<T> &pos) const {	return (T(2.0) * PI / stepPerOneRotation) * (T(0.5)*Ly - pos[1]);	}
	const T getV(const Vector2<T> &pos) const {	return (T(2.0) * PI / stepPerOneRotation) * (pos[0] - T(0.5)*Lx);	}	

	T get_num_steps() const { return stepPerOneRotation; }
	T get_Lx() const { return Lx; }

protected :
	void calc_maxVel() const {
		size_t nX, nY; nX = nY = 200;

		Vector2<T> vel, pos;
		for(size_t i=0;i<nX;++i) { for(size_t j=0;j<nY;++j) { 

			pos[0] = Lx * T(i) / T(nX-1);
			pos[1] = Ly * T(j) / T(nY-1);

			vel = getVel(pos);

			if(abs(maxVel[0]) + abs(maxVel[1]) < abs(vel[0]) + abs(vel[1]))
				maxVel = vel;
		} }
	}

	T	stepPerOneRotation;
	T	Lx, Ly;
};

template<typename T>
class ZalesakVelocityField2D : public RotationVelocityField2D<T> {
public : 
	ZalesakVelocityField2D(std::size_t _step = 1, T L = 1) : RotationVelocityField2D<T>(_step, L, L) {}
};

template<typename T>
class VortexFlowVelocityField2D : public VelocityField2D<T> {
public : 
	VortexFlowVelocityField2D(T _L) : L(_L) {}

	const T getU(const Vector2<T> &pos) const {	return L * SQ( sin(PI*pos[0]/L) ) * ( -sin(2.0f*PI*pos[1]/L) );	}	
	const T getV(const Vector2<T> &pos) const {	return L * SQ( sin(PI*pos[1]/L) ) * (  sin(2.0f*PI*pos[0]/L) );	}

protected : 
	void calc_maxVel() const {
		basis::throwError("Not Implemented _ VortexFlowVelocityField2D");
	}

	T L;
};

};			// fluidspace