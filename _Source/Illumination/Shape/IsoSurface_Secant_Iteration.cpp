#include <stdafx.h>

using namespace std;
using namespace illuminationspace;

inline const real trilerp(const Vector3<real> w, const real v[8]) {
	return (1-w[0]) * ((1-w[1]) * ((1-w[2]) * v[0] + w[2] * v[1])   
						+ w[1]  * ((1-w[2]) * v[2] + w[2] * v[3]))
			+ w[0]  * ((1-w[1]) * ((1-w[2]) * v[4] + w[2] * v[5]) 
						+ w[1]  * ((1-w[2]) * v[6] + w[2] * v[7]));
}

const real IsoSurface::IsoSurfaceValue::secant_iteration(const real v[8], const Vector3<real>& w0, const Vector3<real>& w_icr, const real tol, const bool need_exact_time) const {
	real v_[4];

	if	   (abs(v_[0] = trilerp(w0                               , v)) < tol)	return 0;						
	else if(abs(v_[1] = trilerp(w0 + 0.33333333333333333f * w_icr, v)) < tol)	return 0.33333333333333333f;		
	else if(abs(v_[2] = trilerp(w0 + 0.66666666666666667f * w_icr, v)) < tol)	return 0.66666666666666667f;	
	else if(abs(v_[3] = trilerp(w0 +						w_icr, v)) < tol)	return 1.00000000000000000f;	
	else																		return secant_iteration(v_, tol, need_exact_time);	
}

const real IsoSurface::IsoSurfaceValue::secant_iteration(const real v[4], const real tol, const bool need_exact_time) const {

	real discrim;
	real T0, T1, V0, V1;
	real TM, VM, alpha;	

	const real A = -4.5f * v[0] + 13.5f * v[1] - 13.5f * v[2] + 4.5f * v[3];
	const real B =  9    * v[0] - 22.5f * v[1] + 18    * v[2] - 4.5f * v[3];
	const real C = -5.5f * v[0] +  9    * v[1] -  4.5f * v[2] +        v[3];
	const real D =         v[0];											// coefficients of cubic: A*T^3+B*T^2+C*T+D  (where T=(t-t0)/(t1-t0))

#define SIGNCHANGE(a,b) (((a)>=0 && (b)<=0) || ((a)<=0 && (b)>=0))

	// look for first subinterval of [0,1] with a sign change in the cubic
	discrim=4*B*B-12*A*C;
	if(discrim<=0){ // if the cubic is strictly monotonic
		if(SIGNCHANGE(v[0],v[3])){
			T0=0; T1=1;
			V0=v[0]; V1=v[3];
		}else
			return -1;
	}else{ // we need to divide up the cubic into chunks to test for roots
		real rootDiscrim=sqrt(discrim), q;
		if(B<0) q=-.5*(B-rootDiscrim);
		else q=-.5*(B+rootDiscrim);
		real E1=q/A, E2=C/q;
		if(E2<E1) std::swap(E1,E2);
		// look for the intervals we need to check
		if(E1<0){
			if(E2<0){ // check [0,1]
				if(SIGNCHANGE(v[0],v[3])){
					T0=0; T1=1;
					V0=v[0]; V1=v[3];
				}else
					return -1;
			}else{
				if(E2<1){ // check [0,E2] and [E2,1]
					real E2v=D+E2*(C+E2*(B+E2*A));
					if(abs(E2v)<tol){
						return TM = E2;
					}
					if(SIGNCHANGE(v[0],E2v)){
						T0=0; T1=E2;
						V0=v[0]; V1=E2v;
					}else if(SIGNCHANGE(E2v,v[3])){
						T0=E2; T1=1;
						V0=E2v; V1=v[3];
					}else
						return -1;
				}else{ // check [0,1]
					if(SIGNCHANGE(v[0],v[3])){
						T0=0; T1=1;
						V0=v[0]; V1=v[3];
					}else
						return -1;
				}
			}
		}else{
			if(E2<1){ // check [0,E1] [E1,E2] and [E2,1]
				real E1v=D+E1*(C+E1*(B+E1*A));
				if(abs(E1v)<tol){
					return TM = E1;
				}
				if(SIGNCHANGE(v[0],E1v)){
					T0=0; T1=E1;
					V0=v[0]; V1=E1v;
				}else{
					real E2v=D+E2*(C+E2*(B+E2*A));
					if(abs(E2v)<tol){
						return TM = E2;
					}
					if(SIGNCHANGE(E1v,E2v)){
						T0=E1; T1=E2;
						V0=E1v; V1=E2v;
					}else if(SIGNCHANGE(E2v,v[3])){
						T0=E2; T1=1;
						V0=E2v; V1=v[3];
					}else
						return -1;
				}
			}else{
				if(E1<1){ // check [0,E1] and [E1,1]
					real E1v=D+E1*(C+E1*(B+E1*A));
					if(abs(E1v)<tol){
						return TM = E1;
					}
					if(SIGNCHANGE(v[0],E1v)){
						T0=0; T1=E1;
						V0=v[0]; V1=E1v;
					}else if(SIGNCHANGE(E1v,v[3])){
						T0=E1; T1=1;
						V0=E1v; V1=v[3];
					}else
						return -1;
				}else{ // check [0,1]
					if(SIGNCHANGE(v[0],v[3])){
						T0=0; T1=1;
						V0=v[0]; V1=v[3];
					}else
						return -1;
				}
			}
		}
	}

	// early exit if we don't care exactly where the intersection is
	if(!need_exact_time) return 1;	//return true;

	// now do a few iterations of secant search
	if(V0 < 0) {
		for(int i=0; i<5; ++i) {
			alpha = V1/(V1-V0);
			TM = alpha*T0+(1-alpha)*T1;
			VM = D+TM*(C+TM*(B+TM*A));
			
			if(abs(VM) < tol)	{ break;		}
			else if(VM > 0)		{ T1=TM; V1=VM; }
			else				{ T0=TM; V0=VM;	}
		}
	}
	else {
		for(int i=0; i<5; ++i) {
			alpha = V1/(V1-V0);
			TM = alpha*T0+(1-alpha)*T1;
			VM = D+TM*(C+TM*(B+TM*A));

			if(abs(VM) < tol)	{ break;		}
			else if(VM < 0)		{ T1=TM; V1=VM; }
			else				{ T0=TM; V0=VM; }
		}
	}
	return TM;
}