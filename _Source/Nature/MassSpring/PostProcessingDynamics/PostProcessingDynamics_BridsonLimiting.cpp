#include <stdafx.h>
#include <Nature/MassSpring/PostProcessingDynamics/PostProcessingDynamics.h>

using namespace std;
using namespace springspace;

BridsonLimiting::BridsonLimiting(bool stretch, bool shear, bool bending) {
	minlimitStrain = 0.9f; 
	maxlimitStrain = 1.1f;	maxIter_Strain = 20;
	limitStrainRate = 0.1f;	maxIter_StrainRate = 20;
	modifyStretch = stretch;
	modifyShear = shear;
	modifyBending = bending;
}

void BridsonLimiting::setParameters_StrainLimiting(real _min, real _max, size_t _maxIter) {
	minlimitStrain = _min;
	maxlimitStrain = _max;
	maxIter_Strain = _maxIter;	
}

void BridsonLimiting::setParameters_StrainRateLimiting(real _limit, size_t _maxIter) {
	limitStrainRate = _limit;
	maxIter_StrainRate = _maxIter;
}

void BridsonLimiting::postProcess(real dt, TriDynamicMesh<real>& mesh) const {
	std::vector<Array<Spring<real>>*> springs;

	if(modifyStretch)	springs.push_back(&mesh.stretchSpring);

	real maxStr = strainLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);
	real maxStrRate = strainRateLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);

	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;
}

void BridsonLimiting::postProcess(real dt, QuadDynamicMesh<real>& mesh) const {
	std::vector<Array<Spring<real>>*> springs;

	if(modifyStretch)	springs.push_back(&mesh.stretchSpring);
	if(modifyShear)		springs.push_back(&mesh.shearSpring);
	if(modifyBending)	springs.push_back(&mesh.bendingSpring);

	real maxStr = strainLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);
	real maxStrRate = strainRateLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);

	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;
}

// Gauss-Seidal
const real BridsonLimiting::strainLimiting(Vec3Array<real>& pos, Vec3Array<real>& vel, real dt, const Constraints<real>& ctrt, const vector<Array<Spring<real>>*>& ss) const {
	real minlen, maxlen, len, inclen, maxStrain;
	Vector3<real> v2v;	
	bool isConverge;

	Array<Vector3<real>>	impulse(vel.size());
	Array<int>				nImpulse(vel.size());

	for(size_t iter=0;iter<maxIter_Strain;++iter) {
		isConverge = true;				
		maxStrain = 0;
		impulse.memsetZero();
		nImpulse.memsetZero();

		for(size_t s=0;s<ss.size();++s) {

			const Array<Spring<real>>& springs = (*ss[s]);

			for(size_t i=0;i<springs.size();++i) {
				const Spring<real>& s = springs[i];

				v2v = pos[s.idx1] - pos[s.idx0] + (vel[s.idx1] - vel[s.idx0]) * dt + (impulse[s.idx1] - impulse[s.idx0]);
				len = mag(v2v);
				minlen = s.restlen * minlimitStrain;
				maxlen = s.restlen * maxlimitStrain;

				maxStrain = basis::MAX(maxStrain, abs((len - s.restlen) / s.restlen));

				if(len < minlen) 		inclen = minlen - len;
				else if(maxlen < len)	inclen = maxlen - len;
				else					continue;

				isConverge = false;					

				v2v /= len;

				if(ctrt.isConstrained(s.idx0) && ctrt.isConstrained(s.idx1)) continue;
				else if(ctrt.isConstrained(s.idx0))	{ impulse[s.idx1] += (inclen * v2v); nImpulse[s.idx1]++; }
				else if(ctrt.isConstrained(s.idx1))	{ impulse[s.idx0] -= (inclen * v2v); nImpulse[s.idx0]++; }
				else {
					impulse[s.idx1] += (0.5f * inclen * v2v); nImpulse[s.idx1]++;
					impulse[s.idx0] -= (0.5f * inclen * v2v); nImpulse[s.idx0]++;
				}		
			}
		}		

		if(isConverge) {
//			cout << "Max Strain 01 : " << maxStrain << " , Iter : " << iter << endl;
			return maxStrain;
		}

		for(size_t i=0;i<vel.size();++i) {
			//if(nImpulse[i] > 0) vel[i] += (impulse[i] / (nImpulse[i] * dt));
			vel[i] += (impulse[i] / dt);
		}

//		if(iter==0) cout << "Max Strain 02 : " << maxStrain << endl;
	}

//	cout << "Max Strain 03 : " << maxStrain << endl;

	return maxStrain;
}

// Gauss-Seidal
const real BridsonLimiting::strainRateLimiting(Vec3Array<real>& pos, Vec3Array<real>& vel, real dt, const Constraints<real>& ctrt, const vector<Array<Spring<real>>*>& ss) const {
	real minlen, maxlen, len, currentlen, inclen, maxStrainRate;
	Vector3<real> v2v;	
	bool isConverge;

	Array<Vector3<real>>	impulse(vel.size());
	Array<int>				nImpulse(vel.size());

	for(size_t iter=0;iter<maxIter_Strain;++iter) {
		isConverge = true;				
		maxStrainRate = 0;
		impulse.memsetZero();
		nImpulse.memsetZero();

		for(size_t s=0;s<ss.size();++s) {

			const Array<Spring<real>>& springs = (*ss[s]);

			for(size_t i=0;i<springs.size();++i) {
				const Spring<real>& s = springs[i];

				v2v = pos[s.idx1] - pos[s.idx0];
				currentlen = mag(v2v);

				v2v += (vel[s.idx1] - vel[s.idx0]) * dt + (impulse[s.idx1] - impulse[s.idx0]);
				len = mag(v2v);

				minlen = currentlen * (1 - limitStrainRate);
				maxlen = currentlen * (1 + limitStrainRate);

				maxStrainRate = basis::MAX(maxStrainRate, abs((len - currentlen) / currentlen));

				if(len < minlen) 		inclen = minlen - len;
				else if(maxlen < len)	inclen = maxlen - len;
				else					continue;

				isConverge = false;			

				v2v /= len;

				if(ctrt.isConstrained(s.idx0) && ctrt.isConstrained(s.idx1)) continue;
				else if(ctrt.isConstrained(s.idx0))	{ impulse[s.idx1] += (inclen * v2v); nImpulse[s.idx1]++; }
				else if(ctrt.isConstrained(s.idx1))	{ impulse[s.idx0] -= (inclen * v2v); nImpulse[s.idx0]++; }
				else {
					impulse[s.idx1] += (0.5f * inclen * v2v); nImpulse[s.idx1]++;
					impulse[s.idx0] -= (0.5f * inclen * v2v); nImpulse[s.idx0]++;
				}		
			}
		}		

		if(isConverge) {
//			cout << "Max Strain Rate 01 : " << maxStrainRate << " , Iter : " << iter << endl;
			return maxStrainRate;
		}

		for(size_t i=0;i<vel.size();++i) {
			//if(nImpulse[i] > 0) vel[i] += (impulse[i] / (nImpulse[i] * dt));
			vel[i] += (impulse[i] / dt);
		}

//		if(iter==0) cout << "Max Strain Rate 02 : " << maxStrainRate << endl;
	}

//	cout << "Max Strain Rate 03 : " << maxStrainRate << endl;

	return maxStrainRate;
}