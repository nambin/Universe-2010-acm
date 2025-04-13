#include <stdafx.h>
#include <Nature/MassSpring/PostProcessingDynamics/PostProcessingDynamics.h>

using namespace std;
using namespace springspace;

ProvotLimiting::ProvotLimiting(bool stretch, bool shear, bool bending) {
	minlimit = 0.9f;
	maxlimit = 1.1f;
	maxIter = 200;
	modifyStretch = stretch;
	modifyShear = shear;
	modifyBending = bending;
}

void ProvotLimiting::setParameters_StrainLimiting(springspace::real _min, springspace::real _max, std::size_t _maxIter) {
	minlimit = _min;
	maxlimit = _max;
	maxIter = _maxIter;
}

void ProvotLimiting::postProcess(real dt, TriDynamicMesh<real>& mesh) const {
	std::vector<Array<Spring<real>>*> springs;

	if(modifyStretch)	springs.push_back(&mesh.stretchSpring);

	real maxStr = strainLimiting(mesh.pos, mesh.constraint, springs);

	real dtInv = 1.0f / dt;
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.vel[i] = (mesh.pos[i] - mesh.pos0[i]) * dtInv;
}

void ProvotLimiting::postProcess(real dt, QuadDynamicMesh<real>& mesh) const {
	std::vector<Array<Spring<real>>*> springs;

	if(modifyStretch)	springs.push_back(&mesh.stretchSpring);
	if(modifyShear)		springs.push_back(&mesh.shearSpring);
	if(modifyBending)	springs.push_back(&mesh.bendingSpring);

	real maxStr = strainLimiting(mesh.pos, mesh.constraint, springs);

	real dtInv = 1.0f / dt;
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.vel[i] = (mesh.pos[i] - mesh.pos0[i]) * dtInv;

//	cout << "Max Strain : " << maxStr << endl;
}

const real ProvotLimiting::strainLimiting(Vec3Array<real>& pos, const Constraints<real>& constr, const vector<Array<Spring<real>>*>& springss) const {
	real minlen, maxlen, len, inclen, maxStrain;
	Vector3<real> v2v;	
	bool isConverge;

	for(size_t iter=0;iter<maxIter;++iter) {
		isConverge = true;				
		maxStrain = 0;

		for(size_t s=0;s<springss.size();++s) {

			const Array<Spring<real>>& springs = (*springss[s]);

			for(size_t i=0;i<springs.size();++i) {
				const Spring<real>& s = springs[i];

				v2v = pos[s.idx1] - pos[s.idx0];
				len = mag(v2v);
				minlen = s.restlen * minlimit;
				maxlen = s.restlen * maxlimit;

				if(len < minlen) 		inclen = minlen - len;
				else if(maxlen < len)	inclen = maxlen - len;
				else					continue;

				isConverge = false;	
				maxStrain = basis::MAX(maxStrain, abs((len - s.restlen) / s.restlen));

				v2v /= len;

				if(constr.isConstrained(s.idx0) && constr.isConstrained(s.idx1)) continue;
				else if(constr.isConstrained(s.idx0))	pos[s.idx1] += (inclen * v2v);
				else if(constr.isConstrained(s.idx1))	pos[s.idx0] -= (inclen * v2v);
				else {
					pos[s.idx1] += (0.5f * inclen * v2v);
					pos[s.idx0] -= (0.5f * inclen * v2v);
				}
			}
		}		

		if(isConverge)
			return maxStrain;

//		cout << maxStrain << " ";
	}

	return maxStrain;
}
