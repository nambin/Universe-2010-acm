#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/ParticleLevelset_ErrorCorrection.h>
#include <Nature/Fluid/InterfaceTracking/ParticleLevelset_Mnp.h>

using namespace std;
using namespace fluidspace;

void PLS_ErrorCorrection::errorCorrect(InterfaceTracking_PLS2D& pls, bool setNeumann) const {
	PLS_ParticleMnp::detect_EscapedParticles(pls);

	Grid2r phiPos, phiNeg;
	Grid2b isModified;

	phiPos.resize(pls.phi);	phiNeg.resize(pls.phi);
	isModified.resize(pls.phi);

	phiPos = phiNeg = pls.phi;
	isModified.zero();

	vector<Vector2<grid_index_type>> nIdx;	nIdx.resize(4);
	nIdx.push_back(Vector2<grid_index_type>(0,0)); nIdx.push_back(Vector2<grid_index_type>(0,1));
	nIdx.push_back(Vector2<grid_index_type>(1,0)); nIdx.push_back(Vector2<grid_index_type>(1,1));

	getCorrectedPhi(phiPos, phiNeg, isModified, pls, nIdx);

	for(grid_index_type i=1;i<pls.phi.nX-1;++i) {
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) {
//			if(!isModified(i,j))
//				continue;

			pls.phi(i,j) = (abs(phiPos(i,j)) < abs(phiNeg(i,j))) ? phiPos(i,j) : phiNeg(i,j);
		}
	}

	if(setNeumann)
		pls.phi.setNeumannBoundary();
}

void PLS_ErrorCorrection::errorCorrect(InterfaceTracking_PLS3D& pls, bool setNeumann) const {
	PLS_ParticleMnp::detect_EscapedParticles(pls);

	Grid3r phiPos, phiNeg;
	Grid3b isModified;

	phiPos.resize(pls.phi);	phiNeg.resize(pls.phi);
	isModified.resize(pls.phi);

	phiPos = phiNeg = pls.phi;
	isModified.zero();

	vector<Vector3<grid_index_type>> nIdx;	nIdx.resize(8);
	nIdx.push_back(Vector3<grid_index_type>(0,0,0)); nIdx.push_back(Vector3<grid_index_type>(0,1,0));
	nIdx.push_back(Vector3<grid_index_type>(1,0,0)); nIdx.push_back(Vector3<grid_index_type>(1,1,0));
	nIdx.push_back(Vector3<grid_index_type>(0,0,1)); nIdx.push_back(Vector3<grid_index_type>(0,1,1));
	nIdx.push_back(Vector3<grid_index_type>(1,0,1)); nIdx.push_back(Vector3<grid_index_type>(1,1,1));

	getCorrectedPhi(phiPos, phiNeg, isModified, pls, nIdx);

	for(grid_index_type i=1;i<pls.phi.nX-1;++i) {
		for(grid_index_type j=1;j<pls.phi.nY-1;++j) {
			for(grid_index_type k=1;k<pls.phi.nZ-1;++k) {
//				if(!isModified(i,j,k))
//					continue;

				pls.phi(i,j,k) = (abs(phiPos(i,j,k)) < abs(phiNeg(i,j,k))) ? phiPos(i,j,k) : phiNeg(i,j,k);
			}			
		}
	}

	if(setNeumann)
		pls.phi.setNeumannBoundary();
}