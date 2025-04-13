#include <stdafx.h>
#include <Principle/Geometry/SRL/SRL_Derivative.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_SRLGrid.h>

using namespace std;
using namespace geometryspace;
using namespace illuminationspace;

IsoSurface_SRLGrid::IsoSurface_SRLGridValue::IsoSurface_SRLGridValue(const illuminationspace::Transform& o2w, const SRL3D<real>& srl, const bool isBdrySurface[3][2]) 
: IsoSurface_Grid::IsoSurface_GridValue(o2w), interp(true)  
{
	init();
	set(srl, isBdrySurface);
}

void IsoSurface_SRLGrid::IsoSurface_SRLGridValue::init() {
	srl.release();
	normals.clear();
	surfaceMarker.release();
}

void IsoSurface_SRLGrid::IsoSurface_SRLGridValue::set(const SRL3D<real>& _srl, const bool bdrySurface[3][2]) {
	IsoSurface_GridValue::set(_srl.phi.dim, _srl.phi.getAABB(), bdrySurface);

	srl.set(_srl);

	setNormals();
	setSurfaceMarker();
}

const real IsoSurface_SRLGrid::IsoSurface_SRLGridValue::get_phi(const Vector3<real> &p) const { 
	int i, j, k;
	Vector3<real> x;

	srl.get_idx_normalized_coord(i, j, k, x, p);

	if(srl.is_exist_node(i,j,k))	return lerp(srl, i, j, k, x);
	else							return coarse_lerp(srl, i, j, k, x);
}

const Vector3<real> IsoSurface_SRLGrid::IsoSurface_SRLGridValue::get_normal(const Vector3<real>& p) const {	
	int i, j, k;
	Vector3<real> x;

	srl.get_idx_normalized_coord(i, j, k, x, p);

	return unitVector(lerp(srl, normals[srl.indices(i,j,k)], x));							
}

void IsoSurface_SRLGrid::IsoSurface_SRLGridValue::setNormals() {
	srl.getNormals_SubGrid(normals);	
}

void IsoSurface_SRLGrid::IsoSurface_SRLGridValue::setSurfaceMarker() {
	surfaceMarker.init(celldim[0], celldim[1], celldim[2],0,0,0);
	surfaceMarker.zero();

	for(int i=0; i<celldim[0]; ++i) {
		for(int j=0; j<celldim[1]; ++j) {			
			for(int k=0; k<celldim[2]; ++k) {
				
				if(srl.is_exist_node(i,j,k) == false)
					continue;

				surfaceMarker(i,j,k) = srl.node(i,j,k).has_interface();
			}
		}
	}
}

void IsoSurface_SRLGrid::IsoSurface_SRLGridValue::getShadingGeometry(const illuminationspace::Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const {
	dgShading = dg;
}

const bool IsoSurface_SRLGrid::IsoSurface_SRLGridValue::hit(const Ray& r, real& tHit, DifferentialGeometry& dg) const {
	Ray ray;
	WorldToObject.transformRay(r, &ray);
	return hitCore(*this, ray, &tHit, &dg);
}

const bool IsoSurface_SRLGrid::IsoSurface_SRLGridValue::shadowHit(const Ray& r) const {
	Ray ray;
	WorldToObject.transformRay(r, &ray);
	return hitCore(*this, ray);
}

const bool IsoSurface_SRLGrid::IsoSurface_SRLGridValue::hitCell(const Ray& r, const int voxel[3], real t0, real t1, real* tHit, DifferentialGeometry* dg) const {
	const real v[8] = {	srl.phi(voxel[0]  ,voxel[1]  ,voxel[2]), srl.phi(voxel[0]  ,voxel[1]  ,voxel[2]+1),
						srl.phi(voxel[0]  ,voxel[1]+1,voxel[2]), srl.phi(voxel[0]  ,voxel[1]+1,voxel[2]+1),
						srl.phi(voxel[0]+1,voxel[1]  ,voxel[2]), srl.phi(voxel[0]+1,voxel[1]  ,voxel[2]+1),
						srl.phi(voxel[0]+1,voxel[1]+1,voxel[2]), srl.phi(voxel[0]+1,voxel[1]+1,voxel[2]+1)};

	const real tol = 1e-4 * ( abs(v[0])+abs(v[1])+abs(v[2])+abs(v[3]) + abs(v[4])+abs(v[5])+abs(v[6])+abs(v[7]) );

	const Vector3<real> w0 = (r(t0) - VoxelToPos(voxel)) * InvWidth;
	const Vector3<real> w_icr = (r.direction() * InvWidth) * (t1 - t0);

	const real TM = secant_iteration(v, w0, w_icr, tol, tHit != NULL);

	if(TM < - tol)
		return false;

	if(tHit == NULL) 
		return true;

	(*tHit) = t0 + (t1 - t0) * TM;

	const Vector3<real> pt = r(*tHit);
	const Vector3<real> snormal = get_normal(pt);

	dg->hit_p = ObjectToWorld.transformPoint(pt);
	dg->shadingNormal = ObjectToWorld.transformNormal(snormal);
	dg->uvw.initFromUnitVector(dg->shadingNormal);

	return true;
}