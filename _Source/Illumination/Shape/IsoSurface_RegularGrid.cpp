#include <stdafx.h>
#include <Principle/Geometry/MacGrid/MacGridLibrary.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_RegularGrid.h>

using namespace std;
using namespace illuminationspace;

IsoSurface_RegularGrid::IsoSurface_RegularGridValue::IsoSurface_RegularGridValue(const Transform& o2w, const Grid3<real,real>& phi, const bool isBdrySurface[3][2]) 
: IsoSurface_Grid::IsoSurface_GridValue(o2w) 
{
	init();
	set(phi, isBdrySurface);
}

void IsoSurface_RegularGrid::IsoSurface_RegularGridValue::init() {
	phi.release();
	normals.release();
	surfaceMarker.release();
}

void IsoSurface_RegularGrid::IsoSurface_RegularGridValue::set(const Grid3<real,real>& grid, const bool bdrySurface[3][2]) {
	IsoSurface_GridValue::set(grid.dim, grid.getAABB(), bdrySurface);

	phi.resize(grid);	
	phi.set(grid);	

	setNormals();
	setSurfaceMarker();
}

void IsoSurface_RegularGrid::IsoSurface_RegularGridValue::setNormals() {
	normals.resize(phi);
	geometryspace::MacGridLibrary<real>::getNormals(normals, phi);
}

void IsoSurface_RegularGrid::IsoSurface_RegularGridValue::setSurfaceMarker() {
	surfaceMarker.init(celldim[0], celldim[1], celldim[2],0,0,0);
	surfaceMarker.zero();

	for(int i=0; i<celldim[0]; ++i) {
		for(int j=0; j<celldim[1]; ++j) {			
			for(int k=0; k<celldim[2]; ++k) {
				if((phi(i,j,k)<0 && phi(i+1,j,k)<0 && phi(i,j+1,k)<0 && phi(i+1,j+1,k)<0 && phi(i,j,k+1)<0 && phi(i+1,j,k+1)<0 && phi(i,j+1,k+1)<0 && phi(i+1,j+1,k+1)<0)
				|| (phi(i,j,k)>0 && phi(i+1,j,k)>0 && phi(i,j+1,k)>0 && phi(i+1,j+1,k)>0 && phi(i,j,k+1)>0 && phi(i+1,j,k+1)>0 && phi(i,j+1,k+1)>0 && phi(i+1,j+1,k+1)>0))
					; // do nothing -- already false
				else
					surfaceMarker(i,j,k) = true;
			}
		}
	}
}

void IsoSurface_RegularGrid::IsoSurface_RegularGridValue::getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const {
	dgShading = dg;
}

const bool IsoSurface_RegularGrid::IsoSurface_RegularGridValue::hit(const Ray& r, real& tHit, DifferentialGeometry& dg) const {
	Ray ray;
	WorldToObject.transformRay(r, &ray);
	return hitCore(*this, ray, &tHit, &dg);
}

const bool IsoSurface_RegularGrid::IsoSurface_RegularGridValue::shadowHit(const Ray& r) const {
	Ray ray;
	WorldToObject.transformRay(r, &ray);
	return hitCore(*this, ray);
}

const bool IsoSurface_RegularGrid::IsoSurface_RegularGridValue::hitCell(const Ray &r, const int voxel[3], real t0, real t1, real *tHit, DifferentialGeometry *dg) const {
	const real v[8] = {	phi(voxel[0]  ,voxel[1]  ,voxel[2]), phi(voxel[0]  ,voxel[1]  ,voxel[2]+1),
						phi(voxel[0]  ,voxel[1]+1,voxel[2]), phi(voxel[0]  ,voxel[1]+1,voxel[2]+1),
						phi(voxel[0]+1,voxel[1]  ,voxel[2]), phi(voxel[0]+1,voxel[1]  ,voxel[2]+1),
						phi(voxel[0]+1,voxel[1]+1,voxel[2]), phi(voxel[0]+1,voxel[1]+1,voxel[2]+1)};

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
