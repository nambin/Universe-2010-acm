#include <stdafx.h>
#include <Nature/Fluid/DensityTracking/DensityTracking.h>

using namespace std;
using namespace basis;
using namespace geometryspace;
using namespace fluidspace;

void DensityTracking2D::set(grid_index_type _nX, grid_index_type _nY, real _h) {
	MacGrid2D<real>::set(_nX, _nY, _h);

	const Vector2<real> cOffset(-0.5f, -0.5f);

	density.init(nX+2, nY+2, h, h, cOffset);
	density0.init(nX+2, nY+2, h, h, cOffset);
}

void DensityTracking2D::set(const Grid2r& field) {
	assert(density.nX == field.nX && density.nY == field.nY);
	density = density0 = field;
}

void DensityTracking2D::set(const DistanceField2D<real>& field) {
	Vector2<real> pos;
	real dist, alpha;
	grid_index_type nX = density.nX, nY = density.nY;
	real dx = density.dx, dy = density.dy;
	real hInv = MIN(1.0f/dx, 1.0f/dy);

	for(grid_index_type x = 0; x < nX; x++) {
		for(grid_index_type y = 0; y < nY; y++) {
			pos = Vector2<real>(dx * x, dy * y);
			dist = field.getPhi(pos);
			alpha = MacGridLibrary<real>::heaviside(dist * hInv, real(0), real(1));	

			density(x,y) = alpha;
		}
	}

	density.setNeumannBoundary();
	density0 = density;
}

const real DensityTracking2D::getTotalDensity() const {
	real _sum = 0;
	for(index_type i=1;i<density.nX-1;++i) 
		for(index_type j=1;j<density.nY-1;++j)
			_sum += density(i,j);
	return _sum * density.dx * density.dy;
}

void DensityTracking3D::set(grid_index_type _nX, grid_index_type _nY, grid_index_type _nZ, real _h) {
	MacGrid3D<real>::set(_nX, _nY, _nZ, _h);

	const Vector3<real> cOffset(-0.5f, -0.5f, -0.5f);

	density.init(nX+2, nY+2, nZ+2, h, h, h, cOffset);
	density0.init(nX+2, nY+2, nZ+2, h, h, h, cOffset);
	isUpdated = true;
}

void DensityTracking3D::set(const Grid3r& field) {
	assert(density.nX == field.nX && density.nY == field.nY && density.nZ == field.nZ);
	density = density0 = field;
	isUpdated = true;
}

void DensityTracking3D::set(const DistanceField3D<real>& field) {
	Vector3<real> pos;
	real dist, alpha;
	grid_index_type nX = density.nX, nY = density.nY, nZ = density.nZ;
	real dx = density.dx, dy = density.dy, dz = density.dz;
	real hInv = MIN3(1.0f/dx, 1.0f/dy, 1.0f/dz);

	for(grid_index_type x = 0; x < nX; x++) {
		for(grid_index_type y = 0; y < nY; y++) {
			for(grid_index_type z = 0; z < nZ; z++) {
				pos = Vector3<real>(dx * x, dy * y, dz * z);
				dist = field.getPhi(pos);
				alpha = MacGridLibrary<real>::heaviside(dist * hInv, real(0), real(1));	

				density(x,y,z) = alpha;
			}
		}
	}

	density.setNeumannBoundary();
	density0 = density;
}

const real DensityTracking3D::getTotalDensity() const {
	real _sum = 0;
	for(index_type i=1;i<density.nX-1;++i) 
		for(index_type j=1;j<density.nY-1;++j)
			for(index_type k=1;k<density.nZ-1;++k)
				_sum += density(i,j,k);
	return _sum * density.dx * density.dy * density.dz;
}