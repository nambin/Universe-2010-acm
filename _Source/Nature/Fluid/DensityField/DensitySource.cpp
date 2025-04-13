#include <stdafx.h>
#include <Nature/Fluid/DensityTracking/DensitySource.h>

using namespace std;
using namespace fluidspace;

void DensitySource_DistanceField::add(DensityTracking* tracking, real dt) const {

	if(field->getDimension() != tracking->getDimension())
		throwError("Error in DensitySource_DistanceField<FIELD>::add(DensityTracking* density, fluidspace::real dt)");

	if(t_start <= time && time < t_end) {
		if(tracking->getDimension() == 2) {
			Grid2r& density = static_cast<DensityTracking2D*>(tracking)->density;

			Vector2<real> pos;
			real dist, alpha;
			grid_index_type nX = density.nX, nY = density.nY;
			real dx = density.dx, dy = density.dy;
			real hInv = MIN(1.0f/dx, 1.0f/dy);

			for(grid_index_type x = 1; x < nX-1; x++) {
				for(grid_index_type y = 1; y < nY-1; y++) {
					pos = Vector2<real>(dx * x, dy * y);
					dist = field->getPhi(pos);
					alpha = MacGridLibrary<real>::heaviside(dist * hInv, real(0), real(1));	

					density(x,y) += (alpha * magnitude * dt);
				}
			}

			density.setNeumannBoundary();
		}
		else {
			Grid3r& density = static_cast<DensityTracking3D*>(tracking)->getDensity();

			Vector3<real> pos;
			real dist, alpha;
			grid_index_type nX = density.nX, nY = density.nY, nZ = density.nZ;
			real dx = density.dx, dy = density.dy, dz = density.dz;
			real hInv = MIN3(1.0f/dx, 1.0f/dy, 1.0f/dz);

			for(grid_index_type x = 1; x < nX-1; x++) {
				for(grid_index_type y = 1; y < nY-1; y++) {
					for(grid_index_type z = 1; z < nZ-1; z++) {
						pos = Vector3<real>(dx * x, dy * y, dz * z);
						dist = field->getPhi(pos);
						alpha = MacGridLibrary<real>::heaviside(dist * hInv, real(0), real(1));	

						density(x,y,z) += (alpha * magnitude * dt);
					}
				}
			}

			density.setNeumannBoundary();
		}
	}

	time += dt;
}