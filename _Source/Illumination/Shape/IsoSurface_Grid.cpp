#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_Grid.h>

using namespace std;
using namespace illuminationspace;

void IsoSurface_Grid::IsoSurface_GridValue::init() {
	for(size_t i=0;i<3;++i) {
		dim[i] = 0; celldim[i] = -1;
	}
	bounds.init();
	Width = InvWidth = Vector3<real>(0);

	for(size_t i=0;i<3;++i)
		for(size_t j=0;j<2;++j)
			isBoundarySurface[i][j] = false;
}

void IsoSurface_Grid::IsoSurface_GridValue::set(const int _dim[3], const BBox _bounds, const bool bdrySurface[3][2]) {
	for(size_t i=0;i<3;++i) {
		dim[i] = _dim[i]; celldim[i] = _dim[i] - 1;		
	}

	bounds = _bounds;
	Width = bounds.extent() / Vector3<real>(celldim[0], celldim[1], celldim[2]);
	InvWidth = static_cast<real>(1.0) / Width;

	for(size_t i=0;i<3;++i)
		for(size_t j=0;j<2;++j)
			isBoundarySurface[i][j] = bdrySurface[i][j];
}

