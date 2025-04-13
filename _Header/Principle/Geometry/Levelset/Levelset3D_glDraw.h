#pragma once

namespace geometryspace {

template<typename T>
void Levelset3D<T>::glDrawSurfaceMesh(const Material& mat) const {
	mesh.glDraw_SmoothSurface(mat);
}

};			// geometryspace