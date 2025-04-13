#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_Levelset2D.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void InterfaceTracking_Levelset2D::set(index_type nX, index_type nY, real h) {
	Levelset2D<real>::set(nX, nY, h);
	phi0.resize(phi);

	log_info << "InterfaceTracking_Levelset2D::set()" << endl;
	log_info << "(nX, nY, h) = " << nX << ", " << nY << ", " << h << endl;
}

void InterfaceTracking_Levelset2D::setPhi(const Grid2r& _phi) {
	phi = phi0 = _phi;
}

void InterfaceTracking_Levelset2D::setPhi(const DistanceField2D<real>& field) {
	Levelset2D<real>::setPhi(field);
	phi0 = phi;
}

void InterfaceTracking_Levelset2D::createDrawableInterface() {
	Levelset2D<real>::createRenderingPrimitive();
}

void InterfaceTracking_Levelset2D::createDrawableInterface(const Grid2b& isbdry) {
	Levelset2D<real>::createRenderingPrimitive(isbdry);
}

void InterfaceTracking_Levelset2D::glDrawInterface(const Material& mat, bool filled) const {
	if(filled)	Levelset2D<real>::glDrawLiquid(mat.getDiffuse());
	else		Levelset2D<real>::glDrawInterfaceEdges(mat.getDiffuse());
}

void InterfaceTracking_Levelset2D::glDrawPhi(const RGBA &colorP, const RGBA &colorN, real band, GLfloat pointsize) const {
	Levelset2D<real>::glDrawPhi(colorP, colorN, band, pointsize);
}

void InterfaceTracking_Levelset2D::glDrawGrid(const RGBA& c, GLfloat linewidth) const {
	Levelset2D<real>::glDrawUniformGrid(c, linewidth);
}

void InterfaceTracking_Levelset2D::glDrawBoundingBox(const RGBA& c, GLfloat width) const {
	Levelset2D<real>::glDrawBoundingBox(c, width);
}

void InterfaceTracking_Levelset2D::getNormal(Grid2<Vector2<real>,real>& _normal) const {
	if(_normal.nX == 0 && _normal.nY == 0) {
		_normal.resize(phi);
		MacGridLibrary<real>::getNormals(_normal, phi);
		return;
	}

	if(phi.nX == _normal.nX && phi.nY == _normal.nY) {		
		MacGridLibrary<real>::getNormals(_normal, phi);
		return;
	}

	throwError("Not Implemented _ InterfaceTracking_Levelset2D::getNormals(Grid2<real,real>& _normal)");
	return;
}

void InterfaceTracking_Levelset2D::getCurvature(Grid2<real,real>& _curvature) const {
	if(_curvature.nX == 0 && _curvature.nY == 0) {
		_curvature.resize(phi);
		MacGridLibrary<real>::getCurvatures(_curvature, phi);
		return;
	}

	if(phi.nX == _curvature.nX && phi.nY == _curvature.nY) {		
		MacGridLibrary<real>::getCurvatures(_curvature, phi);
		return;
	}

	throwError("Not Implemented _ InterfaceTracking_Levelset2D::getCurvature(Grid2<real,real>& _curvature)");
	return;
}

void InterfaceTracking_Levelset2D::getPhi(Grid2<real,real>& _phi) const {
	if(_phi.nX == 0 && _phi.nY == 0) {
		_phi.resize(phi);
		_phi = phi;
		return;
	}
	if(phi.nX == _phi.nX && phi.nY == _phi.nY) {
		_phi = phi;
		return;
	}

#ifdef USE_OPENMP
	int my_end = int(_phi.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = _phi.nX-1;
	index_type i;
#endif
	for(i=1;i<my_end;++i) 	
		for(index_type j=1;j<_phi.nY-1;++j) 
			_phi(i,j) = phi(_phi.pos(i,j));			

	_phi.setNeumannBoundary();
}