#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_Levelset3D.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void InterfaceTracking_Levelset3D::set(index_type nX, index_type nY, index_type nZ, real h) {
	Levelset3D<real>::set(nX, nY, nZ, h);
	phi0.resize(phi);

	eye = Vector3f(0.5f*Lx, 0.5f*Ly, Lz + MAX(Lx, Ly) * 1.3f);
	ref = Vector3f(0.5f*Lx, 0.5f*Ly, 0.5f*Lz);	
	up  = Vector3f(0.0f, 1.0f, 0.0f);
	theta = 45.0f;	

	log_info << "InterfaceTracking_Levelset3D::set()" << endl;
	log_info << "(nX, nY, nZ, h) = " << nX << ", " << nY << ", " << nZ << ", " << h << endl;
}

void InterfaceTracking_Levelset3D::setPhi(const Grid3r& _phi) {
	phi = phi0 = _phi;
}

void InterfaceTracking_Levelset3D::setPhi(const DistanceField3D<real>& field) {
	Levelset3D<real>::setPhi(field);
	phi0 = phi;
}

void InterfaceTracking_Levelset3D::createDrawableInterface() {
	Levelset3D<real>::createRenderingPrimitive();
}

void InterfaceTracking_Levelset3D::createDrawableInterface(const Grid3b& isbdry) {
	Levelset3D<real>::createRenderingPrimitive(isbdry);
}

void InterfaceTracking_Levelset3D::glDrawInterface(const Material& mat, bool filled) const {
	if(filled) Levelset3D<real>::glDrawSurfaceMesh(mat);
}

void InterfaceTracking_Levelset3D::glDrawBoundingBox(const RGBA& c, GLfloat width) const {
	Levelset3D<real>::glDrawBoundingBox(c, width);
}

void InterfaceTracking_Levelset3D::getNormal(Grid3<Vector3<real>,real>& _normal) const {
	if(_normal.nX == 0 && _normal.nY == 0 && _normal.nZ == 0) {
		_normal.resize(phi);
		MacGridLibrary<real>::getNormals(_normal, phi);
		return;
	}

	if(phi.nX == _normal.nX && phi.nY == _normal.nY && phi.nZ == _normal.nZ) {		
		MacGridLibrary<real>::getNormals(_normal, phi);
		return;
	}

	throwError("Not Implemented _ InterfaceTracking_Levelset3D::getNormal(Grid3<real,real>& _normal)");
	return;
}

void InterfaceTracking_Levelset3D::getCurvature(Grid3<real,real>& _curvature) const {
	if(_curvature.nX == 0 && _curvature.nY == 0 && _curvature.nZ == 0) {
		_curvature.resize(phi);
		MacGridLibrary<real>::getCurvatures(_curvature, phi);
		return;
	}

	if(phi.nX == _curvature.nX && phi.nY == _curvature.nY && phi.nZ == _curvature.nZ) {		
		MacGridLibrary<real>::getCurvatures(_curvature, phi);
		return;
	}

	throwError("Not Implemented _ InterfaceTracking_Levelset3D::getCurvature(Grid3<real,real>& _curvature)");
	return;
}

void InterfaceTracking_Levelset3D::getPhi(Grid3<real,real>& _phi) const {
	if(_phi.nX == 0 && _phi.nY == 0 && _phi.nZ == 0) {
		_phi.resize(phi);
		_phi = phi;
		return;
	}

	if(phi.nX == _phi.nX && phi.nY == _phi.nY && phi.nZ == _phi.nZ) {
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
			for(index_type k=1;k<_phi.nZ-1;++k) 
				_phi(i,j,k) = phi(_phi.pos(i,j,k));

	_phi.setNeumannBoundary();
}

const bool InterfaceTracking_Levelset3D::load_compressed(const string& filename) { 
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(XmlElementName().compare(name) != 0)
		return false;

	phi.importBinary(ifs);

	return true;
}

const bool InterfaceTracking_Levelset3D::save_compressed(const string& filename) const {

	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);
	phi.exportBinary(ofs);

	return true;
}