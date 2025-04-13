#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_SRL3D.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

InterfaceTracking_SRL3D::InterfaceTracking_SRL3D(): srl(NULL), srl0(NULL) {}

void InterfaceTracking_SRL3D::release() { 
	if(srl) delete srl; srl = NULL;
	if(srl0) delete srl0; srl0 = NULL;

	narrow_region.clear();
}

void InterfaceTracking_SRL3D::set(index_type nX, index_type nY, index_type nZ, real h, real narrowband, bool uniform_quadrature) {
	release();

	srl = new SRL3D<real>();	srl->set(nX, nY, nZ, h, narrowband, uniform_quadrature);
	srl0 = new SRL3D<real>();	srl0->set(nX, nY, nZ, h, narrowband, uniform_quadrature);

	calc_narrow_region(narrowband);

	log_info << "InterfaceTracking_SRL3D::set()" << endl;
	log_info << "(nX, nY, nZ, h) = " << nX << ", " << nY << ", " << nZ << ", " << h << endl;
	log_info << "NarrowBand : " << narrowband << " _ Uniform Quadrature : " << uniform_quadrature << endl;

	eye = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, srl->Lz + MAX(srl->Lx,srl->Ly) * 1.3f); 
	ref = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, 0.5f*srl->Lz);		
	up = Vector3f(0.0f, 1.0f, 0.0f);								
	theta = 45.0f;		
}

void InterfaceTracking_SRL3D::getPhi(Grid3<real,real>& _phi) const {
	if(_phi.nX == 0 && _phi.nY == 0 && _phi.nZ == 0) {
		_phi.resize(srl->phi);
		_phi = srl->phi;
		return;
	}

	if(srl->phi.nX == _phi.nX && srl->phi.nY == _phi.nY && srl->phi.nZ == _phi.nZ) {
		_phi = srl->phi;
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
				_phi(i,j,k) = srl->phi(_phi.pos(i,j,k));

	_phi.setNeumannBoundary();
}

void InterfaceTracking_SRL3D::getCurvature(Grid3<real,real>& _curvature) const {
	if(_curvature.nX == 0 && _curvature.nY == 0 && _curvature.nZ == 0) {
		_curvature.resize(srl->phi);
		MacGridLibrary<real>::getCurvatures(_curvature, srl->phi);
		return;
	}

	if(srl->phi.nX == _curvature.nX && srl->phi.nY == _curvature.nY && srl->phi.nZ == _curvature.nZ) {		
		MacGridLibrary<real>::getCurvatures(_curvature, srl->phi);
		return;
	}

	throwError("Not Implemented _ InterfaceTracking_SRL3D::getCurvature(Grid3<real,real>& _curvature)");
	return;
}


void InterfaceTracking_SRL3D::calc_narrow_region(real narrowband) {

	typedef Grid3<bool,real>::index_type index_type;

	narrow_region.clear();

	index_type narrow = narrowband + 2;
	index_type _size = narrow * 2 + 1;

	for(int i=-narrow;i<=narrow;++i) {
		for(int j=-narrow;j<=narrow;++j) {
			for(int k=-narrow;k<=narrow;++k) {
				if( ( real(SQ(i-0) + SQ(j-0) + SQ(k-0)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-1) + SQ(j-0) + SQ(k-0)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-0) + SQ(j-1) + SQ(k-0)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-1) + SQ(j-1) + SQ(k-0)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-0) + SQ(j-0) + SQ(k-1)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-1) + SQ(j-0) + SQ(k-1)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-0) + SQ(j-1) + SQ(k-1)) <= narrowband * narrowband ) 
				 || ( real(SQ(i-1) + SQ(j-1) + SQ(k-1)) <= narrowband * narrowband )  )

				{
					narrow_region.push_back(i);
					narrow_region.push_back(j);
					narrow_region.push_back(k);
				}
			}			
		}
	}

	//	cout << "InterfaceTracking_SRL3D _ Narrow Region for " << narrowband << " narrow band" << endl;
	//	for(size_t i=0;i<narrow_region.size();i+=3)
	//		cout << narrow_region[i] << " " << narrow_region[i+1] << " " << narrow_region[i+2] << endl;
}

void InterfaceTracking_SRL3D::createDrawableInterface() {
	srl->createRenderingPrimitive(isbdry);
}

void InterfaceTracking_SRL3D::createDrawableInterface(const Grid3b& isbdry) {
	srl->createRenderingPrimitive(isbdry);
}

void InterfaceTracking_SRL3D::glDrawInterface(const Material& mat, bool filled) const {
	if(filled)	srl->glDrawSurfaceMesh(mat);	
	else		srl->glDrawLineStrip(mat);	
}

void InterfaceTracking_SRL3D::glDrawInnerBoundary(const virtualspace::Material& mat) const {

	if(isbdry.is_empty())
		return;

	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glDisable(GL_COLOR_MATERIAL);
	mat.glBindMaterial();	

	Vector3<real> org;

	for(index_type i=1;i<isbdry.nX-1;++i) {
		for(index_type j=1;j<isbdry.nY-1;++j) {
			for(index_type k=1;k<isbdry.nZ-1;++k) {

				if(isbdry(i,j,k) == false)
					continue;

				org = isbdry.pos(i,j,k);
				
				glTranslatef(org[0], org[1], org[2]);
				glutSolidCube(isbdry.dx);
				glTranslatef(-org[0], -org[1], -org[2]);
			}			
		}
	}
}

void InterfaceTracking_SRL3D::glDrawPhi(const RGBA &cp, const RGBA &cn, real band, GLfloat ptsize) const {
	srl->glDrawPhi(cp, cn, band, ptsize);
}

void InterfaceTracking_SRL3D::glDrawPhi_SubGrid(const RGBA& cp, const RGBA &cn, real band, GLfloat ptsize) const {
	srl->glDrawPhi_SubGrid(cp, cn, band, ptsize);
}

void InterfaceTracking_SRL3D::glDrawGrid(const RGBA& c, GLfloat width) const {
	srl->glDrawGrid(c, width);
}

void InterfaceTracking_SRL3D::glDrawBoundingBox(const RGBA& c, GLfloat width) const {
	srl->glDrawBoundingBox(c, width);
}

void InterfaceTracking_SRL3D::glDrawAuxiliary() const {
	srl->glDrawSubGrid(RGBA::white(), 1.5f);
}

const bool InterfaceTracking_SRL3D::load(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(XmlElementName().compare(name) != 0)
		return false;

	if(srl == NULL) srl = new SRL3D<real>();
	if(srl0 == NULL) srl0 = new SRL3D<real>();

	isbdry.importBinary(ifs);
	srl->importBinary(ifs);
	srl0->set(*srl);

	calc_narrow_region(srl->narrowband);

	eye = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, srl->Lz + MAX(srl->Lx,srl->Ly) * 1.3f); 
	ref = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, 0.5f*srl->Lz);		
	up = Vector3f(0.0f, 1.0f, 0.0f);								
	theta = 45.0f;	

	return true;
}

const bool InterfaceTracking_SRL3D::save(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;
	
	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);

	isbdry.exportBinary(ofs);
	srl->exportBinary(ofs);

	return true;
}

const bool InterfaceTracking_SRL3D::load_compressed(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(XmlElementName().compare(name) != 0)
		return false;

	if(srl == NULL) srl = new SRL3D<real>();
	srl->importBinary_Compressed(ifs);
	isbdry.importBinary(ifs);

	((InterfaceTracking_SRL3D*)(this))->createDrawableInterface();
	
	eye = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, srl->Lz + MAX(srl->Lx,srl->Ly) * 1.3f); 
	ref = Vector3f(0.5f*srl->Lx, 0.5f*srl->Ly, 0.5f*srl->Lz);		
	up = Vector3f(0.0f, 1.0f, 0.0f);								
	theta = 45.0f;	

	return true;
}

const bool InterfaceTracking_SRL3D::save_compressed(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);	
	srl->exportBinary_Compressed(ofs);
	isbdry.exportBinary(ofs);

	return true;
}

// Need to Parallelize SRL
void InterfaceTracking_SRL3D::update_InternalBoundary(const InternalBoundary_Implicit3D* boundary) {
	if(boundary)
		isbdry.init(0,0,0,1,1,1);

	isbdry.resize(boundary->phi);
	isbdry.set_boundary(true);

#ifdef USE_OPENMP	
	int my_end = int(isbdry.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = isbdry.nX-1;
	index_type i;
#endif
	for( i = 1 ; i < my_end; ++i ) {
		for(index_type j=1;j<isbdry.nY-1;++j) {
			for(index_type k=1;k<isbdry.nZ-1;++k) {
				isbdry(i,j,k) = (boundary->phi(i,j,k) > 0);
			}			
		}
	}
}