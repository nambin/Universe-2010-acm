#include <stdafx.h>
#include <fstream>
#include <Basis/vectorLibrary/vectorLibrary.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_SRL2D.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void InterfaceTracking_SRL2D::release() { 
	if(srl) delete srl; srl = NULL;
	if(srl0) delete srl0; srl0 = NULL;

	narrow_region.clear();	
}

void InterfaceTracking_SRL2D::set(index_type nX, index_type nY, real h, real narrowband, bool uniform_quadrature) {
	release();

	srl = new SRL2D<real>();	srl->set(nX, nY, h, narrowband, uniform_quadrature);
	srl0 = new SRL2D<real>();	srl0->set(nX, nY, h, narrowband, uniform_quadrature);

	calc_narrow_region(narrowband);

	log_info << "InterfaceTracking_SRL2D::set()" << endl;
	log_info << "(nX, nY, h) = " << nX << ", " << nY << ", " << h << endl;
	log_info << "NarrowBand : " << narrowband << " _ Uniform Quadrature : " << uniform_quadrature << endl;
}

void InterfaceTracking_SRL2D::getCurvature(Grid2<real,real>& curvature) const {
	if(curvature.nX == 0 && curvature.nY == 0) {
		curvature.resize(srl->phi);		
		MacGridLibrary<real>::getCurvatures(curvature, srl->phi);
		return;
	}
	if(srl->phi.nX == curvature.nX && srl->phi.nY == curvature.nY) {
		MacGridLibrary<real>::getCurvatures(curvature, srl->phi);
		return;
	}
}

void InterfaceTracking_SRL2D::getPhi(Grid2<real,real>& _phi) const {
	if(_phi.nX == 0 && _phi.nY == 0) {
		_phi.resize(srl->phi);
		_phi = srl->phi;
		return;
	}
	if(srl->phi.nX == _phi.nX && srl->phi.nY == _phi.nY) {
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
			_phi(i,j) = srl->phi(_phi.pos(i,j));			

	_phi.setNeumannBoundary();
}

void InterfaceTracking_SRL2D::calc_narrow_region(real narrowband) {

	typedef Grid2<bool,real>::index_type index_type;

	narrow_region.clear();

	index_type narrow = narrowband + 2;
	index_type _size = narrow * 2 + 1;

	for(int i=-narrow;i<=narrow;++i) {
		for(int j=-narrow;j<=narrow;++j) {
			if( ( real(SQ(i-0) + SQ(j-0)) <= narrowband * narrowband ) 
			 || ( real(SQ(i-1) + SQ(j-0)) <= narrowband * narrowband ) 
			 || ( real(SQ(i-0) + SQ(j-1)) <= narrowband * narrowband ) 
			 || ( real(SQ(i-1) + SQ(j-1)) <= narrowband * narrowband ) )

			{
				narrow_region.push_back(i);
				narrow_region.push_back(j);
			}
		}
	}

//	cout << "InterfaceTracking_SRL2D _ Narrow Region for " << narrowband << " narrow band" << endl;
//	for(size_t i=0;i<narrow_region.size();i+=2)
//		cout << narrow_region[i] << " " << narrow_region[i+1] << endl;
}

void InterfaceTracking_SRL2D::createDrawableInterface() {
	srl->createRenderingPrimitive(isbdry);
}

void InterfaceTracking_SRL2D::createDrawableInterface(const Grid2b& isbdry) {
	srl->createRenderingPrimitive(isbdry);
}

void InterfaceTracking_SRL2D::glDrawInterface(const Material& mat, bool filled) const {
	if(filled)	srl->glDrawLiquid(mat.getDiffuse());
	else		srl->glDrawInterfaceEdges(mat.getDiffuse());
}

void InterfaceTracking_SRL2D::glDrawPhi(const RGBA &cp, const RGBA &cn, real band, GLfloat ptsize) const {
	srl->glDrawPhi(cp, cn, band, ptsize);
}

void InterfaceTracking_SRL2D::glDrawPhi_SubGrid(const RGBA& cp, const RGBA &cn, real band, GLfloat ptsize) const {
	srl->glDrawPhi_SubGrid(cp, cn, band, ptsize);
}

void InterfaceTracking_SRL2D::glDrawGrid(const RGBA& c, GLfloat width) const {
	srl->glDrawGrid(c, width);
}

void InterfaceTracking_SRL2D::glDrawBoundingBox(const RGBA& c, GLfloat width) const {
	srl->glDrawBoundingBox(c, width);
}

void InterfaceTracking_SRL2D::glDrawAuxiliary() const {
	srl->glDrawSubGrid(RGBA::white(), 1.5f);
}

const bool InterfaceTracking_SRL2D::load(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(XmlElementName().compare(name) != 0)
		return false;

	if(srl == NULL) srl = new SRL2D<real>();
	if(srl0 == NULL) srl0 = new SRL2D<real>();

	isbdry.importBinary(ifs);
	srl->importBinary(ifs);
	srl0->importBinary(ifs);

	calc_narrow_region(srl->narrowband);

	((InterfaceTracking_SRL2D*)(this))->createDrawableInterface(isbdry);

	return true;
}

const bool InterfaceTracking_SRL2D::save(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);	

	isbdry.exportBinary(ofs);
	srl->exportBinary(ofs);
	srl0->exportBinary(ofs);

	return true;
}

const bool InterfaceTracking_SRL2D::load_compressed(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	if(XmlElementName().compare(name) != 0)
		return false;

	if(srl == NULL) srl = new SRL2D<real>();
	srl->importBinary_Compressed(ifs);
	
	((InterfaceTracking_SRL2D*)(this))->createDrawableInterface();

	return true;
}

const bool InterfaceTracking_SRL2D::save_compressed(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);	
	srl->exportBinary_Compressed(ofs);

	return true;
}

// Need to Parallelize SRL
void InterfaceTracking_SRL2D::update_InternalBoundary(const InternalBoundary_Implicit2D* boundary) {
	if(boundary)
		isbdry.init(0,0,1,1);

	isbdry.resize(boundary->phi);
	isbdry.set_boundary(true);

	for(index_type i = 1 ; i < isbdry.nX-1; ++i ) {
		for(index_type j=1;j<isbdry.nY-1;++j) {
			isbdry(i,j) = (boundary->phi(i,j) > 0);
		}
	}


}