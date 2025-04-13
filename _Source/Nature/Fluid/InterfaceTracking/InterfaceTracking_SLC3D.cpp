#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_SLC3D.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void InterfaceTracking_SLC3D::set(double _lc, double _uc, int _max_lv) {
	lc = _lc; uc = _uc; max_level = _max_lv;
	double L = uc - lc;	

	s = new slcspace::SlcSurface();
	s0 = new slcspace::SlcSurface();	

	MacGrid3D::set(slcspace::pow2(max_level), slcspace::pow2(max_level), slcspace::pow2(max_level), L / slcspace::pow2(max_level));

	eye = Vector3f(0.5f*L, 0.5f*L, L + L * 1.3f);
	ref = Vector3f(0.5f*L, 0.5f*L, 0.5f*L);	
	up  = Vector3f(0.0f, 1.0f, 0.0f);
	theta = 45.0f;	

	log_info << "InterfaceTracking_SLC3D::set()" << endl;
	log_info << "(LC, UC, Max_Level, Dim) = " << lc << ", " << uc << ", " << max_level << ", " << slcspace::pow2(max_level) << endl;
}

void InterfaceTracking_SLC3D::setPhi(const DistanceField3D<real>& field) {
	class DistField3D_SLC {
	public : 
		DistField3D_SLC(const DistanceField3D<real>& f) : field(f) {}
		const double operator()(const slcspace::SmVector3& _x) const { return field.getPhi(arithmetic::Vector3<real>(_x[0], _x[1], _x[2])); }

		const DistanceField3D<real>& field;
	};

	s->tree.buildTree(DistField3D_SLC(field), lc, uc, max_level);
	s->contourTree(DistField3D_SLC(field));
	s->redistance();
}

void InterfaceTracking_SLC3D::glDrawInterface(const Material& mat, bool filled) const {
	
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glDisable(GL_COLOR_MATERIAL);
	mat.glBindMaterial();	

	glBegin( GL_TRIANGLES );
	
	for(unsigned int i=0;i<s->triangles.numTriangles();++i)	{
		unsigned int i0 = s->triangles[i].a;
		unsigned int i1 = s->triangles[i].b;
		unsigned int i2 = s->triangles[i].c;

		const slcspace::SmVector3 *V0, *V1, *V2, *N0, *N1, *N2;
		V0 = &s->meshPts[i0]; V1 = &s->meshPts[i1]; V2 = &s->meshPts[i2];
		N0 = &s->normals[i0]; N1 = &s->normals[i1]; N2 = &s->normals[i2];

		glNormal3f((*N0)[0], (*N0)[1], (*N0)[2]);
		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glNormal3f((*N1)[0], (*N1)[1], (*N1)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glNormal3f((*N2)[0], (*N2)[1], (*N2)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();
}

void InterfaceTracking_SLC3D::glDrawBoundingBox(const RGBA& c, GLfloat width) const {
	MacGrid3D<double>::glDrawBoundingBox(c, width);
}

void InterfaceTracking_SLC3D::getPhi(Grid3<double,double>& _phi) const {
	throwError("Not Implemented _ InterfaceTracking_SLC3D::getPhi(Grid3<double,double>& _phi)");
}

const bool InterfaceTracking_SLC3D::load(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);
	
	if(s != NULL) 
		delete s;		

	s = new slcspace::SlcSurface();
	s->read(ifs);

	int max_level = s->tree.max_level;
	double lc = s->tree.lc;
	double uc = s->tree.uc;
	double L = uc - lc;

	MacGrid3D::set(slcspace::pow2(max_level), slcspace::pow2(max_level), slcspace::pow2(max_level), L / slcspace::pow2(max_level));

	eye = Vector3f(0.5f*L, 0.5f*L, L + L * 1.3f); 
	ref = Vector3f(0.5f*L, 0.5f*L, 0.5f*L);		
	up = Vector3f(0.0f, 1.0f, 0.0f);								
	theta = 45.0f;	

	return true;
}

const bool InterfaceTracking_SLC3D::save(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);

	s->write(ofs);

	return true;
}