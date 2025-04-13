#include <stdafx.h>
#include <Basis/FileIO/FileIO.h>
#include <Nature/Fluid/CFDsolver/CFDsolver2D.h>
#include <Nature/Fluid/FluidLibrary/VelocityFields.h>

using namespace fluidspace;
using namespace basis;

const real CFDsolver2D::getCFL(real dt, real _h) const {
	MacGrid2D_VelocityGetter<real> field(u,v,TYPE_S);

	Vector2<real> vel;
	real max_absVel = 0, absVel = 0;

	for(index_type i=1;i<nX+1;++i) { for(index_type j=1;j<nY+1;++j) { 

		vel = field.getVel(i,j);
		absVel = abs(vel[0]) + abs(vel[1]);

		if(max_absVel < absVel) 
			max_absVel = absVel;
	} }

	return dt * max_absVel / _h;
}

void CFDsolver2D::setSingleVortexFlow() {
	if(Lx != Ly) {
		cout << "ERROR _ Vortex Flow has square domain" << endl;
		return;
	}
	VortexFlowVelocityField2D<real> velField(Lx);
	setVelocityField(velField);
}

void CFDsolver2D::setRotationFlow(size_t nPerOneRotation)	{
	RotationVelocityField2D<real> velField(nPerOneRotation, Lx, Ly);
	setVelocityField(velField);	
}

const bool CFDsolver2D::load(const string& filename) {
	ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ifs)
		return false;

	string name;
	StringLibrary::importString_Binary(name, ifs);

	isKinematicBoundary.importBinary(ifs);
	u.importBinary(ifs);
	v.importBinary(ifs);
	pressure.importBinary(ifs);

	u0.resize(u); u0 = u;
	v0.resize(v); v0 = v;

	MacGrid2D<real>::set(pressure.nX-2, pressure.nY-2, pressure.dx);

	::importBinary(ifs, Lx);
	::importBinary(ifs, Ly);
	::importBinary(ifs, h);
	::importBinary(ifs, hInv);
	::importBinary(ifs, h2);

	return true;
}

const bool CFDsolver2D::save(const string& filename) const {
	ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
	if(!ofs)
		return false;

	StringLibrary::exportString_Binary(XmlElementName().c_str(), ofs);

	isKinematicBoundary.exportBinary(ofs);
	u.exportBinary(ofs);
	v.exportBinary(ofs);
	pressure.exportBinary(ofs);

	::exportBinary(ofs, Lx);
	::exportBinary(ofs, Ly);
	::exportBinary(ofs, h);
	::exportBinary(ofs, hInv);
	::exportBinary(ofs, h2);

	return true;
}