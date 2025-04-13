#pragma once

namespace geometryspace {

template<typename T>
void MacGrid3D<T>::set(index_type _nX, index_type _nY, index_type _nZ, T _h)	{
	nX = _nX;
	nY = _nY;
	nZ = _nZ;

	h = _h;	h2 = h*h;
	hInv = 1.0f / h;

	dx = dy = dz = h;
	dxInv = dyInv = dzInv = hInv;
	dx2 = dy2 = dz2 = h2;

	Lx = h * T(nX);	Ly = h * T(nY);	Lz = h * T(nZ);

	numU = (nX+1)*(nY+2)*(nZ+2);
	numV = (nX+2)*(nY+1)*(nZ+2);
	numW = (nX+2)*(nY+2)*(nZ+1);
	numS = (nX+2)*(nY+2)*(nZ+2);
	numCells = nX*nY*nZ;
	numNodes  = (nX+1)*(nY+1)*(nZ+1);	
	numXEdges = nX*(nY+1)*(nZ+1);
	numYEdges = (nX+1)*nY*(nZ+1);
	numZEdges = (nX+1)*(nY+1)*nZ;
}

template<typename T>
const bool MacGrid3D<T>::importBinary(std::ifstream& ifs) {
	::importBinary(ifs, Lx);		::importBinary(ifs, Ly);		::importBinary(ifs, Lz);
	::importBinary(ifs, nX);		::importBinary(ifs, nY);		::importBinary(ifs, nZ);
	::importBinary(ifs, numU);		::importBinary(ifs, numV);		::importBinary(ifs, numW);	::importBinary(ifs, numS);
	::importBinary(ifs, numNodes);	::importBinary(ifs, numCells);
	::importBinary(ifs, numXEdges);	::importBinary(ifs, numYEdges);	::importBinary(ifs, numZEdges);

	::importBinary(ifs, dx);		::importBinary(ifs, dy);		::importBinary(ifs, dz);
	::importBinary(ifs, dxInv);		::importBinary(ifs, dyInv);		::importBinary(ifs, dzInv);
	::importBinary(ifs, dx2);		::importBinary(ifs, dy2);		::importBinary(ifs, dz2);
	::importBinary(ifs, h);			::importBinary(ifs, hInv);		::importBinary(ifs, h2);

	return true;
}

template<typename T>
const bool MacGrid3D<T>::exportBinary(std::ofstream& ofs) const {
	::exportBinary(ofs, Lx);		::exportBinary(ofs, Ly);		::exportBinary(ofs, Lz);
	::exportBinary(ofs, nX);		::exportBinary(ofs, nY);		::exportBinary(ofs, nZ);
	::exportBinary(ofs, numU);		::exportBinary(ofs, numV);		::exportBinary(ofs, numW);	::exportBinary(ofs, numS);
	::exportBinary(ofs, numNodes);	::exportBinary(ofs, numCells);
	::exportBinary(ofs, numXEdges);	::exportBinary(ofs, numYEdges);	::exportBinary(ofs, numZEdges);

	::exportBinary(ofs, dx);		::exportBinary(ofs, dy);		::exportBinary(ofs, dz);
	::exportBinary(ofs, dxInv);		::exportBinary(ofs, dyInv);		::exportBinary(ofs, dzInv);
	::exportBinary(ofs, dx2);		::exportBinary(ofs, dy2);		::exportBinary(ofs, dz2);
	::exportBinary(ofs, h);			::exportBinary(ofs, hInv);		::exportBinary(ofs, h2);

	return true;
}

template<typename T>
void MacGrid3D<T>::glDrawBoundingBox( const RGBA &color, GLfloat width ) const {
	virtualspace::GlDrawLibrary::glDrawCube_LineStrip(Lx, Ly, Lz, color, width);
}

};		// geometryspace