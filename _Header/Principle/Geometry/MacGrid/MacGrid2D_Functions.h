#pragma once

namespace geometryspace {

template<typename T>
void MacGrid2D<T>::set(index_type _nX, index_type _nY, T _h)	{
	nX = _nX;
	nY = _nY;

	h = _h;	h2 = h*h;
	hInv = 1.0f / h;

	Lx = h * nX;	Ly = h * nY;

	numU = (nX+1) * (nY+2);
	numV = (nX+2) * (nY+1);
	numS = (nX+2) * (nY+2);
	numNodes = (nX+1) * (nY+1);
	numEdges = nX * (nY+1) + (nX+1) * nY;	
}

template<typename T>
const bool MacGrid2D<T>::importBinary(std::ifstream& ifs) {
	::importBinary(ifs, Lx);		::importBinary(ifs, Ly);
	::importBinary(ifs, nX);		::importBinary(ifs, nY);
	::importBinary(ifs, numU);		::importBinary(ifs, numV);		::importBinary(ifs, numS);
	::importBinary(ifs, numNodes);	::importBinary(ifs, numEdges);
	::importBinary(ifs, h);			::importBinary(ifs, hInv);		::importBinary(ifs, h2);

	return true;
}

template<typename T>
const bool MacGrid2D<T>::exportBinary(std::ofstream& ofs) const {
	::exportBinary(ofs, Lx);		::exportBinary(ofs, Ly);
	::exportBinary(ofs, nX);		::exportBinary(ofs, nY);
	::exportBinary(ofs, numU);		::exportBinary(ofs, numV);		::exportBinary(ofs, numS);
	::exportBinary(ofs, numNodes);	::exportBinary(ofs, numEdges);
	::exportBinary(ofs, h);			::exportBinary(ofs, hInv);		::exportBinary(ofs, h2);

	return true;
}

template<typename T>
void MacGrid2D<T>::glDrawBoundingBox( const RGBA &color, GLfloat width ) const {

	glLineWidth(width);
	glColor4f( color[0], color[1], color[2], color[3] );

	glBegin( GL_LINES );
	glVertex3f( (GLfloat)0 , (GLfloat)0  , 0.001f);
	glVertex3f( (GLfloat)0 , (GLfloat)Ly , 0.001f);
	glVertex3f( (GLfloat)0 , (GLfloat)Ly , 0.001f);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly , 0.001f);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly , 0.001f);
	glVertex3f( (GLfloat)Lx, (GLfloat)0  , 0.001f);
	glVertex3f( (GLfloat)Lx, (GLfloat)0  , 0.001f);
	glVertex3f( (GLfloat)0 , (GLfloat)0  , 0.001f);
	glEnd();		

}

template<typename T>
void MacGrid2D<T>::glDrawUniformGrid( const RGBA &color, GLfloat width ) const	{

	glLineWidth(width);
	glColor4f( color[0], color[1], color[2], color[3] );

	for(index_type i=0;i<nX+1;++i)	{
		glBegin( GL_LINES );
		glVertex3f( (GLfloat)(i*h), (GLfloat)0  , 0.0f);
		glVertex3f( (GLfloat)(i*h), (GLfloat)Ly , 0.0f);
		glEnd();		
	}
	for(index_type i=0;i<nY+1;++i)	{
		glBegin( GL_LINES );
		glVertex3f( (GLfloat)0 , (GLfloat)(i*h) , 0.0f);
		glVertex3f( (GLfloat)Lx, (GLfloat)(i*h) , 0.0f);
		glEnd();		
	}
}

template<typename T>
void MacGrid2D<T>::glDrawNode( index_type i, index_type j, const RGBA &color, GLfloat size ) const {

	glPointSize(size);
	glColor4f( color[0], color[1], color[2], color[3] );

	const Vector2<T> pos = posNode(i,j);
	glBegin( GL_POINTS );
	glVertex3f( (GLfloat)pos[0], (GLfloat)pos[1], 0.0f );
	glEnd();
}

template<typename T>
void MacGrid2D<T>::glDrawNode( index_type iNode, const RGBA &color, GLfloat size ) const {
	glDrawNode(iNode / (nY+1), iNode % (nY+1), color, size);
}

template<typename T>
void MacGrid2D<T>::glDrawEdge( index_type iEdge, const RGBA &color, GLfloat width) const	{

	glLineWidth(width);
	glColor4f( color[0], color[1], color[2], color[3] );

	glBegin(GL_LINES);
	for(index_type i=0;i<2;i++)	{
		index_type node = idxNodeEdge(iEdge, i);
		Vector2<T> pos = posNode(node / (nY+1),node % (nY+1));
		glVertex3f( (GLfloat)pos[0], (GLfloat)pos[1], 0.0f );
	}
	glEnd();
}

template<typename T>
void MacGrid2D<T>::glDrawCell(index_type i, index_type j, const RGBA &color) const	{

	index_type iNode[4];
	Vector2<T> pos[4];
	for(index_type k=0;k<4;++k)	{
		iNode[k] = idxNodeCell(i,j,k);
		pos[k] = posNode(iNode[k]);
	}

	glColor4f(color[0], color[1], color[2], color[3]);	
	glBegin(GL_QUADS);
		glVertex3f(pos[0][0], pos[0][1], 0.0f);
		glVertex3f(pos[1][0], pos[1][1], 0.0f);
		glVertex3f(pos[2][0], pos[2][1], 0.0f);
		glVertex3f(pos[3][0], pos[3][1], 0.0f);
	glEnd();
}

template<typename T>
void MacGrid2D<T>::glDrawCell( index_type cell, const RGBA &color) const	{
	glDrawCell(cell / (nY+2), cell % (nY+2), color);
}

};		// geometryspace