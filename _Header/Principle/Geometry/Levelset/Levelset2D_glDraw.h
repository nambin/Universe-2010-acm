#pragma once

namespace geometryspace {

template<typename T> 
void Levelset2D<T>::glDrawLiquid(const RGBA &color) const {
	glColor4f(color[0], color[1], color[2], color[3]);

	glBegin(GL_TRIANGLES);
	for(size_t i=0;i<renderingTriangles.size();i++)	{
		glVertex3f(renderingTriangles[i].v[0][0], renderingTriangles[i].v[0][1], 0);
		glVertex3f(renderingTriangles[i].v[1][0], renderingTriangles[i].v[1][1], 0);
		glVertex3f(renderingTriangles[i].v[2][0], renderingTriangles[i].v[2][1], 0);
	}
	glEnd();
}

template<typename T> 
void Levelset2D<T>::glDrawInterfaceEdges(const RGBA &color, GLfloat width) const {
	glColor4f(color[0], color[1], color[2], color[3]);
	glLineWidth(width);

	glBegin(GL_LINES);
	for(std::size_t i=0;i<renderingSegments.size();++i) {
		glVertex3f(renderingSegments[i].v[0][0], renderingSegments[i].v[0][1], 0);
		glVertex3f(renderingSegments[i].v[1][0], renderingSegments[i].v[1][1], 0);
	}

	glEnd();
}

template<typename T> 
void Levelset2D<T>::glDrawPhi(const RGBA &colorPos, const RGBA &colorNeg, T band, GLfloat pointsize) const {

	glPointSize(pointsize);

	glBegin(GL_POINTS);
	for(index_type i=1;i<phi.nX-1;i++)	{
		for(index_type j=1;j<phi.nY-1;j++)	{

			T scale = (band - abs(phi(i,j)) * hInv) / band;

			if(scale < 0)			continue;
			else if(phi(i,j) > 0)	glColor4f(colorPos[0] * scale, colorPos[1] * scale, colorPos[2] * scale, colorPos[3]);
			else					glColor4f(colorNeg[0] * scale, colorNeg[1] * scale, colorNeg[2] * scale, colorNeg[3]);

			Vector2<T> pos = phi.pos(i,j);
			glVertex3f(pos[0], pos[1], 0.001f);

		}
	}	
	glEnd();
}

};		// geometryspace