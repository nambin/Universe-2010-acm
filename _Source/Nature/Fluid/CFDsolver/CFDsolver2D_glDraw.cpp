#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDsolver2D.h>

using namespace std;
using namespace fluidspace;

void CFDsolver2D::glDrawVelocity( real scale, const RGBA &color, const RGBA &colorP, real pointsize ) const {

	glColor4f( color[0], color[1], color[2], 0.8f );

	for(index_type i=1;i<=nX;++i)	{
		for(index_type j=1;j<=nY;++j)	{
			Vector2<real> vel, pos;
			vel[0] = 0.5f * (u[iS2U(i,j,0)] + u[iS2U(i,j,1)]);
			vel[1] = 0.5f * (v[iS2V(i,j,0)] + v[iS2V(i,j,1)]);

			pos = posS(i,j);
			glBegin(GL_LINES);
				glVertex3f(pos[0], pos[1], 0.001f);
				glVertex3f(pos[0]+vel[0]*scale, pos[1]+vel[1]*scale, 0.001f);
			glEnd();
		}
	}

	glPointSize(pointsize);
	glColor4f( colorP[0], colorP[1], colorP[2], 0.8f );

	glBegin(GL_POINTS);
	for(index_type i=1;i<=nX;++i)	{
		for(index_type j=1;j<=nY;++j)	{
			Vector2<real> pos = posS(i,j);
			glVertex3f(pos[0], pos[1], 0.001f);
		}
	}
	glEnd();
}

void CFDsolver2D::glDrawPressure( const RGBA &color ) const {

	real maxPressure = -LARGE, minPressure = LARGE;	
	for(index_type i=1;i<pressure.nX-1;++i)	{
		for(index_type j=1;j<pressure.nY-1;++j)	{			
			maxPressure = (maxPressure > -pressure(i,j)) ? maxPressure : -pressure(i,j);
			minPressure = (minPressure < -pressure(i,j)) ? minPressure : -pressure(i,j);
		}
	}

	for(index_type i=1;i<pressure.nX-1;++i)	{
		for(index_type j=1;j<pressure.nY-1;++j)	{	
			real ratio = (-pressure(i,j) - minPressure) / (maxPressure - minPressure);

			glColor4f( color[0] * ratio, color[1] * ratio, color[2] * ratio, color[3] );

			Vector2<real> pos[4];
			pos[0] = posNode(i-1,j-1);
			pos[1] = posNode(i,j-1);
			pos[2] = posNode(i,j);
			pos[3] = posNode(i-1,j);

			glBegin(GL_QUADS);
				glVertex3f(pos[0][0], pos[0][1], 0.001f);
				glVertex3f(pos[1][0], pos[1][1], 0.001f);
				glVertex3f(pos[2][0], pos[2][1], 0.001f);
				glVertex3f(pos[3][0], pos[3][1], 0.001f);
			glEnd();
		}
	}

	//cout << "Max Pressure : " << maxPressure << endl;
	//cout << "Min Pressure : " << minPressure << endl;
	//cout << "(1,1) Pressure : " << - pressure(1,1) << endl;
	//cout << "Center Pressure : " << - pressure(pressure.nX/2, pressure.nY/2) << endl;

	//real _sum = 0; index_type num = 0;
	//for(index_type i=1;i<pressure.nX-1;++i) {
	//	for(index_type j=1;j<pressure.nY-1;++j)	{	
	//		if(pressure(i,j) < 0) {
	//			_sum += pressure(i,j); num++;
	//		}
	//	}
	//}
	//cout << "Average Pressure : " << - _sum / num << endl;
}

void CFDsolver2D::glDrawKinematicBoundary( const RGBA &color ) const {
	for(index_type i=0;i<isKinematicBoundary.nX;++i)	{
		for(index_type j=0;j<isKinematicBoundary.nY;++j)	{				

			glColor4f( color[0], color[1], color[2], color[3] );

			if(isKinematicBoundary(i,j) == false )
				continue;

			Vector2<real> pos[4];
			pos[0] = posNode(i-1,j-1);
			pos[1] = posNode(i,j-1);
			pos[2] = posNode(i,j);
			pos[3] = posNode(i-1,j);

			glBegin(GL_QUADS);
				glVertex3f(pos[0][0], pos[0][1], 0.001f);
				glVertex3f(pos[1][0], pos[1][1], 0.001f);
				glVertex3f(pos[2][0], pos[2][1], 0.001f);
				glVertex3f(pos[3][0], pos[3][1], 0.001f);
			glEnd();
		}
	}
}