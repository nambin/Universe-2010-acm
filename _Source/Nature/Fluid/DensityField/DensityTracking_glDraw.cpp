#include <stdafx.h>
#include <Nature/Fluid/DensityTracking/DensityTracking.h>
#include <Illumination/RealtimeShader/ShaderViewer/ShaderViewer.h>

using namespace std;
using namespace fluidspace;

void DensityTracking2D::glDraw(const RGBA& c) const {
	glDraw(density, c);
}

void DensityTracking2D::glDraw(const Grid2r& density, const RGBA& color) {
	Vector2<real> posC;
	Vector2<real> pos[4];
	for(Grid2r::index_type i=1;i<density.nX-1;++i)	{
		for(Grid2r::index_type j=1;j<density.nY-1;++j)	{
			posC = density.pos(i,j);
			pos[0] = posC + Vector2<real>(-0.5f * density.dx, -0.5f * density.dy);
			pos[1] = posC + Vector2<real>(+0.5f * density.dx, -0.5f * density.dy);
			pos[2] = posC + Vector2<real>(+0.5f * density.dx, +0.5f * density.dy);
			pos[3] = posC + Vector2<real>(-0.5f * density.dx, +0.5f * density.dy);

			real frac = density(i,j);

			if( frac > 1.0f )	glColor4f( color[0], color[1], color[2], color[3] );
			else				glColor4f( color[0]*frac, color[1]*frac, color[2]*frac, color[3] );

			glBegin(GL_QUADS);
				glVertex3f(pos[0][0], pos[0][1], 0.0f);
				glVertex3f(pos[1][0], pos[1][1], 0.0f);
				glVertex3f(pos[2][0], pos[2][1], 0.0f);
				glVertex3f(pos[3][0], pos[3][1], 0.0f);
			glEnd();
		}
	}
}

void DensityTracking3D::glDraw(const RGBA& c) const {
	if(!shader)
		throwError("Error in DensityTracking3D::glDraw() _ Shader should be initialized first");
		
	if(isUpdated) {
		shader->setVolumeData(density, c);
		isUpdated = false;
	}

	shader->glDraw();
}

void DensityTracking3D::initialize_Shader(unsigned int width, unsigned int height) {
	if(!shader)
		shader = new cgShader_VolumeRendering();
	shader->resizeScreen(width, height);
}