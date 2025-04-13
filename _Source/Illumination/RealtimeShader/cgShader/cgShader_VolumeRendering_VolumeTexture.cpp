#include <stdafx.h>
#include <Illumination/RealtimeShader/cgShader/cgShader_VolumeRendering.h>

using namespace std;
using namespace basis;
using namespace arithmetic;
using namespace geometryspace;
using namespace shaderspace;

void cgShader_VolumeRendering::createVolumeTexture(GLsizei width, GLsizei height, GLsizei depth, GLubyte* data) {
	if(volume_texture != -1)
		glDeleteTextures(1, &volume_texture);

	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glGenTextures(1, &volume_texture);
	glBindTexture(GL_TEXTURE_3D, volume_texture);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, width, height, depth, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

	cout << "volume texture created _ " << volume_texture << endl;
}

void cgShader_VolumeRendering::setVolumeData_Default() {
#define VOLUME_TEX_SIZE 128

	int size = VOLUME_TEX_SIZE*VOLUME_TEX_SIZE*VOLUME_TEX_SIZE* 4;
	GLubyte *data = new GLubyte[size];

	for(int x = 0; x < VOLUME_TEX_SIZE; x++)
	{for(int y = 0; y < VOLUME_TEX_SIZE; y++)
	{for(int z = 0; z < VOLUME_TEX_SIZE; z++)
	{
		data[(x*4)   + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = z%250;
		data[(x*4)+1 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = y%250;
		data[(x*4)+2 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
		data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 230;

		Vector3<float> p =	Vector3<float>(x,y,z)- Vector3<float>(VOLUME_TEX_SIZE-20,VOLUME_TEX_SIZE-30,VOLUME_TEX_SIZE-30);
		bool test = (p.mag() < 42);
		if(test)
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 0;

		p =	Vector3<float>(x,y,z)- Vector3<float>(VOLUME_TEX_SIZE/2,VOLUME_TEX_SIZE/2,VOLUME_TEX_SIZE/2);
		test = (p.mag() < 24);
		if(test)
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 0;


		if(x > 20 && x < 40 && y > 0 && y < VOLUME_TEX_SIZE && z > 10 &&  z < 50)
		{

			data[(x*4)   + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 100;
			data[(x*4)+1 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
			data[(x*4)+2 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = y%100;
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
		}

		if(x > 50 && x < 70 && y > 0 && y < VOLUME_TEX_SIZE && z > 10 &&  z < 50)
		{

			data[(x*4)   + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
			data[(x*4)+1 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
			data[(x*4)+2 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = y%100;
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
		}

		if(x > 80 && x < 100 && y > 0 && y < VOLUME_TEX_SIZE && z > 10 &&  z < 50)
		{

			data[(x*4)   + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
			data[(x*4)+1 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 70;
			data[(x*4)+2 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = y%100;
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 250;
		}

		p =	Vector3<float>(x,y,z)- Vector3<float>(24,24,24);
		test = (p.mag() < 40);
		if(test)
			data[(x*4)+3 + (y * VOLUME_TEX_SIZE * 4) + (z * VOLUME_TEX_SIZE * VOLUME_TEX_SIZE * 4)] = 0;


	}}}

	createVolumeTexture(VOLUME_TEX_SIZE, VOLUME_TEX_SIZE, VOLUME_TEX_SIZE, data);

	delete [] data;	
}