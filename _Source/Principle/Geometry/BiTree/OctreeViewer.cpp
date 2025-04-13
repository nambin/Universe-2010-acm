#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/BiTree/OctreeViewer.h>
#include <Principle/Geometry/DistanceFieldLibrary/DistanceFields.h>

using namespace virtualspace;
using namespace geometryspace;

OctreeViewer::OctreeViewer() {
	ioHandler = new OctreeViewer_IOHandler(this);

	octree.init(8, Vector3<double>(0), Vector3<double>(1));
	octree.makeTree(*octree.getSplitProc());

	//for(int i=0;i<100;++i) {
	//	pos[i] = Vector3<double>(rand(), rand(), rand());
	//	pos[i] /= RAND_MAX;
	//}

	//for(int i=0;i<50;++i) {
	//	double random = (rand() / (double)RAND_MAX) * 0.05f;
	//	int cellIdx = (int)(octree.cells.size() * random);
	//	pos[i*2] = octree.getPos(octree.cells[cellIdx]->locCode) + octree.getCellSize(octree.cells[cellIdx]->level) * pos[i*2];
	//	pos[i*2+1] = octree.getPos(octree.cells[cellIdx]->locCode) + octree.getCellSize(octree.cells[cellIdx]->level) * pos[i*2+1];

	//	Vector2<double> _min = minVec(pos[i*2], pos[i*2+1]), _max = maxVec(pos[i*2], pos[i*2+1]);
	//	pos[i*2  ] = _min;		
	//	pos[i*2+1] = _max;
	//}

	idxCell = 0;
	drawTree = drawRenderingPrimitive = true;
	drawOnlyMinimumLeaf = false;

	octree.setPhi(SphereDistanceField3D<double>(Vector3<double>(0.3f,0.7f,0.45f), 0.117f));	
	//octree.setPhi(CubeDistanceField3D<double>(Vector3<double>(0.1f,0.2f,0.13f), Vector3<double>(0.42f,0.81f,0.74f)));
}

void OctreeViewer::glBind() const {

	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	if(drawTree) {
		if(drawOnlyMinimumLeaf)	octree.glDrawOnlyMinimumLeaf(RGBA::white());
		else					octree.glDrawOctree(RGBA::white());
	}
	if(drawRenderingPrimitive)	{
		state.enable_Lighting();
		state.enable_ColorMaterial();
		state.glBind();
		octree.glDraw(RGBA::blue());
	}

	//for(unsigned int i=0;i<40;i+=2) {
	//	glPointSize(5.0f);	
	//	glColor3f(1,0,0);
	//	glBegin(GL_POINTS);
	//		glVertex3f(pos[i*2][0], pos[i*2][1], 0.001f);
	//	glEnd();		

	//	glDrawCell(octree.getLocateCell(pos[i*2]), &octree, RGBA::blue());
	//}

	//for(unsigned int i=2*idxCell;i<2*idxCell+2;++i) {
	//	glPointSize(5.0f);	
	//	glColor3f(1,0,0);
	//	glBegin(GL_POINTS);
	//		glVertex3f(pos[i][0], pos[i][1], 0.001f);
	//	glEnd();		
	//}
	//glDrawCell(octree.getLocateRegion(pos[2*idxCell], pos[2*idxCell+1]), RGBA::blue(), &octree);

	OctreeLibrary::glDrawCell(octree.cells[idxCell], &octree, RGBA::red(), 7.0f);
	OctreeLibrary::glDrawCellVertices(octree.cells[idxCell], &octree, RGBA::blue(), 15.0f);
}

void OctreeViewer::idle() {
}