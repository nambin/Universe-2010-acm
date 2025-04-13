#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/BiTree/QuadTreeViewer.h>
#include <Principle/Geometry/DistanceFieldLibrary/DistanceFields.h>

using namespace virtualspace;
using namespace geometryspace;

QuadTreeViewer::QuadTreeViewer() {
	ioHandler = new QuadTreeViewer_IOHandler(this);

	quadtree.init(8, Vector2<double>(0), Vector2<double>(1));
	quadtree.makeTree(*quadtree.getSplitProc());
	
	for(int i=0;i<100;++i) {
		pos[i] = Vector2<double>(rand(), rand());
		pos[i] /= RAND_MAX;
	}

	for(int i=0;i<50;++i) {
		double random = (rand() / (double)RAND_MAX) * 0.05f;
		int cellIdx = (int)(quadtree.cells.size() * random);
		pos[i*2] = quadtree.getPos(quadtree.cells[cellIdx]->locCode) + quadtree.getCellSize(quadtree.cells[cellIdx]->level) * pos[i*2];
		pos[i*2+1] = quadtree.getPos(quadtree.cells[cellIdx]->locCode) + quadtree.getCellSize(quadtree.cells[cellIdx]->level) * pos[i*2+1];

		Vector2<double> _min = minVec(pos[i*2], pos[i*2+1]);
		Vector2<double> _max = maxVec(pos[i*2], pos[i*2+1]);
//		pos[i*2  ] = _min;	pos[i*2+1] = _max;
	}

	idxCell = 0;
	drawTree = drawRenderingPrimitive = true;

	//quadtree.setPhi(CircleDistanceField2D<double>(Vector2<double>(0.5f,0.5f), 0.1997f));
	//quadtree.setPhi(RectangleDistanceField2D<double>(Vector2<double>(0.1f,0.2f), Vector2<double>(0.42f,0.81f)));
	quadtree.setPhi(ZalesakDistanceField2D<double>(1.0f));
}

void QuadTreeViewer::glBind() const {

	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	if(drawTree)				quadtree.glDrawQuadTree(RGBA::white());
	if(drawRenderingPrimitive)	quadtree.glDraw(RGBA::blue());
	
	//for(unsigned int i=0;i<100;i+=1) {
	//	glPointSize(5.0f);	
	//	glColor3f(1,0,0);
	//	glBegin(GL_POINTS);
	//		glVertex3f(pos[i][0], pos[i][1], 0.001f);
	//	glEnd();		

	//	QuadTreeLibrary::glDrawCell(quadtree.getLocateCell(pos[i]), &quadtree, RGBA::blue());
	//}

	//for(unsigned int i=2*idxCell;i<2*idxCell+2;++i) {
	//	glPointSize(5.0f);	
	//	glColor3f(1,0,0);
	//	glBegin(GL_POINTS);
	//		glVertex3f(pos[i][0], pos[i][1], 0.001f);
	//	glEnd();		
	//}
	//QuadTreeLibrary::glDrawCell(quadtree.getLocateRegion(pos[2*idxCell], pos[2*idxCell+1]), RGBA::blue(), &quadtree);

	QuadTreeLibrary::glDrawCell(quadtree.cells[idxCell], &quadtree, RGBA::red());
	QuadTreeLibrary::glDrawCellVertices(quadtree.cells[idxCell], &quadtree, RGBA::yellow(), 7.0f);
//	QuadTreeLibrary::glDrawCellNeighbors(quadtree.cells[idxCell], &quadtree, RGBA::red(), 1.0f);	
}

void QuadTreeViewer::idle() {
}
