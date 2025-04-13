#include <stdafx.h>
#include "Header/TestFunction.h"

using namespace std;
using namespace basis;
using namespace arithmetic;
using namespace geometryspace;
using namespace illuminationspace;

//RandomClass::rand;

void testFunction() {

	//testFunction_XML();
	//testFunction_Geometry();
	//basis::TestTBB();

//	SVDTest();
//	VectorMatrixPCATest();
//	testPCA();
//	Vector3<real> p0, p1, p2, bmin, bmax;

//	IntersectionTestLibrary::testTriangleAABB(p0, p1, p2, bmin, bmax);
//	DistanceLibrary::SqDistPointAABB(p0, bmin, bmax);

	//{
	//	int val,n;

	//	istringstream iss;
	//	string strvalues = "32 240 2 1450";

	//	iss.str (strvalues);

	//	for (n=0; n<4; n++)	{
	//		iss >> val;
	//		cout << val+1 << endl;
	//	}
	//}

//	int iii = 0;
//	int kkk0 = 1;
//	char str[256];
//	
//
//#define VARIABLE(_name, _i) (_name##_i)
//#define VARIABLE01(_name, _i) (_name)
//#define VARIABLE02(_name, _i) (#_i)
//
//	//cout << VARIABLE("kkk", iii) << endl;
//	cout << VARIABLE01("kkk") << endl;
//	cout << VARIABLE02("kkk", iii) << endl;

	//basis::Timer timer;

	//DataPath datapath;
	//datapath.putPath("../Fluid");
	//datapath.putPath("../Mesh/OBJs");
	//cout << datapath.getPath("body.obj") << endl;
	//cout << datapath.getFile("OBJs/body.obj") << endl;

	//timer.printTime("");

}

void testFunction_Geometry() {

	illuminationspace::ONB onb;
	cout << onb;

	geometryspace::Transform<double> t[10];
	Matrix4<double> mat;

	t[0] = t[1] = t[2];
	t[6] = t[5].translate(Vector3<double>(0.0f));
	//t[0].set(mat);
	t[0].set(mat);

	for(int i=0;i<10;++i) {
		if(i % 5 == 0) 
			cout << endl;
		cout << t[i].getTransformValue()->getRefCount() << " ";		
	}
}

#ifdef _USE_XML_
void testFunction_XML() {

	//string fname("color01.xml");
	//string ap("../Color/");
	//string ap("../Color/Data/Color/");
	//string path = MakeAbsolutePath(ap, fname);
	//cout << path << endl;

	RGBA color0(1.0f, 0.0f, 0.5f, 1.0f), color1;
	XmlLibrary::exportXmlFile(&color0, "../Color/Data/Color/color01.xml");
	XmlLibrary::importXmlFile(&color1, "../Color/Data/Color/color01.xml");
	cout << color0 << endl << color1 << endl;

	virtualspace::Material mat0, mat1;
	mat0.setAmbient(RGBA::red());
	mat0.setDiffuse(RGBA::green());
	mat0.setSpecular(RGBA::blue());
	mat0.setShininess(5.0f);
	XmlLibrary::exportXmlFile(&mat0, "../Color/Data/Material/material01.xml");
	XmlLibrary::importXmlFile(&mat1, "../Color/Data/Material/material01.xml");

	cout << mat0 << mat1 << endl;
}
#endif