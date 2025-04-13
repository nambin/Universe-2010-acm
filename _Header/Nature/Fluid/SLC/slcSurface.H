//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// SLC Surface Tracker
//   -- Header File
//  
// Primary Author: Adam Bargteil (adamb@cs.berkeley.edu)
// 
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Copyright (c) 2003-2005, Regents of the University of California.  All
// rights reserved.
//
// This software is part of the Berkeley Fluid Animation & Simulation
// Toolkit.  The name "Berkeley Fluid Animation & Simulation Toolkit" is
// a trademark of the Regents of the University of California.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//   Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//
//  Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//
//  Redistributions in binary form as an executable program, or part of
//  an executable program must, when run for the first time by a given
//  user, prominently display the above copyright notice, this list of
//  conditions and the following disclaimer.
//
//  Neither the name of the University of California, Berkeley nor the
//  names of its contributors may be used to endorse or promote products
//  derived from this software without specific prior written
//  permission.
//
//  ** Animations, still images, or other works created using this
//  ** software must clearly indicate in the list of credits that this
//  ** software was used.  The software shall be referred to as the
//  ** Berkeley Fluid Animation & Simulation Toolkit.  If the software
//  ** is included, either in whole or in part, within another software
//  ** package, credit must still be given to the Berkeley Fluid
//  ** Animation & Simulation Toolkit.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

#ifndef slcSurface_H
#define slcSurface_H

#include <fstream>
#include <iostream>
//
//#include <Nature/Fluid/SLC/smVector.H>
//#include <Nature/Fluid/SLC/abUtil.H>
//#include <Nature/Fluid/SLC/distTree.H>

#include "smVector.H"
#include "abUtil.H"
#include "distTree.H"

namespace slcspace {

#pragma warning(disable:4101)
#pragma warning(disable:4018)

class DtFMMHeap;

class __slc_hash_unsigned_long_long {
public:
	inline std::size_t operator()(const unsigned long long __s) const {
		unsigned long __h = (__s >> 32) + (__s & 0xFFFFFFFF);
		return std::size_t(__h);
	};

public:
};

class SlcSurface {
public:
	SlcSurface();
	~SlcSurface();

	DtTree tree;
	AbTriList triangles; // triangles in polygon mesh
	AbVecList meshPts; // points in polygon mesh
	AbVecList normals; // normals (at points) in polygon mesh
	AbVecList psuedoNormals; // angle weighted normals (at points) in polygon mesh
	AbVecList faceNormals; // face normals in polygon mesh

	// whenever phi is evaluated, the returned result will be min(minFunc, phi)
	double (*minFunc)(const double *x);

	// whenever phi is evaluated, the returned result will be max(maxFunc, phi)
	double (*maxFunc)(const double *x);

	bool builtFromTriangles; // was the tree constructed from triangles or a phi func.

#ifdef _WIN32
#if _MSC_VER <= 1200
	std::set<int> visitedTris;
  std::map<unsigned long long,unsigned int, __slc_hash_unsigned_long_long> edgeList;
#else
	stdext::hash_set<int> visitedTris;
  //stdext::hash_map<unsigned long long, unsigned int, __slc_hash_unsigned_long_long> edgeList;
	stdext::hash_map<unsigned long long, unsigned int> edgeList;
#endif
#else
	__gnu_cxx::hash_set<int> visitedTris;
  __gnu_cxx::hash_map<unsigned long long, unsigned int, 
											__slc_hash_unsigned_long_long> edgeList;
#endif

	// Marching cubes functions.
	// vertid creates new vertices or finds already created ones.
	// dotet creates triangles for a tet through which the surface passes.
	// contourTree drives the contouring.

	int vertid(double (*evalPhi)(const SmVector3 &x), int a, int b, SmVector3 x, SmVector3 y, bool apos, bool bpos, double aphi, double bphi);
	void dotet(double (*evalPhi)(const SmVector3 &x), DtCell *cell, int a, int b, int c, int d, SmVector3 w, SmVector3 x, SmVector3 y, SmVector3 z);
	void contourTree(double (*evalPhi)(const SmVector3 &x));

	template<typename FUNC>
	int vertid(const FUNC&, int a, int b, SmVector3 x, SmVector3 y, bool apos, bool bpos, double aphi, double bphi);

	template<typename FUNC>
	void dotet(const FUNC&, DtCell *cell, int a, int b, int c, int d, SmVector3 w, SmVector3 x, SmVector3 y, SmVector3 z);

	template<typename FUNC>
	void contourTree(const FUNC&);

	// end marching cubes functions

	void redistance(); // redistance the tree.

	void computeNormals();
	void computePsuedoNormals();

	// a variety of signed-distance evaluation routines.  x is the
	// evaluation point.  If interpolate is true, the evaluation is
	// gotten from linear interpolation of the values stored at
	// octree vertices.  s & t are barycentric coordinates and face is the 
	// triangle of the nearest point in the mesh.  creval fits a cubic
	// catmull-rom spline and evaluates that.

	double eval(const double *x, bool interpolate = false);
	double eval(const double *x, double &s, double &t, int &face, bool interpolate = false);
	double eval(const SmVector3 &x, bool interpolate = false);
	double eval(const SmVector3 &x, double &s, double &t, int &face, bool interpolate = false);
	double creval(const SmVector3 &x);
	double creval(const double *x);


	// evaluate the gradient at x.
	SmVector3 grad(const SmVector3 &x);
	SmVector3 grad(const double *x);

	// distance from a point to a triangle.  This searches the entire triangle
	// list of cell.  The second version is more
	// complicated because it is meant to be called for each cell in the concentric
	// triple and has to deal with sharp features, etc.
	double pointToTri(DtCell *c, const SmVector3 &x, double &phi);
	double pointToTri(DtCell *cell, const SmVector3 &x, double &phi, int &sgn,
											double &s, double &t, int &face, double &dp, int edge[2]);

	// i/o functions
	void read(std::ifstream &in);
	void write(std::ofstream &out);
	void objDump(char *fname);
	void ribDump(char *fname);
	void bvtDump(char *fname);
	void bvtRead(char *fname);
	void bdtDump(char *fname);
	void bdtRead(char *fname);

#ifdef INVENTOR
	SoSeparator *ivTree();
	void ivDump(char *fname);
#endif

private:	
	void addToFMMHeap(int i, int j, int k, DtFMMHeap &heap,
										short *levels, bool *accepted, AbDoubleList &newPhi, int foo);
	double exacteval(DtCell *c, const SmVector3 &x, double &s, double &t, 
									 int &face);
	double interpolate(DtCell *c, const SmVector3 &x);
};

void readObjFile(char *fname, AbVecList &pts, AbTriList &triangles);

#ifdef INVENTOR
#include <Inventor/SoOutput.h>
#include <Inventor/actions/SoWriteAction.h>
int SmSoWriteIVFile(const char *filename, SoNode *node);
#endif

// trace back in time (semi-Lagrangian)
inline SmVector3 traceBack(const SmVector3 &x, SmVector3 (*evalVel)(const SmVector3 &x),
													 double dt) {
	SmVector3 v = evalVel(x);
	SmVector3 y = x-0.5*dt*v;
	v = evalVel(y);
	y = x-dt*v;
	return y;
}

// trace forward in time (fully Lagrangian)
inline SmVector3 traceForward(const SmVector3 &x, SmVector3 (*evalVel)(const SmVector3 &x),
															double dt) {
	SmVector3 v = evalVel(x);
	SmVector3 y = x+0.5*dt*v;
	v = evalVel(y);
	y = x+dt*v;
	return y;
}

// trace back in time (semi-Lagrangian)
template<typename VEL>
inline SmVector3 traceBack(const SmVector3 &x, const VEL& evalVel, double dt) {
	SmVector3 v = evalVel(x);
	SmVector3 y = x-0.5*dt*v;
	v = evalVel(y);
	y = x-dt*v;
	return y;
}

// trace forward in time (fully Lagrangian)
template<typename VEL>
inline SmVector3 traceForward(const SmVector3 &x, const VEL& evalVel, double dt) {
	SmVector3 v = evalVel(x);
	SmVector3 y = x+0.5*dt*v;
	v = evalVel(y);
	y = x+dt*v;
	return y;
}



// Fast marching method code.
class DtFMMHeapEntry {
 public:
  DtFMMHeapEntry(){};
  DtFMMHeapEntry(int x, int y, int z, double phi) {i = x, j = y, k = z, this->phi = phi;};
  DtFMMHeapEntry &operator=(const DtFMMHeapEntry &e){i=e.i,j=e.j,k=e.k,phi=e.phi; return *this;};
  int i,j,k;
  double phi;
};
  
class DtFMMHeap {
 public:
  DtFMMHeap(int dim, int dim2);
  ~DtFMMHeap();
  void addToHeap(int i, int j, int k, double phi);
  DtFMMHeapEntry removeFromHeap();
  inline unsigned int numEntries() const {return _numEntries;};
  void setNumEntries(unsigned int num, int copyOld = 0, int exactSize = 1);
  DtFMMHeapEntry *findEntry(int i, int j, int k);
 private:
  int dim, dim2;
  void swap(int i, int j);
  int _numEntries;
  DtFMMHeapEntry *_entries;
  int _numMalloc;

#ifdef _WIN32
#if _MSC_VER <= 1200
  std::map<int,int> mymap;
#else
  stdext::hash_map<int,int> mymap;
#endif
#else
  __gnu_cxx::hash_map<int, int> mymap;
#endif
};

////////////////////////////////////////////////////////////////


template<typename FUNC>
void SlcSurface::contourTree(const FUNC& evalPhi) {
	int i,j,k;
	SmVector3 x0, x1, x2, x3, x4, x5, x6, x7;
	bool apos, bpos, cpos, dpos;	
	AbIndexList cubes;
	DtCell *c;
	edgeList.clear();

	// create list of cubes to look at
	c = &(tree.cells[0]);
	for (i=0; i<tree.cells.numCells(); i++, c++) {
		if (c->level != tree.max_level) continue;
		int foo = sign(tree.phi[c->vertices[0]]);
		if ((foo != sign(tree.phi[c->vertices[1]])) ||
			(foo != sign(tree.phi[c->vertices[2]])) ||
			(foo != sign(tree.phi[c->vertices[3]])) ||
			(foo != sign(tree.phi[c->vertices[4]])) ||
			(foo != sign(tree.phi[c->vertices[5]])) ||
			(foo != sign(tree.phi[c->vertices[6]])) ||
			(foo != sign(tree.phi[c->vertices[7]])))
			cubes.add(i);
	}

	for (i=0; i<cubes.getNum(); i++) {
		c = &(tree.cells[cubes[i]]);
		x0[0] = c->lc[0]; x0[1] = c->lc[1]; x0[2] = c->lc[2];
		x1[0] = c->lc[0]; x1[1] = c->lc[1]; x1[2] = c->uc[2];
		x2[0] = c->lc[0]; x2[1] = c->uc[1]; x2[2] = c->lc[2];
		x3[0] = c->lc[0]; x3[1] = c->uc[1]; x3[2] = c->uc[2];
		x4[0] = c->uc[0]; x4[1] = c->lc[1]; x4[2] = c->lc[2];
		x5[0] = c->uc[0]; x5[1] = c->lc[1]; x5[2] = c->uc[2];
		x6[0] = c->uc[0]; x6[1] = c->uc[1]; x6[2] = c->lc[2];
		x7[0] = c->uc[0]; x7[1] = c->uc[1]; x7[2] = c->uc[2];

		dotet(evalPhi,c,c->vertices[0],c->vertices[2],c->vertices[4],c->vertices[1],
			x0,x2,x4,x1);
		dotet(evalPhi,c,c->vertices[6],c->vertices[2],c->vertices[1],c->vertices[4],
			x6,x2,x1,x4);
		dotet(evalPhi,c,c->vertices[6],c->vertices[2],c->vertices[3],c->vertices[1],
			x6,x2,x3,x1);
		dotet(evalPhi,c,c->vertices[6],c->vertices[4],c->vertices[1],c->vertices[5],
			x6,x4,x1,x5);
		dotet(evalPhi,c,c->vertices[6],c->vertices[1],c->vertices[3],c->vertices[5],
			x6,x1,x3,x5);
		dotet(evalPhi,c,c->vertices[6],c->vertices[3],c->vertices[7],c->vertices[5],
			x6,x3,x7,x5);
	}

	edgeList.clear();

	computeNormals();
	computePsuedoNormals();
}

template<typename FUNC>
int SlcSurface::vertid(const FUNC& evalPhi, int a, int b, SmVector3 x, 
					   SmVector3 y, bool apos, bool bpos, double aphi, double bphi) {
						   int i;
						   double TOL = 1e-10;
						   double MAXITER = 100;
						   double p,sum;
						   bool cpos;
						   int iter = 0;
						   unsigned long long ei;
						   int vi;
						   SmVector3 z,v,zprime;

						   if (a > b) { 
							   int ti;
							   bool tb;
							   double td;
							   SmVector3 tv;
							   ti = a; a = b; b = ti; 
							   tv = x; x = y; y = tv;
							   tb = apos; apos = bpos; bpos = tb;
							   td = aphi; aphi = bphi; bphi = td;
						   }

						   ei = a;
						   ei = ei << 32;
						   ei+=b;
						   if (edgeList.count(ei))
							   return edgeList[ei];

						   sum = fabs(aphi-bphi);
						   z = (fabs(bphi)/sum)*x+(fabs(aphi)/sum)*y;
						   while (fabs(p=evalPhi(z)) > TOL && iter++ < MAXITER) {
							   if ((cpos = (p > 0.0)) == apos) {
								   x = z;
								   aphi = p;
							   } else {
								   y = z;
								   bphi = p;
							   }
							   sum = fabs(aphi-bphi);
							   z = (fabs(bphi)/sum)*x+(fabs(aphi)/sum)*y;
						   }

						   vi = meshPts.numVecs();
						   edgeList[ei] = vi;
						   meshPts.setNumVecs(meshPts.numVecs()+1,1,0);
						   meshPts[vi] = z;

						   return vi;
}

template<typename FUNC>
void SlcSurface::dotet(const FUNC& evalPhi, DtCell *cell, 
					   int a, int b, int c, int d,
					   SmVector3 w, SmVector3 x, SmVector3 y, SmVector3 z) {
						   bool apos, bpos, cpos, dpos;
						   double aphi, bphi, cphi, dphi;
						   int index=0, e1, e2, e3, e4, e5, e6;
						   SmVector3 foo;
						   if (apos = ((aphi=tree.phi[a]) > 0.0)) index += 8;
						   if (bpos = ((bphi=tree.phi[b]) > 0.0)) index += 4;
						   if (cpos = ((cphi=tree.phi[c]) > 0.0)) index += 2;
						   if (dpos = ((dphi=tree.phi[d]) > 0.0)) index += 1;

						   if (apos != bpos) e1 = vertid(evalPhi,a,b,w,x,apos,bpos,aphi,bphi);
						   if (apos != cpos) e2 = vertid(evalPhi,a,c,w,y,apos,cpos,aphi,cphi);
						   if (apos != dpos) e3 = vertid(evalPhi,a,d,w,z,apos,dpos,aphi,dphi);
						   if (bpos != cpos) e4 = vertid(evalPhi,b,c,x,y,bpos,cpos,bphi,cphi);
						   if (bpos != dpos) e5 = vertid(evalPhi,b,d,x,z,bpos,dpos,bphi,dphi);
						   if (cpos != dpos) e6 = vertid(evalPhi,c,d,y,z,cpos,dpos,cphi,dphi);

						   switch (index) {
	case 1: cell->triangles.add(triangles.add(e5,e3,e6)); break;
	case 2: cell->triangles.add(triangles.add(e2,e4,e6)); break;
	case 3: cell->triangles.add(triangles.add(e3,e4,e5)); 
		cell->triangles.add(triangles.add(e3,e2,e4)); break;
	case 4: cell->triangles.add(triangles.add(e1,e5,e4)); break;
	case 5: cell->triangles.add(triangles.add(e3,e4,e1)); 
		cell->triangles.add(triangles.add(e3,e6,e4)); break;
	case 6: cell->triangles.add(triangles.add(e1,e6,e2)); 
		cell->triangles.add(triangles.add(e1,e5,e6)); break;
	case 7: cell->triangles.add(triangles.add(e1,e3,e2)); break;
	case 8: cell->triangles.add(triangles.add(e1,e2,e3)); break;
	case 9: cell->triangles.add(triangles.add(e1,e6,e5)); 
		cell->triangles.add(triangles.add(e1,e2,e6)); break;
	case 10: cell->triangles.add(triangles.add(e1,e6,e3)); 
		cell->triangles.add(triangles.add(e1,e4,e6)); break;
	case 11: cell->triangles.add(triangles.add(e1,e4,e5)); break;
	case 12: cell->triangles.add(triangles.add(e3,e4,e2)); 
		cell->triangles.add(triangles.add(e3,e5,e4)); break;
	case 13: cell->triangles.add(triangles.add(e6,e4,e2)); break;
	case 14: cell->triangles.add(triangles.add(e5,e6,e3)); break;
						   }
}




//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: slcSurface.H,v $
// Revision 1.2  2005/06/10 20:08:47  adamb
// Various bug fixes, added iuc and ilc to cell output, rebuild vertMap after readining in a surface, fixed argc checking in the drivers.
//
// Revision 1.1.1.1  2005/05/25 05:41:04  adamb
// Initial Revision
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

};

#endif