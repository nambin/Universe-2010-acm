//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Distance Tree
//   -- Header File
//  
// Primary Author: Adam Bargteil (adamb@cs.berkeley.edu)
//
// The distance tree is simply a balanced octree annotated with
// signed distance values.  This tree has two build functions,
// one takes a triangle mesh and splits the tree near the 
// triangles.  The other takes a function which returns a 
// signed distance at any query point.
//
// The individual octree cell's additionally contain a list
// of indices of triangles which intersect the cell.
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


#ifndef distTree_H
#define distTree_H

#include <fstream>
#include <iostream>
#include "smVector.H"
#include "abUtil.H"
//#include <Nature/Fluid/SLC/smVector.H>
//#include <Nature/Fluid/SLC/abUtil.H>

#ifdef _WIN32
#if _MSC_VER >= 1300
#include <hash_map>
#else
#include <map>
#endif
#else
#include <ext/hash_map>
#endif

namespace slcspace {

#pragma warning(disable:4101)
#pragma warning(disable:4018)

class DtTree;

////////////////////////////////////////////////////////////////

class DtCell {
public:
  int level;          // depth in the octree of this cell											// Octree level
  SmVector3 lc, uc;   // coordinates of the lower and upper corners of the cell						// vertex coordinate for locCode
	int ilc[3], iuc[3]; // integer coordinates of the lower and upper corner						// locCode
  int vertices[8];    // indices into the octree's pts array of the vertices of this cell			// vertices indices
  int parent;         // index of the parent of this cell											// parent index
  int child[8];       // indices of the children of this cell										// child index
  double dist;        // (approx) distance from the center of the cell to the surface
	AbIndexList triangles; // list of triangles intersecting this cell
  DtCell() {};
  DtCell &operator=(const DtCell &that);
  DtCell(const DtCell &that);
  bool contains(const double *x); // does x lie in this cell?
	int clamp(double *x);           // make x lie in this cell.

	void read(std::ifstream &in);
	void write(std::ofstream &out);
};

NATUREDLL std::istream &operator>>(std::istream &strm, DtCell &c);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const DtCell &c);

///////////////////////////////////////////////////////////////
class DtCellList {
 public:
  DtCellList();
  DtCellList(const DtCellList &that);
  ~DtCellList() { delete [] _cells;};

  DtCellList &operator=(const DtCellList &that);

  inline DtCell &index(int i) {return _cells[i];};
  inline const DtCell &index(int i) const {return _cells[i];};
  inline DtCell &operator[](int i) {return _cells[i];};
  inline const DtCell &operator[](int i) const {return _cells[i];};

  inline unsigned int numCells() const {return _numCells;};
  void setNumCells(unsigned int num, int copyOld = 0, int exactSize = 1);

	void read(std::ifstream &in);
	void write(std::ofstream &out);

 private:
  int _numCells;
  int _numMalloc;
  DtCell *_cells;
  AbIndexList invalid;
};

NATUREDLL std::istream &operator>>(std::istream &strm, DtCellList &p);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const DtCellList &p);
  
/////////////////////////////////////////////////////////////////

class DtTree {
public:
  double lc, uc; // upper and lower coordinates of the tree 
	               // (the domain is (lc,lc,lc)x(uc,uc,uc)
	int dim, dim2; // integer dimension and dimension squared
  int max_level; // how deep can the tree go?
  DtCellList cells; // list of cells
	AbVecList pts; // oct-tree vertices
  AbDoubleList phi; // phi values at oct-tree points
	
  DtTree() {max_level=9; vertMap.clear();};
  ~DtTree() {vertMap.clear();};

	// functions to build the tree.  One takes a phi function, the other
	// a list of triangles & vertices.
	template<typename FUNC> void buildTree(const FUNC& evalPhi, double input_lc, double input_uc, int max_level);

	void buildTree(double (*evalPhi)(const SmVector3 &x), 
								 double input_lc, double input_uc, int max_level);
	void buildTree(AbVecList &meshPts, AbTriList &triangles, AbVecList &faceNormals, 
								 double input_lc, double input_uc, int max_level);

  NATUREDLL friend std::istream &operator>>(std::istream &strm, DtTree &t);
  NATUREDLL friend std::ostream &operator<<(std::ostream &strm, const DtTree &t);

	// force x into the domain
	int clamp(double *x);
	int clamp(SmVector3 &x);

	// find the cell containing x (maybe cell at level level).
	DtCell* findCell(const SmVector3 &x);
	DtCell* findCell(const SmVector3 &x, int &level);

	void read(std::ifstream &in);
	void write(std::ofstream &out);

#ifdef _WIN32
#if _MSC_VER <= 1200
  std::map<int,int> vertMap;
#else
  stdext::hash_map<int,int> vertMap;
#endif
#else
  __gnu_cxx::hash_map<int, int> vertMap;
#endif

private:
	int nvertices;
};

NATUREDLL std::istream &operator>>(std::istream &strm, DtTree &t);
NATUREDLL std::ostream &operator<<(std::ostream &strm,const DtTree &t);

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: distTree.H,v $
// Revision 1.1.1.1  2005/05/25 05:41:01  adamb
// Initial Revision
//
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

template<typename FUNC>
void DtTree::buildTree(const FUNC& evalPhi, double input_lc, double input_uc, int max_level) {
	const double width = 3;

   int i,j,l,n,p,vert[19],ilc[3],iuc[3],ih,c=0;
   double h;
   SmVector3 x,lc,uc;
   this->lc = input_lc;
   this->uc = input_uc;
   this->max_level = max_level;
   dim = (int)pow2(max_level)+1;
   dim2 = dim*dim;
   nvertices = dim2*dim;
   DtCell *cell;
   SmVector3 v;
   phi.setNumValues(8);
   pts.setNumVecs(8);
   cells.setNumCells(1);
   cells[0].level = 0;
   cells[0].lc[0] = cells[0].lc[1] = cells[0].lc[2] = this->lc;
   cells[0].uc[0] = cells[0].uc[1] = cells[0].uc[2] = this->uc;
   cells[0].ilc[0] = cells[0].ilc[1] = cells[0].ilc[2] = 0;
   cells[0].iuc[0] = cells[0].iuc[1] = cells[0].iuc[2] = (int)pow2(max_level);
   cells[0].vertices[0] = 0;
   cells[0].vertices[1] = 1;
   cells[0].vertices[2] = 2;
   cells[0].vertices[3] = 3;
   cells[0].vertices[4] = 4;
   cells[0].vertices[5] = 5;
   cells[0].vertices[6] = 6;
   cells[0].vertices[7] = 7;
   cells[0].parent = -1;
   cells[0].child[0] = cells[0].child[1] = cells[0].child[2] = cells[0].child[3] =
	   cells[0].child[4] = cells[0].child[5] = cells[0].child[6] = cells[0].child[7] = -1;

   pts[0][0] = cells[0].lc[0];
   pts[0][1] = cells[0].lc[1];
   pts[0][2] = cells[0].lc[2];
   phi[0] = evalPhi(pts[0]);
   vertMap[0] = 0;

   pts[1][0] = cells[0].lc[0];
   pts[1][1] = cells[0].lc[1];
   pts[1][2] = cells[0].uc[2];
   phi[1] = evalPhi(pts[1]);
   vertMap[(dim-1)*dim2] = 1;

   pts[2][0] = cells[0].lc[0];
   pts[2][1] = cells[0].uc[1];
   pts[2][2] = cells[0].lc[2];
   phi[2] = evalPhi(pts[2]);
   vertMap[(dim-1)*dim] = 2;

   pts[3][0] = cells[0].lc[0];
   pts[3][1] = cells[0].uc[1];
   pts[3][2] = cells[0].uc[2];
   phi[3] = evalPhi(pts[3]);
   vertMap[(dim-1)*dim+(dim-1)*dim2] = 3;

   pts[4][0] = cells[0].uc[0];
   pts[4][1] = cells[0].lc[1];
   pts[4][2] = cells[0].lc[2];
   phi[4] = evalPhi(pts[4]);
   vertMap[dim-1] = 4;

   pts[5][0] = cells[0].uc[0];
   pts[5][1] = cells[0].lc[1];
   pts[5][2] = cells[0].uc[2];
   phi[5] = evalPhi(pts[5]);
   vertMap[(dim-1)+(dim-1)*dim2] = 5;

   pts[6][0] = cells[0].uc[0];
   pts[6][1] = cells[0].uc[1];
   pts[6][2] = cells[0].lc[2];
   phi[6] = evalPhi(pts[6]);
   vertMap[(dim-1)+(dim-1)*dim] = 6;

   pts[7][0] = cells[0].uc[0];
   pts[7][1] = cells[0].uc[1];
   pts[7][2] = cells[0].uc[2];
   phi[7] = evalPhi(pts[7]);
   vertMap[(dim-1)+(dim-1)*dim+(dim-1)*dim2] = 7;

   x = pts[0]+0.5*(this->uc-this->lc);
   cells[0].dist = evalPhi(x);

   while (c < cells.numCells()) {
	   cell = &(cells[c]);
	   h = pow2(-cell->level-1)*(this->uc-this->lc);
	   if (cell->level < max_level && (width)*h > fabs(cell->dist)) {
		   ih = (int)pow2(max_level-cell->level-1);
		   l = cell->level+1;
		   lc[0] = cell->lc[0]; lc[1] = cell->lc[1]; lc[2] = cell->lc[2];
		   uc[0] = cell->uc[0]; uc[1] = cell->uc[1]; uc[2] = cell->uc[2];
		   ilc[0] = cell->ilc[0]; ilc[1] = cell->ilc[1]; ilc[2] = cell->ilc[2];
		   iuc[0] = cell->iuc[0]; iuc[1] = cell->iuc[1]; iuc[2] = cell->iuc[2];

		   n = cells.numCells();
		   cells.setNumCells(n+8,1,0);
		   cell = &(cells[c]);
		   cell->child[0] = n;
		   cell->child[1] = n+1;
		   cell->child[2] = n+2;
		   cell->child[3] = n+3;
		   cell->child[4] = n+4;
		   cell->child[5] = n+5;
		   cell->child[6] = n+6;
		   cell->child[7] = n+7;	
		   p = phi.numValues();

		   i = (ilc[0])+dim*(ilc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[0] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[0] = vertMap[i];


		   i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]);
		   if (!vertMap.count(i)) {
			   vert[1] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[1] = vertMap[i];

		   i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[2] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[2] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(ilc[2]);
		   if (!vertMap.count(i)) {
			   vert[3] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = lc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[3] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[4] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[4] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(ilc[2]);
		   if (!vertMap.count(i)) {
			   vert[5] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[5] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[6] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[6] = vertMap[i];

		   i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(iuc[2]);
		   if (!vertMap.count(i)) {
			   vert[7] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[7] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(iuc[2]);
		   if (!vertMap.count(i)) {
			   vert[8] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = uc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[8] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(iuc[2]);
		   if (!vertMap.count(i)) {
			   vert[9] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[9] = vertMap[i];

		   i = (ilc[0])+dim*(iuc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[10] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[10] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(ilc[2]);
		   if (!vertMap.count(i)) {
			   vert[11] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = lc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[11] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[12] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[12] = vertMap[i];

		   i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(iuc[2]);
		   if (!vertMap.count(i)) {
			   vert[13] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = uc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[13] = vertMap[i];

		   i = (iuc[0])+dim*(ilc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[14] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = uc[0]; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[14] = vertMap[i];

		   i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]);
		   if (!vertMap.count(i)) {
			   vert[15] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[15] = vertMap[i];

		   i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[16] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[16] = vertMap[i];

		   i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(iuc[2]);
		   if (!vertMap.count(i)) {
			   vert[17] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[17] = vertMap[i];

		   i = (iuc[0])+dim*(iuc[1])+dim2*(ilc[2]+ih);
		   if (!vertMap.count(i)) {
			   vert[18] = p;
			   vertMap[i] = p;
			   pts.setNumVecs(p+1,1,0);
			   phi.setNumValues(p+1,1,0);
			   pts[p][0] = uc[0]; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
			   phi[p] = evalPhi(pts[p]);
			   p++;
		   } else vert[18] = vertMap[i];

		   cell = &(cells[n]);
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0];
		   cell->lc[1] = lc[1];
		   cell->lc[2] = lc[2];
		   cell->uc[0] = lc[0]+h;
		   cell->uc[1] = lc[1]+h;
		   cell->uc[2] = lc[2]+h;
		   cell->ilc[0] = ilc[0];
		   cell->ilc[1] = ilc[1];
		   cell->ilc[2] = ilc[2];
		   cell->iuc[0] = ilc[0]+ih;
		   cell->iuc[1] = ilc[1]+ih;
		   cell->iuc[2] = ilc[2]+ih;
		   cell->vertices[0] = cells[c].vertices[0];
		   cell->vertices[1] = vert[0];
		   cell->vertices[2] = vert[1];
		   cell->vertices[3] = vert[2];
		   cell->vertices[4] = vert[3];
		   cell->vertices[5] = vert[4];
		   cell->vertices[6] = vert[5];
		   cell->vertices[7] = vert[6];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0];
		   cell->lc[1] = lc[1];
		   cell->lc[2] = lc[2]+h;
		   cell->uc[0] = lc[0]+h;
		   cell->uc[1] = lc[1]+h;
		   cell->uc[2] = uc[2];
		   cell->ilc[0] = ilc[0];
		   cell->ilc[1] = ilc[1];
		   cell->ilc[2] = ilc[2]+ih;
		   cell->iuc[0] = ilc[0]+ih;
		   cell->iuc[1] = ilc[1]+ih;
		   cell->iuc[2] = iuc[2];
		   cell->vertices[0] = vert[0];
		   cell->vertices[1] = cells[c].vertices[1];
		   cell->vertices[2] = vert[2];
		   cell->vertices[3] = vert[7];
		   cell->vertices[4] = vert[4];
		   cell->vertices[5] = vert[8];
		   cell->vertices[6] = vert[6];
		   cell->vertices[7] = vert[9];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0];
		   cell->lc[1] = lc[1]+h;
		   cell->lc[2] = lc[2];
		   cell->uc[0] = lc[0]+h;
		   cell->uc[1] = uc[1];
		   cell->uc[2] = lc[2]+h;
		   cell->ilc[0] = ilc[0];
		   cell->ilc[1] = ilc[1]+ih;
		   cell->ilc[2] = ilc[2];
		   cell->iuc[0] = ilc[0]+ih;
		   cell->iuc[1] = iuc[1];
		   cell->iuc[2] = ilc[2]+ih;
		   cell->vertices[0] = vert[1];
		   cell->vertices[1] = vert[2];
		   cell->vertices[2] = cells[c].vertices[2];
		   cell->vertices[3] = vert[10];
		   cell->vertices[4] = vert[5];
		   cell->vertices[5] = vert[6];
		   cell->vertices[6] = vert[11];
		   cell->vertices[7] = vert[12];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0];
		   cell->lc[1] = lc[1]+h;
		   cell->lc[2] = lc[2]+h;
		   cell->uc[0] = lc[0]+h;
		   cell->uc[1] = uc[1];
		   cell->uc[2] = uc[2];
		   cell->ilc[0] = ilc[0];
		   cell->ilc[1] = ilc[1]+ih;
		   cell->ilc[2] = ilc[2]+ih;
		   cell->iuc[0] = ilc[0]+ih;
		   cell->iuc[1] = iuc[1];
		   cell->iuc[2] = iuc[2];
		   cell->vertices[0] = vert[2];
		   cell->vertices[1] = vert[7];
		   cell->vertices[2] = vert[10];
		   cell->vertices[3] = cells[c].vertices[3];
		   cell->vertices[4] = vert[6];
		   cell->vertices[5] = vert[9];
		   cell->vertices[6] = vert[12];
		   cell->vertices[7] = vert[13];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0]+h;
		   cell->lc[1] = lc[1];
		   cell->lc[2] = lc[2];
		   cell->uc[0] = uc[0];
		   cell->uc[1] = lc[1]+h;
		   cell->uc[2] = lc[2]+h;
		   cell->ilc[0] = ilc[0]+ih;
		   cell->ilc[1] = ilc[1];
		   cell->ilc[2] = ilc[2];
		   cell->iuc[0] = iuc[0];
		   cell->iuc[1] = ilc[1]+ih;
		   cell->iuc[2] = ilc[2]+ih;
		   cell->vertices[0] = vert[3];
		   cell->vertices[1] = vert[4];
		   cell->vertices[2] = vert[5];
		   cell->vertices[3] = vert[6];
		   cell->vertices[4] = cells[c].vertices[4];
		   cell->vertices[5] = vert[14];
		   cell->vertices[6] = vert[15];
		   cell->vertices[7] = vert[16];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0]+h;
		   cell->lc[1] = lc[1];
		   cell->lc[2] = lc[2]+h;
		   cell->uc[0] = uc[0];
		   cell->uc[1] = lc[1]+h;
		   cell->uc[2] = uc[2];
		   cell->ilc[0] = ilc[0]+ih;
		   cell->ilc[1] = ilc[1];
		   cell->ilc[2] = ilc[2]+ih;
		   cell->iuc[0] = iuc[0];
		   cell->iuc[1] = ilc[1]+ih;
		   cell->iuc[2] = iuc[2];
		   cell->vertices[0] = vert[4];
		   cell->vertices[1] = vert[8];
		   cell->vertices[2] = vert[6];
		   cell->vertices[3] = vert[9];
		   cell->vertices[4] = vert[14];
		   cell->vertices[5] = cells[c].vertices[5];
		   cell->vertices[6] = vert[16];
		   cell->vertices[7] = vert[17];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0]+h;
		   cell->lc[1] = lc[1]+h;
		   cell->lc[2] = lc[2];
		   cell->uc[0] = uc[0];
		   cell->uc[1] = uc[1];
		   cell->uc[2] = lc[2]+h;
		   cell->ilc[0] = ilc[0]+ih;
		   cell->ilc[1] = ilc[1]+ih;
		   cell->ilc[2] = ilc[2];
		   cell->iuc[0] = iuc[0];
		   cell->iuc[1] = iuc[1];
		   cell->iuc[2] = ilc[2]+ih;
		   cell->vertices[0] = vert[5];
		   cell->vertices[1] = vert[6];
		   cell->vertices[2] = vert[11];
		   cell->vertices[3] = vert[12];
		   cell->vertices[4] = vert[15];
		   cell->vertices[5] = vert[16];
		   cell->vertices[6] = cells[c].vertices[6];
		   cell->vertices[7] = vert[18];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell++;
		   cell->parent = c;
		   cell->level = l;
		   cell->lc[0] = lc[0]+h;
		   cell->lc[1] = lc[1]+h;
		   cell->lc[2] = lc[2]+h;
		   cell->uc[0] = uc[0];
		   cell->uc[1] = uc[1];
		   cell->uc[2] = uc[2];
		   cell->ilc[0] = ilc[0]+ih;
		   cell->ilc[1] = ilc[1]+ih;
		   cell->ilc[2] = ilc[2]+ih;
		   cell->iuc[0] = iuc[0];
		   cell->iuc[1] = iuc[1];
		   cell->iuc[2] = iuc[2];
		   cell->vertices[0] = vert[6];
		   cell->vertices[1] = vert[9];
		   cell->vertices[2] = vert[12];
		   cell->vertices[3] = vert[13];
		   cell->vertices[4] = vert[16];
		   cell->vertices[5] = vert[17];
		   cell->vertices[6] = vert[18];
		   cell->vertices[7] = cells[c].vertices[7];
		   cell->child[0] = cell->child[1] = cell->child[2] = 
			   cell->child[3] = cell->child[4] = cell->child[5] = 
			   cell->child[6] = cell->child[7] = -1;

		   cell = &(cells[n]);
		   for (i=0; i<8; i++, cell++) {
			   x = pts[cell->vertices[0]]+(0.5*h);
			   cell->dist = evalPhi(x);
		   }
	   }
	   c++;
   }
}

};


#endif