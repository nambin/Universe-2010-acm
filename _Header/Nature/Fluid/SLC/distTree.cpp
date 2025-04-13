//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Distance Tree
//   -- Code File
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

#include <stdafx.h>
//#include "distTree.H"
#include <Nature/Fluid/SLC/distTree.H>
#include "float.h"
#include <cstdlib>
#ifdef PCUBE
extern "C" {
 int
fast_polygon_intersects_cube(int nverts, const double verts[][3],
			const double polynormal[3],
			int already_know_verts_are_outside_cube,
			int already_know_edges_are_outside_cube);
}
#endif

#pragma warning(disable:4101)
#pragma warning(disable:4018)

//namespace slcspace {

using namespace slcspace;

static const double EPS = 1.0e-15;
using namespace std;

const double width = 3;

//////////////////////////////////////////////////////////////////
// DtCell
//////////////////////////////////////////////////////////////////
bool DtCell::contains(const double *x) {
  if (x[0] < lc[0]-EPS || x[0] > uc[0]+EPS) return false;
  if (x[1] < lc[1]-EPS || x[1] > uc[1]+EPS) return false;
  if (x[2] < lc[2]-EPS || x[2] > uc[2]+EPS) return false;
  return true;
}

int DtCell::clamp(double *x) {
	int ret = 0;
  if (x[0] < lc[0]) { x[0] = lc[0]; ret += 1;}
  if (x[0] > uc[0]) { x[0] = uc[0]; ret += 2;}
  if (x[1] < lc[1]) { x[1] = lc[1]; ret += 4;}
  if (x[1] > uc[1]) { x[1] = uc[1]; ret += 8;}
  if (x[2] < lc[2]) { x[2] = lc[2]; ret += 16;}
  if (x[2] > uc[2]) { x[2] = uc[2]; ret += 32;}
  return ret;
}

namespace slcspace {

std::istream &operator>>(std::istream &strm, DtCell &c) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatStr("[",strm);
  strm >> c.level;
  eatChar(',',strm);
  strm >> c.lc[0];
  eatChar(',',strm);
  strm >> c.lc[1];
  eatChar(',',strm);
  strm >> c.lc[2];
  eatChar(',',strm);
  strm >> c.uc[0];
  eatChar(',',strm);
  strm >> c.uc[1];
  eatChar(',',strm);
  strm >> c.uc[2];
  eatChar(',',strm);
  strm >> c.vertices[0];
  eatChar(',',strm);
  strm >> c.vertices[1];
  eatChar(',',strm);
  strm >> c.vertices[2];
  eatChar(',',strm);
  strm >> c.vertices[3];
  eatChar(',',strm);
  strm >> c.vertices[4];
  eatChar(',',strm);
  strm >> c.vertices[5];
  eatChar(',',strm);
  strm >> c.vertices[6];
  eatChar(',',strm);
  strm >> c.vertices[7];
  eatChar(',',strm);
  strm >> c.parent;
  eatChar(',',strm);
  strm >> c.child[0];
  eatChar(',',strm);
  strm >> c.child[1];
  eatChar(',',strm);
  strm >> c.child[2];
  eatChar(',',strm);
  strm >> c.child[3];
  eatChar(',',strm);
  strm >> c.child[4];
  eatChar(',',strm);
  strm >> c.child[5];
  eatChar(',',strm);
  strm >> c.child[6];
  eatChar(',',strm);
  strm >> c.child[7];
  eatChar(',',strm);
  strm >> c.dist;
  eatChar(',',strm);
  strm >> c.triangles;
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

std::ostream &operator<<(std::ostream &strm, const DtCell &c) {
  strm << "[";
  strm << c.level;
  strm << ",";
  strm << c.lc[0];
  strm << ",";
  strm << c.lc[1];
  strm << ",";
  strm << c.lc[2];
  strm << ",";
  strm << c.uc[0];
  strm << ",";
  strm << c.uc[1];
  strm << ",";
  strm << c.uc[2];
  strm << ",";
  strm << c.vertices[0];
  strm << ",";
  strm << c.vertices[1];
  strm << ",";
  strm << c.vertices[2];
  strm << ",";
  strm << c.vertices[3];
  strm << ",";
  strm << c.vertices[4];
  strm << ",";
  strm << c.vertices[5];
  strm << ",";
  strm << c.vertices[6];
  strm << ",";
  strm << c.vertices[7];
  strm << ",";
  strm << c.parent;
  strm << ",";
  strm << c.child[0];
  strm << ",";
  strm << c.child[1];
  strm << ",";
  strm << c.child[2];
  strm << ",";
  strm << c.child[3];
  strm << ",";
  strm << c.child[4];
  strm << ",";
  strm << c.child[5];
  strm << ",";
  strm << c.child[6];
  strm << ",";
  strm << c.child[7];
  strm << ",";
  strm << c.dist;
  strm << ",";
  strm << c.triangles;
  strm << "]";
  return strm;
}  

};	// slcspace

void DtCell::read(ifstream &in) {
	in.read((char*) &level, sizeof(int));
  in.read((char*) &lc[0], sizeof(double));
  in.read((char*) &lc[1], sizeof(double));
  in.read((char*) &lc[2], sizeof(double));
  in.read((char*) &uc[0], sizeof(double));
  in.read((char*) &uc[1], sizeof(double));
  in.read((char*) &uc[2], sizeof(double));
  in.read((char*) &ilc[0], sizeof(double));
  in.read((char*) &ilc[1], sizeof(double));
  in.read((char*) &ilc[2], sizeof(double));
  in.read((char*) &iuc[0], sizeof(double));
  in.read((char*) &iuc[1], sizeof(double));
  in.read((char*) &iuc[2], sizeof(double));
  in.read((char*) &vertices[0], sizeof(int));
  in.read((char*) &vertices[1], sizeof(int));
  in.read((char*) &vertices[2], sizeof(int));
  in.read((char*) &vertices[3], sizeof(int));
  in.read((char*) &vertices[4], sizeof(int));
  in.read((char*) &vertices[5], sizeof(int));
  in.read((char*) &vertices[6], sizeof(int));
  in.read((char*) &vertices[7], sizeof(int));
  in.read((char*) &parent, sizeof(int));
  in.read((char*) &child[0], sizeof(int));
  in.read((char*) &child[1], sizeof(int));
  in.read((char*) &child[2], sizeof(int));
  in.read((char*) &child[3], sizeof(int));
  in.read((char*) &child[4], sizeof(int));
  in.read((char*) &child[5], sizeof(int));
  in.read((char*) &child[6], sizeof(int));
  in.read((char*) &child[7], sizeof(int));
	in.read((char*) &dist, sizeof(double));
  triangles.read(in);
}

void DtCell::write(ofstream &out) {
	out.write((char*) &level, sizeof(int));
  out.write((char*) &lc[0], sizeof(double));
  out.write((char*) &lc[1], sizeof(double));
  out.write((char*) &lc[2], sizeof(double));
  out.write((char*) &uc[0], sizeof(double));
  out.write((char*) &uc[1], sizeof(double));
  out.write((char*) &uc[2], sizeof(double));
  out.write((char*) &ilc[0], sizeof(double));
  out.write((char*) &ilc[1], sizeof(double));
  out.write((char*) &ilc[2], sizeof(double));
  out.write((char*) &iuc[0], sizeof(double));
  out.write((char*) &iuc[1], sizeof(double));
  out.write((char*) &iuc[2], sizeof(double));
  out.write((char*) &vertices[0], sizeof(int));
  out.write((char*) &vertices[1], sizeof(int));
  out.write((char*) &vertices[2], sizeof(int));
  out.write((char*) &vertices[3], sizeof(int));
  out.write((char*) &vertices[4], sizeof(int));
  out.write((char*) &vertices[5], sizeof(int));
  out.write((char*) &vertices[6], sizeof(int));
  out.write((char*) &vertices[7], sizeof(int));
  out.write((char*) &parent, sizeof(int));
  out.write((char*) &child[0], sizeof(int));
  out.write((char*) &child[1], sizeof(int));
  out.write((char*) &child[2], sizeof(int));
  out.write((char*) &child[3], sizeof(int));
  out.write((char*) &child[4], sizeof(int));
  out.write((char*) &child[5], sizeof(int));
  out.write((char*) &child[6], sizeof(int));
  out.write((char*) &child[7], sizeof(int));
	out.write((char*) &dist, sizeof(double));
  triangles.write(out);
}

DtCell &DtCell::operator=(const DtCell &that) {
  level = that.level;
  lc[0] = that.lc[0];   
  lc[1] = that.lc[1];   
  lc[2] = that.lc[2];
  uc[0] = that.uc[0];   
  uc[1] = that.uc[1];   
  uc[2] = that.uc[2];
  ilc[0] = that.ilc[0];   
  ilc[1] = that.ilc[1];   
  ilc[2] = that.ilc[2];
  iuc[0] = that.iuc[0];   
  iuc[1] = that.iuc[1];   
  iuc[2] = that.iuc[2];
  vertices[0] = that.vertices[0];
  vertices[1] = that.vertices[1];
  vertices[2] = that.vertices[2];
  vertices[3] = that.vertices[3];
  vertices[4] = that.vertices[4];
  vertices[5] = that.vertices[5];
  vertices[6] = that.vertices[6];
  vertices[7] = that.vertices[7];
  parent = that.parent;
  child[0] = that.child[0];
  child[1] = that.child[1];
  child[2] = that.child[2];
  child[3] = that.child[3];
  child[4] = that.child[4];
  child[5] = that.child[5];
  child[6] = that.child[6];
  child[7] = that.child[7];
  dist = that.dist;
	triangles = that.triangles;
	return (*this);
}

DtCell::DtCell(const DtCell &that) {
  level = that.level;
  lc[0] = that.lc[0];   
  lc[1] = that.lc[1];   
  lc[2] = that.lc[2];
  uc[0] = that.uc[0];   
  uc[1] = that.uc[1];   
  uc[2] = that.uc[2];
  vertices[0] = that.vertices[0];
  vertices[1] = that.vertices[1];
  vertices[2] = that.vertices[2];
  vertices[3] = that.vertices[3];
  vertices[4] = that.vertices[4];
  vertices[5] = that.vertices[5];
  vertices[6] = that.vertices[6];
  vertices[7] = that.vertices[7];
  parent = that.parent;
  child[0] = that.child[0];
  child[1] = that.child[1];
  child[2] = that.child[2];
  child[3] = that.child[3];
  child[4] = that.child[4];
  child[5] = that.child[5];
  child[6] = that.child[6];
  child[7] = that.child[7];
  dist = that.dist;
}
  
  
//////////////////////////////////////////////////////////////////
// DtCellList
//////////////////////////////////////////////////////////////////
DtCellList::DtCellList() {
  _numCells = 0;
  _numMalloc = 0;
  _cells = NULL;
}

DtCellList::DtCellList(const DtCellList &that) {
  if (_cells) delete [] _cells;
  _numCells = that._numCells;
  _numMalloc = that._numMalloc;
  _cells = new DtCell[_numMalloc];
  memcpy (_cells, that._cells, _numCells*sizeof(DtCell));
}

DtCellList &DtCellList::operator=(const DtCellList &that) {
  if (_cells) delete [] _cells;
  _numCells = that._numCells;
  _numMalloc = that._numMalloc;
  _cells = new DtCell[_numMalloc];
  memcpy(_cells, that._cells, _numCells*sizeof(DtCell));
  return (*this);
}

void DtCellList::setNumCells(unsigned int num, int copyOld, int exactSize) {
  unsigned int oldMalloc = _numMalloc;
  unsigned int newMalloc = num;

  if (exactSize == 0) {
    if ((newMalloc < (oldMalloc))&&(newMalloc > (oldMalloc/2))) {
      newMalloc = oldMalloc;
    }
    if ((newMalloc > (oldMalloc))&&(newMalloc < (oldMalloc*2))) {
      newMalloc = oldMalloc*2;
    }
  }

  if (newMalloc != oldMalloc) {
    unsigned int oldSize = _numCells;
    DtCell *oldList = _cells;
    _numCells = num;
    _numMalloc = newMalloc;
    _cells = new DtCell[_numMalloc];
    if (copyOld) 
			for (int i=0; i<((num>oldSize)?(oldSize):(num)); i++)
				_cells[i] = oldList[i];
    delete [] oldList;
  } else {
    _numCells = num;
  }
}

namespace slcspace {

istream &operator>>(std::istream &strm, DtCellList &l) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> newNum;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.setNumCells(newNum,0);
    for(i=0;i<newNum;i++) {
      strm >> l[i];
      if (i != (newNum-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const DtCellList &l) {
  strm << "[";
  strm << l.numCells();
  strm << ":";
  unsigned int i;
  for(i=0;i<l.numCells();i++) {
    strm << l[i];
    if (i != (l.numCells()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

};

void DtCellList::read(ifstream &in) {
	int i,l;
	in.read((char *)&l, sizeof(int));
	setNumCells(l);
	for (i=0; i<l; i++) _cells[i].read(in);
}

void DtCellList::write(ofstream &out) {
	int i;
	out.write((char *)&_numCells, sizeof(int));
	for (i=0; i<_numCells; i++) _cells[i].write(out);
}

//////////////////////////////////////////////////////////////////
// DtTree
//////////////////////////////////////////////////////////////////

int DtTree::clamp(SmVector3 &x) {
	return clamp(&(x[0]));
}

int DtTree::clamp(double *x) {
	int ret = 0;
  if (x[0] < lc) { x[0] = lc; ret += 1;}
  if (x[0] > uc) { x[0] = uc; ret += 2;}
  if (x[1] < lc) { x[1] = lc; ret += 4;}
  if (x[1] > uc) { x[1] = uc; ret += 8;}
  if (x[2] < lc) { x[2] = lc; ret += 16;}
  if (x[2] > uc) { x[2] = uc; ret += 32;}
  return ret;
}

DtCell* DtTree::findCell(const SmVector3 &x) {
  int i;
  double mid[3];
	DtCell *c=&(cells[0]);

  while (c->child[0] != -1) {
    mid[0] = (c->lc[0]+c->uc[0])/2.0;
    mid[1] = (c->lc[1]+c->uc[1])/2.0;
    mid[2] = (c->lc[2]+c->uc[2])/2.0;
    i = 0;
    if (x[0] >= mid[0]) i+=4;
    if (x[1] >= mid[1]) i+=2;
    if (x[2] >= mid[2]) i+=1;
    c = &(cells[c->child[i]]);
  }
	return c;
}

DtCell* DtTree::findCell(const SmVector3 &x, int &level) {
  int i;
  double mid[3];
	DtCell *c=&(cells[0]);

  while (c->child[0] != -1 && c->level < level) {
    mid[0] = (c->lc[0]+c->uc[0])/2.0;
    mid[1] = (c->lc[1]+c->uc[1])/2.0;
    mid[2] = (c->lc[2]+c->uc[2])/2.0;
    i = 0;
    if (x[0] >= mid[0]) i+=4;
    if (x[1] >= mid[1]) i+=2;
    if (x[2] >= mid[2]) i+=1;
    c = &(cells[c->child[i]]);
  }
	level = c->level;
	return c;
}

namespace slcspace {

istream &operator>>(std::istream &strm, DtTree &t) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatStr("[",strm);
  strm >> t.max_level;
  eatChar(',',strm);
  strm >> t.lc;
  eatChar(',',strm);
  strm >> t.uc;
  eatChar(',',strm);
  strm >> t.cells;
  eatChar(',',strm);
  strm >> t.nvertices;
  eatChar(',',strm);
  strm >> t.phi;
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const DtTree &t) {
  strm << "[";
  strm << t.max_level;
  strm << ",";
  strm << t.lc;
  strm << ",";
  strm << t.uc;
  strm << ",";
  strm << t.cells;
  strm << ",";
  strm << t.nvertices;
  strm << ",";
  strm << t.phi;
  strm << "]";
  return strm;
}

};

void DtTree::read(std::ifstream &in) {
	int version, index;
	in.read((char*)&version, sizeof(int));
	if (in.fail()) {
		cerr<<"can't read from file"<<endl;
		exit(-1);
	}
	if (version == 1) {
		in.read((char*)&max_level, sizeof(int));
		in.read((char*)&lc, sizeof(double));
		in.read((char*)&uc, sizeof(double));
		cells.read(in);
		in.read((char*)&nvertices, sizeof(int));
		phi.read(in);
	} else if (version == 2) {
		in.read((char*)&max_level, sizeof(int));
		in.read((char*)&lc, sizeof(double));
		in.read((char*)&uc, sizeof(double));
		cells.read(in);
		in.read((char*)&nvertices, sizeof(int));
		phi.read(in);
	} else {		
		cerr<<"incompatible bdt file version="<<version<<endl;
	}
	dim = (int)pow2(max_level)+1;
	dim2 = dim*dim;
	pts.setNumVecs(phi.numValues());
	for (int i=0; i<cells.numCells(); i++) {
		pts[cells[i].vertices[0]].set(cells[i].lc[0], cells[i].lc[1], cells[i].lc[2]);
		pts[cells[i].vertices[1]].set(cells[i].lc[0], cells[i].lc[1], cells[i].uc[2]);
		pts[cells[i].vertices[2]].set(cells[i].lc[0], cells[i].uc[1], cells[i].lc[2]);
		pts[cells[i].vertices[3]].set(cells[i].lc[0], cells[i].uc[1], cells[i].uc[2]);
		pts[cells[i].vertices[4]].set(cells[i].uc[0], cells[i].lc[1], cells[i].lc[2]);
		pts[cells[i].vertices[5]].set(cells[i].uc[0], cells[i].lc[1], cells[i].uc[2]);
		pts[cells[i].vertices[6]].set(cells[i].uc[0], cells[i].uc[1], cells[i].lc[2]);
		pts[cells[i].vertices[7]].set(cells[i].uc[0], cells[i].uc[1], cells[i].uc[2]);
		if (!vertMap.count(index=cells[i].ilc[0]+cells[i].ilc[1]*dim+cells[i].ilc[2]*dim2))
			vertMap[index] = cells[i].vertices[0];
		if (!vertMap.count(index=cells[i].ilc[0]+cells[i].ilc[1]*dim+cells[i].iuc[2]*dim2))
			vertMap[index] = cells[i].vertices[1];
		if (!vertMap.count(index=cells[i].ilc[0]+cells[i].iuc[1]*dim+cells[i].ilc[2]*dim2))
			vertMap[index] = cells[i].vertices[2];
		if (!vertMap.count(index=cells[i].ilc[0]+cells[i].iuc[1]*dim+cells[i].iuc[2]*dim2))
			vertMap[index] = cells[i].vertices[3];
		if (!vertMap.count(index=cells[i].iuc[0]+cells[i].ilc[1]*dim+cells[i].ilc[2]*dim2))
			vertMap[index] = cells[i].vertices[4];
		if (!vertMap.count(index=cells[i].iuc[0]+cells[i].ilc[1]*dim+cells[i].iuc[2]*dim2))
			vertMap[index] = cells[i].vertices[5];
		if (!vertMap.count(index=cells[i].iuc[0]+cells[i].iuc[1]*dim+cells[i].ilc[2]*dim2))
			vertMap[index] = cells[i].vertices[6];
		if (!vertMap.count(index=cells[i].iuc[0]+cells[i].iuc[1]*dim+cells[i].iuc[2]*dim2))
			vertMap[index] = cells[i].vertices[7];
	}
}

void DtTree::write(std::ofstream &out) {
	int version = 2;
	out.write((char*)&version, sizeof(int));
  out.write((char*)&max_level, sizeof(int));
  out.write((char*)&lc, sizeof(double));
  out.write((char*)&uc, sizeof(double));
  cells.write(out);
  out.write((char*)&nvertices, sizeof(int));
  phi.write(out);
}


void DtTree::buildTree(AbVecList &meshPts, AbTriList &triangles, AbVecList &faceNormals, 
											 double input_lc, double input_uc, int max_level) {
#ifndef PCUBE
	cerr<<"Distance trees cannot be built from a triangle mesh without the pcube library"<<endl;
#else

  int i,j,l,n,p,vert[19],ilc[3],iuc[3],ih,c=0;
	double h;
	SmVector3 x,lc,uc;
  this->lc = input_lc;
  this->uc = input_uc;
	this->max_level = max_level;
	dim = (int)pow2(max_level)+1;
	dim2 = dim*dim;
	nvertices = dim2*dim;
	DtCell *cell, *parent;
	SmVector3 v, bbl, bbu;
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
	cells[0].triangles.setNum(triangles.numTriangles());
	for (i=0; i<triangles.numTriangles(); i++) {
		cells[0].triangles[i] = i;
	}

	pts[0][0] = cells[0].lc[0];
	pts[0][1] = cells[0].lc[1];
	pts[0][2] = cells[0].lc[2];
	vertMap[0] = 0;

	pts[1][0] = cells[0].lc[0];
	pts[1][1] = cells[0].lc[1];
	pts[1][2] = cells[0].uc[2];
	vertMap[(dim-1)*dim2] = 1;
	
	pts[2][0] = cells[0].lc[0];
	pts[2][1] = cells[0].uc[1];
	pts[2][2] = cells[0].lc[2];
	vertMap[(dim-1)*dim] = 2;

	pts[3][0] = cells[0].lc[0];
	pts[3][1] = cells[0].uc[1];
	pts[3][2] = cells[0].uc[2];
	vertMap[(dim-1)*dim+(dim-1)*dim2] = 3;

	pts[4][0] = cells[0].uc[0];
	pts[4][1] = cells[0].lc[1];
	pts[4][2] = cells[0].lc[2];
	vertMap[dim-1] = 4;

	pts[5][0] = cells[0].uc[0];
	pts[5][1] = cells[0].lc[1];
	pts[5][2] = cells[0].uc[2];
	vertMap[(dim-1)+(dim-1)*dim2] = 5;

	pts[6][0] = cells[0].uc[0];
	pts[6][1] = cells[0].uc[1];
	pts[6][2] = cells[0].lc[2];
	vertMap[(dim-1)+(dim-1)*dim] = 6;

	pts[7][0] = cells[0].uc[0];
	pts[7][1] = cells[0].uc[1];
	pts[7][2] = cells[0].uc[2];
	vertMap[(dim-1)+(dim-1)*dim+(dim-1)*dim2] = 7;

	while (c < cells.numCells()) {
		cell = &(cells[c]);
		h = pow2(-cell->level-1)*(this->uc-this->lc);
		if (cell->level < max_level && cell->triangles.getNum()) {
			ih = (int)pow2(max_level-cell->level-1);
			l = cell->level+1;
			lc[0] = cell->lc[0]; lc[1] = cell->lc[1]; lc[2] = cell->lc[2];
			uc[0] = cell->uc[0]; uc[1] = cell->uc[1]; uc[2] = cell->uc[2];
			ilc[0] = cell->ilc[0]; ilc[1] = cell->ilc[1]; ilc[2] = cell->ilc[2];
			iuc[0] = cell->iuc[0]; iuc[1] = cell->iuc[1]; iuc[2] = cell->iuc[2];

			n = cells.numCells();
			cells.setNumCells(n+8,1,0);
			parent = cell = &(cells[c]);
			p = phi.numValues();
			
			i = (ilc[0])+dim*(ilc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[0] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[0] = vertMap[i];

			i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]);
			if (!vertMap.count(i)) {
				vert[1] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
				p++;
			} else vert[1] = vertMap[i];

			i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[2] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
				p++;
			} else vert[2] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(ilc[2]);
			if (!vertMap.count(i)) {
				vert[3] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = lc[2];
				p++;
			} else vert[3] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[4] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[4] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(ilc[2]);
			if (!vertMap.count(i)) {
				vert[5] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
				p++;
			} else vert[5] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[6] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
				p++;
			} else vert[6] = vertMap[i];

			i = (ilc[0])+dim*(ilc[1]+ih)+dim2*(iuc[2]);
			if (!vertMap.count(i)) {
				vert[7] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
				p++;
			} else vert[7] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1])+dim2*(iuc[2]);
			if (!vertMap.count(i)) {
				vert[8] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]; pts[p][2] = uc[2];
				p++;
			} else vert[8] = vertMap[i];

			i = (ilc[0]+ih)+dim*(ilc[1]+ih)+dim2*(iuc[2]);
			if (!vertMap.count(i)) {
				vert[9] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
				p++;
			} else vert[9] = vertMap[i];

			i = (ilc[0])+dim*(iuc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[10] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[10] = vertMap[i];

			i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(ilc[2]);
			if (!vertMap.count(i)) {
				vert[11] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = lc[2];
				p++;
			} else vert[11] = vertMap[i];

			i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[12] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[12] = vertMap[i];

			i = (ilc[0]+ih)+dim*(iuc[1])+dim2*(iuc[2]);
			if (!vertMap.count(i)) {
				vert[13] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = lc[0]+h; pts[p][1] = uc[1]; pts[p][2] = uc[2];
				p++;
			} else vert[13] = vertMap[i];

			i = (iuc[0])+dim*(ilc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[14] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = uc[0]; pts[p][1] = lc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[14] = vertMap[i];

			i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]);
			if (!vertMap.count(i)) {
				vert[15] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2];
				p++;
			} else vert[15] = vertMap[i];

			i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[16] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = lc[2]+h;
				p++;
			} else vert[16] = vertMap[i];

			i = (iuc[0])+dim*(ilc[1]+ih)+dim2*(iuc[2]);
			if (!vertMap.count(i)) {
				vert[17] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = uc[0]; pts[p][1] = lc[1]+h; pts[p][2] = uc[2];
				p++;
			} else vert[17] = vertMap[i];

			i = (iuc[0])+dim*(iuc[1])+dim2*(ilc[2]+ih);
			if (!vertMap.count(i)) {
				vert[18] = p;
				vertMap[i] = p;
				pts.setNumVecs(p+1,1,0);
				phi.setNumValues(p+1,1,0);
				pts[p][0] = uc[0]; pts[p][1] = uc[1]; pts[p][2] = lc[2]+h;
				p++;
			} else vert[18] = vertMap[i];

			cell->child[0] = n;
      cell->child[1] = n+1;
      cell->child[2] = n+2;
      cell->child[3] = n+3;
      cell->child[4] = n+4;
      cell->child[5] = n+5;
      cell->child[6] = n+6;
      cell->child[7] = n+7;

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

			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);

			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}

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
			for (i=0; i<parent->triangles.getNum(); i++) {
				j = parent->triangles[i];
				SmVector3 center=(cell->lc+cell->uc)/2.0;	
				double lh = (cell->uc[0]-cell->lc[0])*2.00;
				double tri[3][3];
				tri[0][0]=(meshPts[triangles[j].a][0]-center[0])/(lh);
				tri[0][1]=(meshPts[triangles[j].a][1]-center[1])/(lh);
				tri[0][2]=(meshPts[triangles[j].a][2]-center[2])/(lh);
				tri[1][0]=(meshPts[triangles[j].b][0]-center[0])/(lh);
				tri[1][1]=(meshPts[triangles[j].b][1]-center[1])/(lh);
				tri[1][2]=(meshPts[triangles[j].b][2]-center[2])/(lh);
				tri[2][0]=(meshPts[triangles[j].c][0]-center[0])/(lh);
				tri[2][1]=(meshPts[triangles[j].c][1]-center[1])/(lh);
				tri[2][2]=(meshPts[triangles[j].c][2]-center[2])/(lh);
				
				if (fast_polygon_intersects_cube(3,tri,&(faceNormals[j][0]),0,0))
					cell->triangles.add(j);
			}
		} 
		c++;
	} 
#endif
}

void DtTree::buildTree(double (*evalPhi)(const SmVector3 &x), 
											 double input_lc, double input_uc, int max_level) {
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

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: distTree.cpp,v $
// Revision 1.2  2005/06/10 20:08:47  adamb
// Various bug fixes, added iuc and ilc to cell output, rebuild vertMap after readining in a surface, fixed argc checking in the drivers.
//
// Revision 1.1.1.1  2005/05/25 05:41:01  adamb
// Initial Revision
//
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

