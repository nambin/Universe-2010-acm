//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Adam Bargteil's Useful Stuff
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
//#include "abUtil.H"
#include <Nature/Fluid/SLC/abUtil.H>
#include <iostream>
#include <fstream>

#pragma warning(disable:4101)
#pragma warning(disable:4018)

using namespace std;
using namespace slcspace;

namespace slcspace {

double pow2(int p) {
  switch (p) {
	case -20: return 9.53674e-07;
	case -19: return 1.90735e-06;
	case -18: return 3.8147e-06;
	case -17: return 7.62939e-06;
	case -16: return 1.52588e-05;
	case -15: return 3.05176e-05;
	case -14: return 6.10352e-05;
	case -13: return 0.0001220703125;
	case -12: return 0.000244140625;
	case -11: return 0.00048828125;
	case -10: return 0.0009765625;
	case -9: return 0.001953125;
	case -8: return 0.00390625;
	case -7: return 0.0078125;
	case -6: return 0.015625;
	case -5: return 0.03125;
	case -4: return 0.0625;
	case -3: return 0.125;
	case -2: return 0.25;
	case -1: return 0.5;
	case 0: return 1;
	case 1: return 2;
	case 2: return 4;
	case 3: return 8;
	case 4: return 16;
	case 5: return 32;
	case 6: return 64;
	case 7: return 128;
	case 8: return 256;
	case 9: return 512;
	case 10: return 1024;
	case 11: return 2048;
	case 12: return 4096;
	case 13: return 8192;
	case 14: return 16384;
	case 15: return 32768;
	case 16: return 65536;
	case 17: return 131072;
	case 18: return 262144;
	case 19: return 524288;
	case 20: return 1048576;
  default:
    double ret = 1;
	if (::abs(p) == p)
      for (int i=0; i<::abs(p); i++)
				ret *= 2.0;
    else
      for (int i=0; i<::abs(p); i++)
				ret /= 2.0;
    return ret;
  }
}

};

//////////////////////////////////////////////////////////////////
// AbIndexList
//////////////////////////////////////////////////////////////////

AbIndexList::AbIndexList() {
  _numValues = 0;
  _numMalloc = 0;
  _values = NULL;
}

AbIndexList::AbIndexList(const AbIndexList &that) {
  if (_values) delete [] _values;
  _numValues = that._numValues;
  _numMalloc = that._numMalloc;
  _values = new int[_numMalloc];
  memcpy (_values, that._values, _numValues*sizeof(int));
}

AbIndexList &AbIndexList::operator=(const AbIndexList &that) {
  if (_values) delete [] _values;
  _numValues = that._numValues;
  _numMalloc = that._numMalloc;
  _values = new int[_numMalloc];
  memcpy(_values, that._values, _numValues*sizeof(int));
  return (*this);
}

void AbIndexList::setNum(unsigned int num, int copyOld, int exactSize) {
  unsigned int oldMalloc = _numMalloc;
  unsigned int newMalloc = num;

  if (exactSize == 0) {
    if ((newMalloc < (oldMalloc))&&
	(newMalloc > (oldMalloc/2))) {
      newMalloc = oldMalloc;
    }
    if ((newMalloc > (oldMalloc)) &&
	(newMalloc < (oldMalloc*2))) {
      newMalloc = oldMalloc*2;
    }
  }

  if (newMalloc != oldMalloc) {
    unsigned int oldSize = _numValues;
    int *oldList = _values;
    _numValues = num;
    _numMalloc = newMalloc;
    _values = new int[_numMalloc];
    if (copyOld) 
      memcpy(_values, oldList, ((num>oldSize)?(oldSize):(num))*sizeof(int));
    delete [] oldList;
  } else {
    _numValues = num;
  }
}

void AbIndexList::add(int value) {
  setNum(_numValues+1,1,0);
  _values[_numValues-1] = value;
}

void AbIndexList::add(AbIndexList &list) {
  int j,i=_numValues;
  setNum(_numValues+list.getNum(),1,0);
	for (j=0; j<(int)list.getNum(); j++, i++)
		_values[i] = list[j];
}

void AbIndexList::add(int value, int index) {
  int i;
  setNum(_numValues+1,1,0);
  for (i=_numValues-1; i>index; i--)
    _values[i] = _values[i-1];
  _values[index] = value;
}

void AbIndexList::swap(int index1, int index2) {
  int tmp = _values[index1];
  _values[index1] = _values[index2];
  _values[index2] = tmp;
}

void AbIndexList::del(int index) {
  unsigned int i;
  for (i=index; i<(int)_numValues-1; i++)
    _values[i] = _values[i+1];
  setNum(_numValues-1,1,0);
}

void AbIndexList::delValue(int value) {
  unsigned int i;
  for (i=0; i<_numValues; i++)
    if (_values[i] == value)
      del(i);
}

void AbIndexList::clear() {
  if (_numValues != 0) {
    _numValues = 0;
    _numMalloc = 0;
    if (_values) delete [] _values;
    _values = NULL;
  }
}

int AbIndexList::contains(int value) {
  unsigned int i;
  int *v = _values;
  for (i=0; i<_numValues; i++, v++)
    if (*v == value) return 1;
  return 0;
}

int AbIndexList::find(int value) {
  unsigned int i;
  int *v = _values;
  for (i=0; i<_numValues; i++, v++)
    if (*v == value) return i;
  return -1;
}

namespace slcspace {

istream &operator>>(std::istream &strm, AbIndexList &l) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> newNum;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.setNum(newNum,0);
    for(i=0;i<newNum;i++) {
      strm >> l[i];
      if (i != (newNum-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const AbIndexList &l) {
  strm << "[";
  strm << l.getNum();
  strm << ":";
  unsigned int i;
  for(i=0;i<l.getNum();i++) {
    strm << l[i];
    if (i != (l.getNum()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

};

void AbIndexList::read(ifstream &in) {
	int i,l;
	in.read((char *)&l, sizeof(int));
	setNum(l);
	in.read((char*) _values, _numValues*sizeof(int));
}

void AbIndexList::write(ofstream &out) {
	int i;
	out.write((char*)&_numValues, sizeof(int));
	out.write((char*) _values, _numValues*sizeof(int));
}

//////////////////////////////////////////////////////////////////
// AbDoubleList
//////////////////////////////////////////////////////////////////

AbDoubleList::AbDoubleList() {
  _numValues = 0;
  _numMalloc = 0;
  _values = NULL;
}

AbDoubleList::AbDoubleList(const AbDoubleList &that) {
  if (_values) delete [] _values;
  _numValues = that._numValues;
  _numMalloc = that._numMalloc;
  _values = new double[_numMalloc];
  memcpy (_values, that._values, _numValues*sizeof(double));
}

AbDoubleList &AbDoubleList::operator=(const AbDoubleList &that) {
  if (_values) delete [] _values;
  _numValues = that._numValues;
  _numMalloc = that._numMalloc;
  _values = new double[_numMalloc];
  memcpy(_values, that._values, _numValues*sizeof(double));
  return (*this);
}

void AbDoubleList::setNumValues(unsigned int num, int copyOld, int exactSize) {
  unsigned int oldMalloc = _numMalloc;
  unsigned int newMalloc = num;

  if (exactSize == 0) {
    if ((newMalloc < (oldMalloc))&&
	(newMalloc > (oldMalloc/2))) {
      newMalloc = oldMalloc;
    }
    if ((newMalloc > (oldMalloc)) &&
	(newMalloc < (oldMalloc*2))) {
      newMalloc = oldMalloc*2;
    }
  }
  
  if (newMalloc != oldMalloc) {
    unsigned int oldSize = _numValues;
    double *oldList = _values;
    _numValues = num;
    _numMalloc = newMalloc;
    _values = new double[_numMalloc];
		memset(_values, 0, _numMalloc*sizeof(double));
    if (copyOld) 
      memcpy(_values, oldList, ((num>oldSize)?(oldSize):(num))*sizeof(double));
    delete [] oldList;
  } else {
    _numValues = num;
  }
}

namespace slcspace {

istream &operator>>(std::istream &strm, AbDoubleList &l) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> newNum;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.setNumValues(newNum,0);
    for(i=0;i<newNum;i++) {
      strm >> l[i];
      if (i != (newNum-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const AbDoubleList &l) {
  strm << "[";
  strm << l.numValues();
  strm << ":";
  unsigned int i;
  for(i=0;i<l.numValues();i++) {
    strm << l[i];
    if (i != (l.numValues()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

};

void AbDoubleList::read(ifstream &in) {
	int i,l;
	in.read((char *)&l, sizeof(int));
	setNumValues(l);
	in.read((char*) _values, _numValues*sizeof(double));
}

void AbDoubleList::write(ofstream &out) {
	int i;
	out.write((char *)&_numValues, sizeof(int));
	out.write((char*) _values, _numValues*sizeof(double));
}

namespace slcspace {

//////////////////////////////////////////////////////////////////
// AbTri
//////////////////////////////////////////////////////////////////
std::istream &operator>>(std::istream &strm, AbTri &t) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);
  eatStr("[",strm);
  strm >> t.a;
  eatChar(',',strm);
  strm >> t.b;
  eatChar(',',strm);
  strm >> t.c;
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}
std::ostream &operator<<(std::ostream &strm, const AbTri &t) {
  strm << "[";
  strm << t.a;
  strm << ",";
  strm << t.b;
  strm << ",";
  strm << t.c;
  strm << "]";
  return strm;
}

};

void AbTri::read(ifstream &in) {
	in.read((char*) &a, sizeof(int));
	in.read((char*) &b, sizeof(int));
	in.read((char*) &c, sizeof(int));
}
void AbTri::write(ofstream &out) {
	out.write((char*) &a, sizeof(int));
	out.write((char*) &b, sizeof(int));
	out.write((char*) &c, sizeof(int));
}

//////////////////////////////////////////////////////////////////
// AbTriList
//////////////////////////////////////////////////////////////////
AbTriList::AbTriList() {
  _numTriangles = 0;
  _numMalloc = 0;
  _triangles = NULL;
}

AbTriList::AbTriList(const AbTriList &that) {
  if (_triangles) delete [] _triangles;
  _numTriangles = that._numTriangles;
  _numMalloc = that._numMalloc;
  _triangles = new AbTri[_numMalloc];
  memcpy (_triangles, that._triangles, _numTriangles*sizeof(AbTri));
}

AbTriList &AbTriList::operator=(const AbTriList &that) {
  if (_triangles) delete [] _triangles;
  _numTriangles = that._numTriangles;
  _numMalloc = that._numMalloc;
  _triangles = new AbTri[_numMalloc];
  memcpy(_triangles, that._triangles, _numTriangles*sizeof(AbTri));
  return (*this);
}

int AbTriList::add(int a, int b, int c) {
	const int i = numTriangles();
	setNumTriangles(numTriangles()+1,1,0);
	_triangles[i].a = a;
	_triangles[i].b = b;
	_triangles[i].c = c;
	return i;
}

void AbTriList::setNumTriangles(unsigned int num, int copyOld, int exactSize) {
  unsigned int oldMalloc = _numMalloc;
  unsigned int newMalloc = num;

  if (exactSize == 0) {
    if ((newMalloc < (oldMalloc))&&
				(newMalloc > (oldMalloc/2))) {
      newMalloc = oldMalloc;
    }
    if ((newMalloc > (oldMalloc)) &&
				(newMalloc < (oldMalloc*2))) {
      newMalloc = oldMalloc*2;
    }
  }
  
  if (newMalloc != oldMalloc) {
    unsigned int oldSize = _numTriangles;
    AbTri *oldList = _triangles;
    _numTriangles = num;
    _numMalloc = newMalloc;
    _triangles = new AbTri[_numMalloc];
    if (copyOld) 
      memcpy(_triangles, oldList, ((num>oldSize)?(oldSize):(num))*sizeof(AbTri));
    delete [] oldList;
	} else {
    _numTriangles = num;
  }
}

#ifdef INVENTOR
SoIndexedFaceSet *AbTriList::faceSet() {
	int i;
	SoIndexedFaceSet *faces = new SoIndexedFaceSet();
	faces->ref();
	faces->coordIndex.setNum(4*numTriangles());
	int32_t *ptr = faces->coordIndex.startEditing();
	for (i=0; i<numTriangles(); i++) {
		*(ptr++) = _triangles[i].a;
		*(ptr++) = _triangles[i].b;
		*(ptr++) = _triangles[i].c;
		*(ptr++) = -1;
	}
	faces->coordIndex.finishEditing();
	faces->unrefNoDelete();
	return faces;
}
#endif

namespace slcspace {

istream &operator>>(std::istream &strm, AbTriList &l) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> newNum;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.setNumTriangles(newNum,0);
    for(i=0;i<newNum;i++) {
      strm >> l[i];
      if (i != (newNum-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const AbTriList &l) {
  strm << "[";
  strm << l.numTriangles();
  strm << ":";
  unsigned int i;
  for(i=0;i<l.numTriangles();i++) {
    strm << l[i];
    if (i != (l.numTriangles()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

};

void AbTriList::read(ifstream &in) {
	int i,l;
	in.read((char *)&l, sizeof(int));
	setNumTriangles(l);
	for (i=0; i<l; i++) _triangles[i].read(in);
}

void AbTriList::write(ofstream &out) {
	int i;
	out.write((char *)&_numTriangles, sizeof(int));
	for (i=0; i<_numTriangles; i++) _triangles[i].write(out);
}


//////////////////////////////////////////////////////////////////
// AbVecList
//////////////////////////////////////////////////////////////////
AbVecList::AbVecList() {
  _numVecs = 0;
  _numMalloc = 0;
  _vecs = NULL;
}

AbVecList::AbVecList(const AbVecList &that) {
  if (_vecs) delete [] _vecs;
  _numVecs = that._numVecs;
  _numMalloc = that._numMalloc;
  _vecs = new SmVector3[_numMalloc];
  memcpy (_vecs, that._vecs, _numVecs*sizeof(SmVector3));
}

AbVecList &AbVecList::operator=(const AbVecList &that) {
  if (_vecs) delete [] _vecs;
  _numVecs = that._numVecs;
  _numMalloc = that._numMalloc;
  _vecs = new SmVector3[_numMalloc];
  memcpy(_vecs, that._vecs, _numVecs*sizeof(SmVector3));
  return (*this);
}

void AbVecList::clear() {
	if (_vecs) delete [] _vecs;
	_vecs = NULL;
	_numVecs = _numMalloc = 0;
}

void AbVecList::setNumVecs(unsigned int num, int copyOld, int exactSize) {
  unsigned int oldMalloc = _numMalloc;
  unsigned int newMalloc = num;

  if (exactSize == 0) {
    if ((newMalloc < (oldMalloc))&&
				(newMalloc > (oldMalloc/2))) {
      newMalloc = oldMalloc;
    }
    if ((newMalloc > (oldMalloc)) &&
				(newMalloc < (oldMalloc*2))) {
      newMalloc = oldMalloc*2;
    }
  }
  
  if (newMalloc != oldMalloc) {
    unsigned int oldSize = _numVecs;
    SmVector3 *oldList = _vecs;
    _numVecs = num;
    _numMalloc = newMalloc;
    _vecs = new SmVector3[_numMalloc];
    if (copyOld) 
      memcpy(_vecs, oldList, ((num>oldSize)?(oldSize):(num))*sizeof(SmVector3));
    if (oldList) delete [] oldList;
  } else {
    _numVecs = num;
  }
}

#ifdef INVENTOR
SoCoordinate3 *AbVecList::coord() {
	int i;
	SoCoordinate3 *coords = new SoCoordinate3();
	coords->ref();

	if (coords->point.getNum() != numVecs())
		coords->point.setNum(numVecs());

	SbVec3f *ptr = coords->point.startEditing();
	SmVector3 *p = _vecs;
	for (i=0; i<numVecs(); i++,p++) 
		(ptr++)->setValue(float((*p)[0]), float((*p)[1]), float((*p)[2]));
	coords->point.finishEditing();
	coords->unrefNoDelete();
	return coords;
}
#endif

namespace slcspace {

istream &operator>>(std::istream &strm, AbVecList &l) {
  unsigned int i,newNum = 0;
  ios::fmtflags orgFlags = strm.setf(ios::skipws);

  eatStr("[",strm);
  strm >> newNum;
  eatStr(":",strm);
  
  if (strm.good()) {
    l.setNumVecs(newNum,0);
    for(i=0;i<newNum;i++) {
      strm >> l[i];
      if (i != (newNum-1)) eatStr(",",strm);
    }
  }
  eatChar(']',strm);
  strm.flags(orgFlags);
  return strm;
}

ostream &operator<<(ostream &strm,const AbVecList &l) {
  strm << "[";
  strm << l.numVecs();
  strm << ":";
  unsigned int i;
  for(i=0;i<l.numVecs();i++) {
    strm << l[i];
    if (i != (l.numVecs()-1)) strm << ",";
  }
  strm << "]";
  return strm;
}

};

void AbVecList::read(ifstream &in) {
	int i,l;
	in.read((char*)&l, sizeof(int));
	setNumVecs(l);
	for (i=0; i<l; i++) {
		in.read((char*)&_vecs[i][0], sizeof(double));
		in.read((char*)&_vecs[i][1], sizeof(double));
		in.read((char*)&_vecs[i][2], sizeof(double));
	}
}

void AbVecList::write(ofstream &out) {
	int i;
	out.write((char*)&_numVecs, sizeof(int));
	for (i=0; i<_numVecs; i++) {
		out.write((char*)&_vecs[i][0], sizeof(double));
		out.write((char*)&_vecs[i][1], sizeof(double));
		out.write((char*)&_vecs[i][2], sizeof(double));
	}
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: abUtil.cpp,v $
// Revision 1.1.1.1  2005/05/25 05:40:59  adamb
// Initial Revision
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------

