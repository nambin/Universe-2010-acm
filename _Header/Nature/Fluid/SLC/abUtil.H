//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// Adam Bargteil's Useful Stuff
//   -- Header File
//  
// Primary Author: Adam Bargteil (adamb@cs.berkeley.edu)
// 
// This file contains a bunch of useful functions, and classes.
// List of doubles, integers, vectors, tirangles, etc.
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

#ifndef abUtil_H
#define abUtil_H

#include <iostream>
#include "smVector.H"
#include "../../Nature.h"

#ifdef _WIN32
#if _MSC_VER >= 1300
#include <hash_map>
#include <hash_set>
#else
#include <map>
#include <set>
#endif
#else
#include <ext/hash_map>
#include <ext/hash_set>
#endif


#ifdef INVENTOR
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#endif

namespace slcspace {

class NATUREDLL AbIndexList {
 public:
  AbIndexList();
  AbIndexList(const AbIndexList &that);
  ~AbIndexList() { delete [] _values;};

  AbIndexList &operator=(const AbIndexList &that);

  inline int &index(int i) {return _values[i];};
  inline const int &index(int i) const {return _values[i];};
  inline int &operator[](int i) {return _values[i];};
  inline const int &operator[](int i) const {return _values[i];};

  inline unsigned int getNum() const {return _numValues;};
  void setNum(unsigned int num, int copyOld = 0, int exactSize = 1);
  void add(int value);
  void add(AbIndexList &list);
  void add(int value, int index);
  void del(int index);
  void delValue(int value);
  void swap(int index1, int index2);
  void clear();
  int contains(int value);
  inline int pop();
  int find(int value);

	void read(std::ifstream &in);
	void write(std::ofstream &out);

 private:
  int _numValues;
  int _numMalloc;
  int *_values;
};
NATUREDLL std::istream &operator>>(std::istream &strm, AbIndexList &p);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const AbIndexList &p);

///////////////////////////////////////////////////////////////////
  
class NATUREDLL AbDoubleList {
 public:
  AbDoubleList();
  AbDoubleList(const AbDoubleList &that);
  ~AbDoubleList() { if (_values) delete [] _values;};

  AbDoubleList &operator=(const AbDoubleList &that);

  inline double &index(int i) {return _values[i];};
  inline const double &index(int i) const {return _values[i];};
  inline double &operator[](int i) {return _values[i];};
  inline const double &operator[](int i) const {return _values[i];};

  inline unsigned int numValues() const {return _numValues;};
  void setNumValues(unsigned int num, int copyOld = 0, int exactSize = 1);

	void read(std::ifstream &in);
	void write(std::ofstream &out);

 private:
  int _numValues;
  int _numMalloc;
  double *_values;
};
NATUREDLL std::istream &operator>>(std::istream &strm, AbDoubleList &p);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const AbDoubleList &p);
  
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

class NATUREDLL AbTri {
public:
	int a, b, c;
  AbTri &operator=(const AbTri &that) {a=that.a; b=that.b; c=that.c; return *this;};

	void read(std::ifstream &in);
	void write(std::ofstream &out);
};

NATUREDLL std::istream &operator>>(std::istream &strm, AbTri &t);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const AbTri &t);

class NATUREDLL AbTriList {
 public:
  AbTriList();
  AbTriList(const AbTriList &that);
  ~AbTriList() { delete [] _triangles;};

  AbTriList &operator=(const AbTriList &that);

  inline AbTri &index(int i) {return _triangles[i];};
  inline const AbTri &index(int i) const {return _triangles[i];};
  inline AbTri &operator[](int i) {return _triangles[i];};
  inline const AbTri &operator[](int i) const {return _triangles[i];};

  inline unsigned int numTriangles() const {return _numTriangles;};
  void setNumTriangles(unsigned int num, int copyOld = 0, int exactSize = 1);
	
	int add(int a, int b, int c);

#ifdef INVENTOR
	SoIndexedFaceSet *faceSet();
#endif

	void read(std::ifstream &in);
	void write(std::ofstream &out);

 private:
  int _numTriangles;
  int _numMalloc;
  AbTri *_triangles;
};
NATUREDLL std::istream &operator>>(std::istream &strm, AbTriList &p);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const AbTriList &p);
  
/////////////////////////////////////////////////////////////////

class NATUREDLL AbVecList {
 public:
  AbVecList();
  AbVecList(const AbVecList &that);
  ~AbVecList() { delete [] _vecs;};

  AbVecList &operator=(const AbVecList &that);

  inline SmVector3 &index(int i) {return _vecs[i];};
  inline const SmVector3 &index(int i) const {return _vecs[i];};
  inline SmVector3 &operator[](int i) {return _vecs[i];};
  inline const SmVector3 &operator[](int i) const {return _vecs[i];};

  inline unsigned int numVecs() const {return _numVecs;};
  void setNumVecs(unsigned int num, int copyOld = 0, int exactSize = 1);
	void clear();

#ifdef INVENTOR
	SoCoordinate3 *coord();
#endif

	void read(std::ifstream &in);
	void write(std::ofstream &out);

 private:
  int _numVecs;
  int _numMalloc;
  SmVector3 *_vecs;
};
NATUREDLL std::istream &operator>>(std::istream &strm, AbVecList &p);
NATUREDLL std::ostream &operator<<(std::ostream &strm, const AbVecList &p);
  

inline double absmin2 (double x, double y) {
  if (fabs(x) < fabs(y)) return x;
  return y;
}

inline double absmax2 (double x, double y) {
  if (fabs(x) > fabs(y)) return x;
  return y;
}

inline int sign(double x) {
  if (x>0.0) return 1;
	else return -1;
  if (x<0.0) return -1;
  if (x==0.0) return 0;
}

inline double max3 (double x, double y, double z) {
  if (x > y && x > z) return x;
  if (y > z) return y;
  return z;
}

inline double min3 (double x, double y, double z) {
  if (x < y && x < z) return x;
  if (y < z) return y;
  return z;
}

inline int max2 (int x, int y) {
  if (x>y) return x;
  else return y;
}

inline int min2 (int x, int y) {
  if (x<y) return x;
  else return y;
}

inline double min2 (double x, double y) {
  if (x<y) return x;
  else return y;
}

inline double max5(double x1, double x2, double x3, double x4, double x5) {
  double m = x1;
  if (x2 > m) m = x2;
  if (x3 > m) m = x3;
  if (x4 > m) m = x4;
  if (x5 > m) m = x5;
  return m;
}

NATUREDLL double pow2(int p);

inline static
std::istream &eatChar(char c,std::istream &buf) {
  char r;
  buf >> r;
  if (r!=c) {
    buf.clear(buf.rdstate() | std::ios::failbit);
  }
  return buf;
}

inline static
std::istream &eatStr(const char *s,std::istream &buf) {
  while (*s != '\0') {
    eatChar(*s,buf);
    s++;
  }
  return buf;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//
// RCS Revision History
//
// $Log: abUtil.H,v $
// Revision 1.2  2005/06/10 20:08:47  adamb
// Various bug fixes, added iuc and ilc to cell output, rebuild vertMap after readining in a surface, fixed argc checking in the drivers.
//
// Revision 1.1.1.1  2005/05/25 05:40:59  adamb
// Initial Revision
//
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------


};

#endif