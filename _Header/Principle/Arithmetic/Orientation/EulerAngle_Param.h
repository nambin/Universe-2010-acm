#ifndef _ARITHMETIC_EULERANGLE_PARAM_H_
#define _ARITHMETIC_EULERANGLE_PARAM_H_

namespace arithmetic {

enum EulerAxis {X = 0, Y, Z};

// a: axis, p: parity, r: repeat, f: frame
// axis is one of X, Y, or Z
// parity is 0 if X->Y (in static), 1 if X->Z (in static)
// repeat is 0 if ZYX, 1 if XYX 
// frame is 0 if rotation axes are global, 1 if local
#define EulerOrder(a, p, r, f)	((((((a<<1)+p)<<1)+r)<<1)+f)

	// global axes 
#define EA_XYZ_G	EulerOrder(X, 0, 0, 0)	
#define EA_XYX_G	EulerOrder(X, 0, 1, 0)	
#define EA_XZX_G	EulerOrder(X, 1, 1, 0)
#define EA_XZY_G	EulerOrder(X, 1, 0, 0)
#define EA_YXZ_G	EulerOrder(Y, 1, 0, 0)
#define EA_YXY_G	EulerOrder(Y, 1, 1, 0)
#define EA_YZX_G	EulerOrder(Y, 0, 0, 0)
#define EA_YZY_G	EulerOrder(Y, 0, 1, 0)
#define EA_ZXY_G	EulerOrder(Z, 0, 0, 0)
#define EA_ZXZ_G	EulerOrder(Z, 0, 1, 0)
#define EA_ZYX_G	EulerOrder(Z, 1, 0, 0)
#define EA_ZYZ_G	EulerOrder(Z, 1, 1, 0)

// local axes
#define EA_XYZ_L	EulerOrder(Z, 1, 0, 1)
#define EA_XYX_L	EulerOrder(X, 0, 1, 1)
#define EA_XZX_L	EulerOrder(X, 1, 1, 1)
#define EA_XZY_L	EulerOrder(Y, 0, 0, 1)
#define EA_YXZ_L	EulerOrder(Z, 0, 0, 1)
#define EA_YXY_L	EulerOrder(Y, 1, 1, 1)
#define EA_YZX_L	EulerOrder(X, 1, 0, 1)
#define EA_YZY_L	EulerOrder(Y, 0, 1, 1)
#define EA_ZXY_L	EulerOrder(Y, 1, 0, 1)
#define EA_ZXZ_L	EulerOrder(Z, 0, 1, 1)
#define EA_ZYX_L	EulerOrder(X, 0, 0, 1)		// Roll-pitch-yaw method
#define EA_ZYZ_L	EulerOrder(Z, 1, 1, 1)

#define	EA_XY_G		(EA_XYZ_G + 0x18)
#define EA_XZ_G		(EA_XZY_G + 0x18)
#define EA_YX_G		(EA_YXZ_G + 0x18)
#define	EA_YZ_G		(EA_YZX_G + 0x18)
#define	EA_ZX_G		(EA_ZXY_G + 0x18)
#define EA_ZY_G		(EA_ZYX_G + 0x18)
#define	EA_XY_L		(EA_XYX_L + 0x18)
#define EA_XZ_L		(EA_XZX_L + 0x18)
#define EA_YX_L		(EA_YXZ_L + 0x18)
#define	EA_YZ_L		(EA_YZX_L + 0x18)
#define	EA_ZX_L		(EA_ZXY_L + 0x18)
#define EA_ZY_L		(EA_ZYX_L + 0x18)

#define	EA_X		(EA_XYZ_G + 0x30)
#define	EA_Y		(EA_YZX_G + 0x30)
#define EA_Z		(EA_ZXY_G + 0x30)

#define EulSafe	"\000\001\002\000"
#define EulNext	"\001\002\000\001"

};		// arithmetic

#endif