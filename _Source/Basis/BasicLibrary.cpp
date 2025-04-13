#include <stdafx.h>
#include <Basis/Basis/BasicLibrary.h>

namespace basis {

void byteSwap( unsigned char* b, int n ) {
	register int i = 0;
	register int j = n-1;
	while ( i < j ) {
		unsigned char tmp = b[j];
		b[j] = b[i];
		b[i] = tmp;
		i++, j--;
	}
}

};