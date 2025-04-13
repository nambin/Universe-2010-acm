#pragma once

namespace arithmetic {

// Normalized coordinate
// 0 <= x <= 1
// 0 <= y <= 1
template <typename T, typename S> 
void multiCubic2D(T f00, T f10, T f01, T f11, T dx00, T dx10, T dx01, T dx11, T dy00, T dy10, T dy01, T dy11, S x, S y, T &f, T &dx, T &dy) {

	T C31, C13, C30, C21, C12, C03, C20, C11, C02, C10, C01, C00;

	C00 = f00;
	C10 = dx00;
	C01 = dy00;

	C20 = 3.f*(f10-f00) - dx10 - 2.f*dx00;
	C30 = -2.f*(f10-f00) + dx10 + dx00;

	C02 = 3.f*(f01-f00) - dy01 - 2.f*dy00;
	C03 = -2.f*(f01-f00) + dy01 + dy00;

	C21 = 3.f*f11 - 2.f*dx01 - dx11 - C20 - 3.f*(C03 + C02 + C01 + C00);
	C31 = -2.f*f11 + dx01 + dx11 - C30 + 2.f*(C03 + C02 + C01 + C00);

	C12 = 3.f*f11 - 2.f*dy10 - dy11 - C02 - 3.f*(C30 + C20 + C10 + C00);
	C13 = -2.f*f11 + dy10 + dy11 - C03 + 2.f*(C30 + C20 + C10 + C00);

	C11 = dx01 - C13 - C12 - C10;

	f = C31*x*x*x*y + C13*x*y*y*y + C30*x*x*x + C21*x*x*y + C12*x*y*y + C03*y*y*y + C20*x*x + C11*x*y + C02*y*y + C10*x + C01*y + C00;
	dx = 3.f*C31*x*x*y + C13*y*y*y + 3.f*C30*x*x + 2.f*C21*x*y + C12*y*y + 2.f*C20*x + C11*y + C10;
	dy = 3.f*C13*y*y*x + C31*x*x*x + 3.f*C03*y*y + 2.f*C12*y*x + C21*x*x + 2.f*C02*y + C11*x + C01;
}

};		// namespace arithmetic