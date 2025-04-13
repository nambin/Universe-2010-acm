#pragma once

#include "../../Arithmetic/MathLibrary/Interpolation_upwind_adaptive.h"

namespace geometryspace {

template<typename SRL>
typename const Derivative_SRL_2nd<SRL>::T Derivative_SRL_2nd<SRL>::dx(const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_x_left_1(node, p, val, h); val[1] = node(p); s.neighbor_x_right_1(node, p, &val[2], &h[1]);	
	return arithmetic::cd2_adative(val, h);
}

template<typename SRL>
typename const Derivative_SRL_2nd<SRL>::T Derivative_SRL_2nd<SRL>::dy(const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_y_left_1(node, p, val, h); val[1] = node(p); s.neighbor_y_right_1(node, p, &val[2], &h[1]);	
	return arithmetic::cd2_adative(val, h);
}

template<typename SRL>
typename const Derivative_SRL_2nd<SRL>::T Derivative_SRL_2nd<SRL>::dz(const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_z_left_1(node, p, val, h); val[1] = node(p); s.neighbor_z_right_1(node, p, &val[2], &h[1]);	
	return arithmetic::cd2_adative(val, h);
}

};