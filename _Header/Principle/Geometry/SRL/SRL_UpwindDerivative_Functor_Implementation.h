#pragma once

namespace geometryspace {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template<typename SRL>
void Upwind_SRL_1st<SRL>::dx(T& dx_m, T& dx_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	s.neighbor_x_left_1(node, p, val, h);	dx_m = (node(p) - val[0]) / h[0];
	s.neighbor_x_right_1(node, p, val, h);	dx_p = (val[0] - node(p)) / h[0];
}

template<typename SRL>
void Upwind_SRL_1st<SRL>::dy(T& dy_m, T& dy_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	s.neighbor_y_left_1(node, p, val, h);	dy_m = (node(p) - val[0]) / h[0];	
	s.neighbor_y_right_1(node, p, val, h);	dy_p = (val[0] - node(p)) / h[0];	
}

template<typename SRL>
void Upwind_SRL_1st<SRL>::dx(T& dx_m, T& dx_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_x_left_1(node, p, val, h);	dx_m = (node(p) - val[0]) / h[0];
	s.neighbor_x_right_1(node, p, val, h);	dx_p = (val[0] - node(p)) / h[0];
}

template<typename SRL>
void Upwind_SRL_1st<SRL>::dy(T& dy_m, T& dy_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_y_left_1(node, p, val, h);	dy_m = (node(p) - val[0]) / h[0];	
	s.neighbor_y_right_1(node, p, val, h);	dy_p = (val[0] - node(p)) / h[0];	
}

template<typename SRL>
void Upwind_SRL_1st<SRL>::dz(T& dz_m, T& dz_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_z_left_1(node, p, val, h);	dz_m = (node(p) - val[0]) / h[0];	
	s.neighbor_z_right_1(node, p, val, h);	dz_p = (val[0] - node(p)) / h[0];	
}

template<typename SRL>
void Upwind_SRL_2nd<SRL>::dx(T& dx_m, T& dx_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	s.neighbor_x_left_2(node, p, val, h); val[2] = node(p); s.neighbor_x_right_2(node, p, &val[3], &h[2]);	
	eno2nd_adative(dx_m, dx_p, val, h);	
}

template<typename SRL>
void Upwind_SRL_2nd<SRL>::dy(T& dy_m, T& dy_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	s.neighbor_y_left_2(node, p, val, h); val[2] = node(p); s.neighbor_y_right_2(node, p, &val[3], &h[2]);	
	eno2nd_adative(dy_m, dy_p, val, h);	
}

template<typename SRL>
void Upwind_SRL_2nd<SRL>::dx(T& dx_m, T& dx_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_x_left_2(node, p, val, h); val[2] = node(p); s.neighbor_x_right_2(node, p, &val[3], &h[2]);	
	eno2nd_adative(dx_m, dx_p, val, h);	
}

template<typename SRL>
void Upwind_SRL_2nd<SRL>::dy(T& dy_m, T& dy_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_y_left_2(node, p, val, h); val[2] = node(p); s.neighbor_y_right_2(node, p, &val[3], &h[2]);	
	eno2nd_adative(dy_m, dy_p, val, h);	
}

template<typename SRL>
void Upwind_SRL_2nd<SRL>::dz(T& dz_m, T& dz_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	s.neighbor_z_left_2(node, p, val, h); val[2] = node(p); s.neighbor_z_right_2(node, p, &val[3], &h[2]);	
	eno2nd_adative(dz_m, dz_p, val, h);	
}

template<typename SRL>
void Upwind_SRL_WENO3rd<SRL>::dx(T& dx_m, T& dx_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	if(s.neighbor_x_left_2(node, p, val) == false || s.neighbor_x_right_2(node, p, &val[3]) == false) {
		upwind_2nd_srl.dx(dx_m, dx_p, node, p);
		return;
	}

	val[2] = node(p);
	s.weno3rd_adaptive[p[0]].dx(dx_m, dx_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO3rd<SRL>::dy(T& dy_m, T& dy_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	if(s.neighbor_y_left_2(node, p, val) == false || s.neighbor_y_right_2(node, p, &val[3]) == false) {
		upwind_2nd_srl.dy(dy_m, dy_p, node, p);
		return;
	}

	val[2] = node(p);
	s.weno3rd_adaptive[p[1]].dx(dy_m, dy_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO3rd<SRL>::dx(T& dx_m, T& dx_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_x_left_2(node, p, val) == false || s.neighbor_x_right_2(node, p, &val[3]) == false) {		
		upwind_2nd_srl.dx(dx_m, dx_p, node, p);
		return;
	}

	val[2] = node(p);
	s.weno3rd_adaptive[p[0]].dx(dx_m, dx_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO3rd<SRL>::dy(T& dy_m, T& dy_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_y_left_2(node, p, val) == false || s.neighbor_y_right_2(node, p, &val[3]) == false) {		
		upwind_2nd_srl.dy(dy_m, dy_p, node, p);
		return;
	}

	val[2] = node(p);
	s.weno3rd_adaptive[p[1]].dx(dy_m, dy_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO3rd<SRL>::dz(T& dz_m, T& dz_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_z_left_2(node, p, val) == false || s.neighbor_z_right_2(node, p, &val[3]) == false) {		
		upwind_2nd_srl.dz(dz_m, dz_p, node, p);
		return;
	}

	val[2] = node(p);
	s.weno3rd_adaptive[p[2]].dx(dz_m, dz_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO<SRL>::dx(T& dx_m, T& dx_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	if(s.neighbor_x_left_3(node, p, val) == false || s.neighbor_x_right_3(node, p, &val[4]) == false) {
		upwind_2nd_srl.dx(dx_m, dx_p, node, p);
		return;
	}

	val[3] = node(p);
	s.weno_adaptive[p[0]].dx(dx_m, dx_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO<SRL>::dy(T& dy_m, T& dy_p, const SRLNode2D<T>& node, unsigned int p[2]) const {
	if(s.neighbor_y_left_3(node, p, val) == false || s.neighbor_y_right_3(node, p, &val[4]) == false) {
		upwind_2nd_srl.dy(dy_m, dy_p, node, p);
		return;
	}

	val[3] = node(p);
	s.weno_adaptive[p[1]].dx(dy_m, dy_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO<SRL>::dx(T& dx_m, T& dx_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_x_left_3(node, p, val) == false || s.neighbor_x_right_3(node, p, &val[4]) == false) {
		//upwind_2nd_srl.dx(dx_m, dx_p, node, p);
		upwind_1st_srl.dx(dx_m, dx_p, node, p);
		return;
	}

	val[3] = node(p);
	s.weno_adaptive[p[0]].dx(dx_m, dx_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO<SRL>::dy(T& dy_m, T& dy_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_y_left_3(node, p, val) == false || s.neighbor_y_right_3(node, p, &val[4]) == false) {
		//upwind_2nd_srl.dy(dy_m, dy_p, node, p);
		upwind_1st_srl.dy(dy_m, dy_p, node, p);
		return;
	}

	val[3] = node(p);
	s.weno_adaptive[p[1]].dx(dy_m, dy_p, val);
}

template<typename SRL>
void Upwind_SRL_WENO<SRL>::dz(T& dz_m, T& dz_p, const SRLNode3D<T>& node, unsigned int p[3]) const {
	if(s.neighbor_z_left_3(node, p, val) == false || s.neighbor_z_right_3(node, p, &val[4]) == false) {
		//upwind_2nd_srl.dz(dz_m, dz_p, node, p);
		upwind_1st_srl.dz(dz_m, dz_p, node, p);
		return;
	}

	val[3] = node(p);
	s.weno_adaptive[p[2]].dx(dz_m, dz_p, val);
}

#undef __ASSERT

};			// geometryspace