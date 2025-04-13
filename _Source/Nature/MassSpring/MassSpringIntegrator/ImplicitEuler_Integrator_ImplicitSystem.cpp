#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/ImplicitEuler_Integrator.h>

using namespace std;
using namespace springspace;

void ImplicitEuler_Integrator::allocateSystem(const QuadDynamicMesh<real>& mesh) {
	mesh.allocate_SystemMatrix(mat);
	mat.make_MatrixIndexMapping(stretch_matIdx, mesh.stretchSpring);
	mat.make_MatrixIndexMapping(shear_matIdx, mesh.shearSpring);
	mat.make_MatrixIndexMapping(bending_matIdx, mesh.bendingSpring);

	b.resize(mesh.nX * mesh.nY);
	dv.resize(mesh.nX * mesh.nY);
	f.resize(mesh.nX * mesh.nY);
}

void ImplicitEuler_Integrator::allocateSystem(const TriDynamicMesh<real>& mesh) {
	mesh.allocate_SystemMatrix(mat);
	mat.make_MatrixIndexMapping(stretch_matIdx, mesh.stretchSpring);

	b.resize(mesh.size());
	dv.resize(mesh.size());
	f.resize(mesh.size());
}

void ImplicitEuler_Integrator::makeSystem(real dt, const QuadDynamicMesh<real>& mesh) {
	if(mat.n == 0) {
		basis::throwError("Error in ImplicitEuler_Integrator::make_SystemMatrix() _ Should allocate System First");
		return;
	}

	makeSystem_Matrix(dt, mesh);
	makeSystem_RHS(dt, mesh);
}

void ImplicitEuler_Integrator::makeSystem(real dt, const TriDynamicMesh<real>& mesh) {
	if(mat.n == 0) {
		basis::throwError("Error in ImplicitEuler_Integrator::make_SystemMatrix() _ Should allocate System First");
		return;
	}

	makeSystem_Matrix(dt, mesh);
	makeSystem_RHS(dt, mesh);
}

void updateSystem_Matrix(MassSpringSystemMatrix<real>& mat, const Vec3Array<real> &pos0, const Vec3Array<real> &vel0, real dt, 
						const Array<Spring<real>>& spring, const SpringFunctor<real>* func, const Array<std::pair<size_t, size_t>>& matIdx) {
	real dt2 = dt*dt;
	Matrix3<real> dfdv, dfdx, _sum;

	for(size_t i=0;i<spring.size();++i) {
		const Spring<real>& s = spring[i];	

		func->DfDx(dfdx, s, pos0, vel0);
		func->DfDv(dfdv, s, pos0, vel0);

		dfdv *= dt; dfdx *= dt2;
		_sum = dfdv + dfdx;
		mat.diag[s.idx0] -= _sum;				
		mat.diag[s.idx1] -= _sum;				
		mat.value[matIdx[i].first] += _sum;		
		mat.value[matIdx[i].second] += _sum;	
	}
}

void ImplicitEuler_Integrator::makeSystem_Matrix(real dt, const QuadDynamicMesh<real>& mesh) {
	mat.zero();

	Matrix3<real> mTmp;
	for(size_t i=0;i<mat.n;++i)
		mat.diag[i] = mesh.mass * mTmp;

	updateSystem_Matrix(mat, mesh.pos0, mesh.vel0, dt, mesh.stretchSpring, mesh.stretchFunc, stretch_matIdx);
	updateSystem_Matrix(mat, mesh.pos0, mesh.vel0, dt, mesh.shearSpring, mesh.shearFunc, shear_matIdx);
	updateSystem_Matrix(mat, mesh.pos0, mesh.vel0, dt, mesh.bendingSpring, mesh.bendingFunc, bending_matIdx);
}

void ImplicitEuler_Integrator::makeSystem_Matrix(real dt, const TriDynamicMesh<real>& mesh) {
	mat.zero();

	Matrix3<real> mTmp;
	for(size_t i=0;i<mat.n;++i)
		mat.diag[i] = mesh.mass[i] * mTmp;

	updateSystem_Matrix(mat, mesh.pos0, mesh.vel0, dt, mesh.stretchSpring, mesh.stretchFunc, stretch_matIdx);
}

void updateSystem_RHS(Vec3Array<real>& b, const Vec3Array<real> &pos0, const Vec3Array<real> &vel0, real dt, 
					 const Array<Spring<real>>& spring, const SpringFunctor<real>* func) {
	real dt2 = dt * dt;
	Matrix3<real> dfdx;
	Vector3<real> dfdx_vel0, dfdx_vel1;

	for(size_t i=0;i<spring.size();++i) {
		const Spring<real>& s = spring[i];

		func->DfDx(dfdx, s, pos0, vel0);
		dfdx *= dt2;
		dfdx_vel0 = dfdx * vel0[s.idx0];
		dfdx_vel1 = dfdx * vel0[s.idx1];

		b[s.idx0] += dfdx_vel0;
		b[s.idx1] += dfdx_vel1;
		b[s.idx1] -= dfdx_vel0;
		b[s.idx0] -= dfdx_vel1;
	}	
}

void ImplicitEuler_Integrator::makeSystem_RHS(real dt, const QuadDynamicMesh<real>& mesh) {
	for(size_t i=0;i<b.size();++i) 
		b[i] = dt * f[i];

	updateSystem_RHS(b, mesh.pos0, mesh.vel0, dt, mesh.stretchSpring, mesh.stretchFunc);
	updateSystem_RHS(b, mesh.pos0, mesh.vel0, dt, mesh.shearSpring, mesh.shearFunc);
	updateSystem_RHS(b, mesh.pos0, mesh.vel0, dt, mesh.bendingSpring, mesh.bendingFunc);
}

void ImplicitEuler_Integrator::makeSystem_RHS(real dt, const TriDynamicMesh<real>& mesh) {
	for(size_t i=0;i<b.size();++i) 
		b[i] = dt * f[i];

	updateSystem_RHS(b, mesh.pos0, mesh.vel0, dt, mesh.stretchSpring, mesh.stretchFunc);
}