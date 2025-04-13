#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Bridson03_Integrator.h>
#include <Nature/MassSpring/Core/ConstrainedPCG.h>

using namespace std;
using namespace springspace;

void Bridson03_Integrator::velocityUpdate_Implicit(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>& mesh) {
	//velocityUpdate_Implicit_(dt, x, mesh); return;

	Vector3<real> fext(0,g,0);
	for(size_t i=0;i<f.size();++i)
		f[i] = fext * mesh.mass[i];

	if(mesh.stretchSpring.size() != 0) SpringLibrary::addSpringElasticForce(f, mesh.stretchSpring, mesh.stretchFunc, x);

	makeSystem_RHS(dt, mesh);

	real residual;
	size_t nIter;
	makeDiagonalPreconditioner(mesh, dt);
	pcgSolve(dv, b, mesh, dt, residual, nIter);

	for(size_t i=0;i<mesh.vel.size();++i) 
		mesh.vel[i] = mesh.vel0[i] + dv[i];

	//cout << "Iter : " << nIter << " _ Residual : " << residual << endl;
}

void Bridson03_Integrator::velocityUpdate_Implicit_(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>& mesh) {
	Vector3<real> fext(0,g,0);
	for(size_t i=0;i<f.size();++i)
		f[i] = fext * mesh.mass[i];

	if(mesh.stretchSpring.size() != 0) SpringLibrary::addSpringElasticForce(f, mesh.stretchSpring, mesh.stretchFunc, x);

	makeSystem(dt, mesh);

	ConstrainedPCG<real> pcgsolver(mat, mesh.constraint);
	real residual;
	size_t nIter;

	pcgsolver.pcgSolve(dv, b, residual, nIter);

	for(size_t i=0;i<mesh.vel.size();++i) 
		mesh.vel[i] = mesh.vel0[i] + dv[i];

	cout << "Iter : " << nIter << " _ Residual : " << residual << endl;
}

void Bridson03_Integrator::allocateSystem(const TriDynamicMesh<real>& mesh) {
	mesh.allocate_SystemMatrix(mat);
	mat.make_MatrixIndexMapping(stretch_matIdx, mesh.stretchSpring);

	b.resize(mesh.size());
	dv.resize(mesh.size());
	f.resize(mesh.size());

	invDiag.resize(mesh.size());
	z.resize(mesh.size());
	s.resize(mesh.size());
	r.resize(mesh.size());
	tolerance = 1.0e-5;
	maxIter = 100;
}

void Bridson03_Integrator::makeSystem(real dt, const TriDynamicMesh<real>& mesh) {
	if(mat.n == 0) {
		basis::throwError("Error in ImplicitEuler_Integrator::make_SystemMatrix() _ Should allocate System First");
		return;
	}

	makeSystem_Matrix(dt, mesh);
	makeSystem_RHS(dt, mesh);
}

void Bridson03_Integrator::updateSystem_Matrix(MassSpringSystemMatrix<real>& mat, const Vec3Array<real> &pos0, const Vec3Array<real> &vel0, real dt, 
											 const Array<Spring<real>>& spring, const SpringFunctor<real>* func, const Array<std::pair<size_t, size_t>>& matIdx) {
	 real dt2 = dt*dt;
	 Matrix3<real> dfdv, _sum;

	 for(size_t i=0;i<spring.size();++i) {
		 const Spring<real>& s = spring[i];	

		 func->DfDv(dfdv, s, pos0, vel0);

		 _sum = dfdv * dt;
		 mat.diag[s.idx0] -= _sum;				
		 mat.diag[s.idx1] -= _sum;				
		 mat.value[matIdx[i].first] += _sum;		
		 mat.value[matIdx[i].second] += _sum;	
	 }
}

void Bridson03_Integrator::makeSystem_Matrix(real dt, const TriDynamicMesh<real>& mesh) {
	mat.zero();

	Matrix3<real> mTmp;
	for(size_t i=0;i<mat.n;++i)
		mat.diag[i] = mesh.mass[i] * mTmp;

	updateSystem_Matrix(mat, mesh.pos0, mesh.vel0, dt, mesh.stretchSpring, mesh.stretchFunc, stretch_matIdx);

	cout << mat.diag[0] << endl;
}

void Bridson03_Integrator::makeSystem_RHS(real dt, const TriDynamicMesh<real>& mesh) {
	for(size_t i=0;i<b.size();++i) 
		b[i] = dt * f[i];
}