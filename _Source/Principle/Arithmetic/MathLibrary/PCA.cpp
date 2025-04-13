#include <stdafx.h>
#include <Principle/Arithmetic/MathLibrary/PCA.h>

using namespace std;
using namespace basis;
using namespace arithmetic;

namespace arithmetic {

void testPCA() {

	int dimPCA = 8, numData = 6;

	PCA<double> pca(dimPCA);

	std::vector<VectorN<double>> basis(dimPCA);
	std::vector<double> limits(dimPCA);
	std::vector<VectorN<double>> coeff(numData);
	VectorN<double> vec0(dimPCA), vec1(dimPCA);	

	for(int i=0;i<dimPCA;++i) {
		basis[i].resize(dimPCA);
		basis[i][i] = 1.0f;			// set unit vector		
	}
	for(int i=0;i<numData;++i)
		coeff[i].resize(dimPCA);

	// Set Basis
//	basis[0][0] = 0.25088;  basis[0][1] = -0.551301;  basis[0][2] = -0.698255; basis[0][3] = 0.381531;
//	basis[1][0] = 0.389986; basis[1][1] = 0.474038;   basis[1][2] = -0.545798; basis[1][3] = -0.570354;
//	basis[2][0] = 0.530362; basis[2][1] = -0.579799;  basis[2][2] = 0.382941;  basis[2][3] = -0.485701;
//	basis[3][0] = 0.709709; basis[3][1] = 0.36768;	basis[3][2] = 0.260577;  basis[3][3] = 0.541502;

	// Set Data Range
	limits[0] = 100.0f;
	limits[1] = 1.0f;
	limits[2] = 30.0f;
	limits[3] = 0.01f;
	for(int i=4;i<dimPCA;++i) limits[i] = i;	

	// Set Data Average
	MatrixN<double> dataMatrix(dimPCA, numData);
	VectorN<double> avg(dimPCA);
	avg.zero();		avg(0) = 20.0f; avg(1) = 80.0f; avg(2) = 1.0f; avg(3) = 5.0f; for(int i=4;i<dimPCA;++i) avg[i] = i;

	// Set Data Randomly
	cout << endl << "***** Set Data *****" << endl;
	vec0.resize(dimPCA); vec1.resize(dimPCA);
	for(int i=0;i<numData;++i) {
		vec0 = avg;
		for(int j=0;j<dimPCA;++j) {
			vec1 = basis[j] * (limits[j] * (2.0f * ((double)rand() / (double)RAND_MAX) - 1.0f));
			vec0 += vec1;
		}
		pca.push(vec0);
		dataMatrix.setCol(i, vec0);
		cout << i << "th data : " << vec0 << endl;
	}

	// Perform PCA
	cout << endl << "***** PCA Results *****" << endl;
	pca.doPCA();
	cout << pca;

	//// Modify DataMatrix _ Average = 0
	//VectorN<double> data(dimPCA);
	//VectorN<double> average(dimPCA);
	//pca.getAverage(average);
	//for(int i=0;i<numData;++i) {
	//	dataMatrix.getCol(i, data);
	//	data = data - average;
	//	dataMatrix.setCol(i, data);
	//}

	//// Check PCA Basis
	//MatrixN<double> basisMatrix(dimPCA, dimPCA);
	//cout << endl << "***** Check PCA Basis *****" << endl;
	//for(int i=0;i<dimPCA;++i) {
	//	pca.getBasis(i, vec0);
	//	basisMatrix.setCol(i, vec0);
	//	cout << vec0.mag() << " ";
	//	for(int j=i+1;j<dimPCA;++j) {
	//		pca.getBasis(j, vec1);
	//		cout << dot(vec0, vec1) << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	//// Get Deviation Matrix
	//MatrixN<double> deviationMatrix(dimPCA, dimPCA);
	//for(int i=0;i<dimPCA;++i) 
	//	deviationMatrix.set(i,i,pca.getDeviation(i));

	//// Get Coefficient
	//cout << endl << "***** Data Coefficient *****" << endl;
	//for(int i=0;i<numData;++i) {
	//	pca.getCoefficient(i, coeff[i]);
	//	cout << i << "th : " << coeff[i] << endl;
	//}
	//cout << endl;

	//// Data Reconstruction
	//cout << endl << "***** Data Reconstruction *****" << endl;
	//std::vector<VectorN<double>> pcaBasis(dimPCA);
	//for(int i=0;i<dimPCA;++i) { pcaBasis[i].resize(dimPCA); pca.getBasis(i, pcaBasis[i]); }

	//for(int i=0;i<numData;++i) {
	//	pca.getAverage(vec0);
	//	for(int j=0;j<dimPCA;++j) {			
	//		vec1 = pcaBasis[j] * coeff[i][j];
	//		vec0 += vec1;
	//	}

	//	cout << i << "th data : " << vec0 << endl;
	//}

	//// Make V^t for UDV^t
	//cout << endl << "***** Get Vt columns *****" << endl;
	//std::vector<VectorN<double>> columnVt(numData);
	//std::vector<VectorN<double>> rowVt(dimPCA);
	//for(int i=0;i<numData;++i) columnVt[i].resize(dimPCA);
	//for(int i=0;i<dimPCA;++i) rowVt[i].resize(numData);

	//for(int i=0;i<numData;++i) {
	//	for(int j=0;j<dimPCA;++j) {
	//		if(pca.getDeviation(j) > 1.e-4) {
	//			columnVt[i][j] = coeff[i][j] / pca.getDeviation(j);
	//			rowVt[j][i] = columnVt[i][j];
	//		}
	//		else {
	//			columnVt[i][j] = coeff[i][j];
	//			rowVt[j][i] = columnVt[i][j];
	//		}
	//	}
	//	cout << i << "th column : " << columnVt[i] << endl;
	//}

	//// Check Vt Column
	//cout << endl << "***** Check Vt *****" << endl;
	//for(int i=0;i<numData;++i) {
	//	//cout << columnVt[i].mag() << " ";
	//	cout << rowVt[i].mag() << " ";
	//	for(int j=i+1;j<numData;++j) {			
	//		//cout << dot(columnVt[i], columnVt[j]) << " ";
	//		cout << dot(rowVt[i], rowVt[j]) << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	//for(int i=0;i<numData;++i) {
	//	mul(vec0, dataMatrix, rowVt[i]);
	//	//mul(vec0, dataMatrix, coeff[i]);
	//	vec1 = pcaBasis[i] * pca.getDeviation(i);
	//	cout << vec0 << endl;
	//	cout << vec1 << endl;
	//}

	//cout << endl << endl;
	//cout << dataMatrix << endl;
	//cout << basisMatrix << endl;
	//cout << deviationMatrix << endl;

	//MatrixN<double> coeffMatrix(dimPCA, dimPCA), coeffMatrixTranspose(dimPCA, dimPCA);
	//MatrixN<double> basisMatrixTranspose(basisMatrix.dimCol(), basisMatrix.dimRow());
	//	
	//basisMatrixTranspose = transpose(basisMatrix);
	//coeffMatrix = basisMatrixTranspose * dataMatrix;
	//coeffMatrixTranspose = transpose(coeffMatrix);

	//cout << coeffMatrix << endl;
	//cout << coeffMatrixTranspose << endl;

	//MatrixN<double> tmp(dimPCA, dimPCA);
	//tmp = coeffMatrix * coeffMatrixTranspose;

	//tmp /= (double)(numData - 1);
	//cout << tmp << endl;

}

};		// namespace arithmetic