//ver 1.2 7/18/2019
#pragma once
#include "pch.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "Matrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif

int FLAG_DISPLAY_OMP = 0;

//最大値,最小値を返す
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

//絶対値を返す
#define ABS(a) ((a) >= 0 ? (a) : -1 * (a))

using namespace std;

//初期化
Matrix::Matrix(vector<vector<double>>& A) {
	set_size(A);
	Matrix_band = Matrix_size;
	MAX_try = 1000;
	Omega = 1.;
	Error = 0.01;
	Norm_of_b = 0;


	if (FLAG_DISPLAY_OMP == 0) {
		#ifdef _OPENMP
			cout << "OPENMP is ON\n";
		#endif
		#ifndef _OPENMP
			cout << "OPENMP is OFF\n";
		#endif
		FLAG_DISPLAY_OMP++;
	}
}

//初期化
Matrix::Matrix(vector<vector<double>>& A, int i_MAX_try, double i_Error, double i_Omega) {
	set_size(A);
	Matrix_band = Matrix_size;
	MAX_try = i_MAX_try;
	Omega = i_Omega;
	Error = i_Error;
	Norm_of_b = 0; //for CG

	if (FLAG_DISPLAY_OMP == 0) {
		#ifdef _OPENMP
			cout << "OPENMP is ON\n";
		#endif
		#ifndef _OPENMP
			cout << "OPENMP is OFF\n";
		#endif
		FLAG_DISPLAY_OMP++;
	}
}

//行列のサイズを求める　もし正方行列でないならエラーを出す
int Matrix::set_size(vector<vector<double>>& A) {

	try {
		if (A.size() != A.front().size()) throw "Error:Mtrix is not square";
	}
	catch (char* errstr) {
		cout << errstr << "\n";
		exit(1);
	}

	Matrix_size = A.size();

	cout << "Matrix_size is " << Matrix_size << "\n" ;

}

//行列のバンド幅を求める
void Matrix::Matrix_Band(vector<vector<double>>& A) {
 
#ifdef _OPENMP
	//各スレッドでの最大のバンド幅ほ保存するための列
	vector<int> TEMP(omp_get_max_threads(), 0);
	int S = 0;

	//行列のバンドサイズを求める
	#pragma omp parallel for firstprivate(S)
	for (int i = 0; i < A.size(); i++) {
		for (int j = i; j < A.front().size(); j++) {
			if ((A[i][j] != 0) && (j - i > S)) {
				S = j - i;
			}
		}
		TEMP[omp_get_thread_num()] = S;
	}

	//各スレッドでの最大のバンド幅ほ保存するための列
	vector<int> TEMP2(omp_get_max_threads(), 0);
	S = 0;

	//行列のバンドサイズを求める
	#pragma omp parallel for firstprivate(S)
	for (int i = 0; i < A.front().size(); i++) {
		for (int j = i; j < A.size(); j++) {
			if ((A[j][i] != 0) && (i - j > S)) {
				S = i - j;
			}
		}
		TEMP2[omp_get_thread_num()] = S;
	}

	//各スレッドで求めたバンドの最大値が行列のバンド幅となる
	Matrix_band = MAX(*max_element(TEMP.begin(), TEMP.end()), *max_element(TEMP2.begin(), TEMP2.end()));

	cout << "Band size is " << Matrix_band << "\n";
#endif

#ifndef _OPENMP
	int S = 0;

	//行列のバンドサイズを求める
	for (int i = 0; i < A.size(); i++) {
		for (int j = i; j < A.front().size(); j++) {
			if ((A[i][j] != 0) && (j - i > S)) {
				S = j - i;
			}
		}
	}

	//各スレッドでの最大のバンド幅ほ保存するための列
	int S2 = 0;

	//行列のバンドサイズを求める
	for (int i = 0; i < A.front().size(); i++) {
		for (int j = i; j < A.size(); j++) {
			if ((A[j][i] != 0) && (i - j > S)) {
				S2 = i - j;
			}
		}
	}

	Matrix_band = max(S, S2);

	cout << " Band size with no omp is " << Matrix_band << "\n";
#endif
}

//b - Auのノルムを計算する
double Matrix::Error_check(const vector<vector<double>>& A, const vector<double>& u, const vector<double>& b) {
	
	vector<double> v(A.size(), 0);
	int p = 2;
	double tmp = 0, norm = 0;

	//v = Au
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < A.size(); i++) {
			for (int j = MAX(i - Matrix_band, 0); j < MIN(i + Matrix_band + 1, Matrix_size); j++) {
				v[i] += A[i][j] * u[j];
			}
		}

		//v = v - b
		#pragma omp for
		for (int i = 0; i < A.size(); i++) {
			v[i] -= b[i];
		}

		//vのノルムを求める
		#pragma omp for private(tmp)
		for (int i = 0; i < A.size(); i++) {
			tmp = abs(v[i]);
			for (int j = 1; j < p; j++) {
				tmp *= abs(v[i]);
			}
			#pragma omp atomic
			norm += tmp;
		}
	}
	norm = pow(norm, 1. / (double)p);

	//bのノルムを求める
	if (Norm_of_b == 0) {
		tmp = 0;
		for (int i = 0; i < A.size(); i++) {
			tmp = abs(b[i]);
			for (int j = 1; j < p; j++) {
				tmp *= abs(b[i]);
			}
			Norm_of_b += tmp;
		}
		Norm_of_b = pow(Norm_of_b, 1. / (double)p);
	}

	return norm / Norm_of_b;
}

//u - u_1のノルムを計算する
double Matrix::Error_check_C(const vector<double> &u, const vector<double> &u_1) {

	vector<double> v(u.size(), 0);
	int p = 2;
	double tmp = 0, norm = 0;

	//#pragma omp parallel
	{
		//#pragma omp for
		for (int i = 0; i < v.size(); i++) {
			v[i] = u[i] - u_1[i];
		}

		//#pragma omp for private(tmp)
		for (int i = 0; i < v.size(); i++) {
			tmp = abs(v[i]);
			for (int j = 1; j < p; j++) {
				tmp *= abs(v[i]);
			}
		//#pragma omp atomic
			norm += tmp;
		}
	}
	norm = pow(norm, 1. / (double)p);
	return norm;
}

//ベクトルのノルムの2乗を計算する (parallel) 
inline double Matrix::Norm_2(const vector<double> &r) {

	double tmp = 0, norm = 0;

	#pragma omp parallel
	{

		#pragma omp for private(tmp) reduction(+:norm)
		for (int i = 0; i < r.size(); i++) {
			tmp = abs(r[i]);
			tmp *= abs(r[i]);
			norm += tmp;
		}
	}
	return norm;

}

int Matrix::LU(vector<vector<double>>& A, vector<double>& u, vector<double>& b){
	
	double m = 0;

	vector<int> P(Matrix_size, 0);
	for (int i = 0; i < Matrix_size; i++) {
		P[i] = i; //初期化
	}

	for (int k = 0; k < Matrix_size; k++) {

		//Pivot
		int Pivot = k;

		for (int i = k + 1; i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
			if (ABS(A[i][k]) > ABS(A[Pivot][k])) {
				Pivot = i;
				cout << "Pivot";
			}
		}

		if (A[Pivot][k] == 0) {
			cout << "singular matrix";
			return 1;
		}

		if (Pivot != k) {
			cout << "pivot";
			int temp = P[k];
			P[k] = P[Pivot];
			P[Pivot] = temp;
			double tempd;
			for (int j = 0; j < Matrix_size; j++) {
				tempd = A[k][j];
				A[k][j] = A[Pivot][j];
				A[Pivot][j] = tempd;
			}
		}
		
		//以下LU分解		
		m = 1 / A[k][k];
		#pragma omp parallel
		{		
			#pragma omp for
			for (int i = k + 1; i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
				A[i][k] = m * A[i][k];
				for (int j = k + 1; j < MIN(k + Matrix_band + 1, Matrix_size); j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
		cout << "\r" << Matrix_size << " : " << k + 1;
	}

	double TEMP;

	//前進代入
	for (int i = 0; i < Matrix_size; i++) {
		TEMP = 0;
		for (int k = 0; k < i; k++) {
			TEMP += A[i][k] * u[k];
		}
		u[i] = b[P[i]] - TEMP;
	}

	//後進代入
	for (int i = Matrix_size - 1; i >= 0; i--) {
		TEMP = 0;
		for (int k = i + 1; k < Matrix_size; k++) {
			TEMP += A[i][k] * u[k];
		}
		u[i] = (u[i] - TEMP) / A[i][i];
	}

	return 0;
}

vector<int> Matrix::LU_D(vector<vector<double>> &A) {


	double m = 0;

	vector<int> P(Matrix_size, 0);
	for (int i = 0; i < Matrix_size; i++) {
		P[i] = i; //初期化
	}

	for (int k = 0; k < Matrix_size; k++) {

		//Pivot
		int Pivot = k;

		for (int i = k + 1; i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
			if (ABS(A[i][k]) > ABS(A[Pivot][k])) { //k を pivot に変更
				Pivot = i;
				cout << "Pivot";
			}
		}

		if (A[Pivot][k] == 0) {
			cout << "singular matrix\n";
		}

		if (Pivot != k) {
			cout << "a";
			int temp = P[k];
			P[k] = P[Pivot];
			P[Pivot] = temp;
			double tempd;
			for (int j = 0; j < Matrix_size; j++) {
				tempd = A[k][j];
				A[k][j] = A[Pivot][j];
				A[Pivot][j] = tempd;
			}
		}

		//以下LU分解		
		m = 1 / A[k][k];
		#pragma omp parallel
		{
			#pragma omp for
			for (int i = k + 1; i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
				A[i][k] = m * A[i][k];
				for (int j = k + 1; j < MIN(k + Matrix_band + 1, Matrix_size); j++) {
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
		cout << "\r" << Matrix_size << " : " << k + 1;
	}

	return P;

}

int Matrix::LU_SUB(const vector<vector<double>>& A, vector<double>& u, vector<double>& b, const vector<int>& P) {

	double TEMP;

	//前進代入
	for (int i = 0; i < Matrix_size; i++) {
		TEMP = 0;
		#pragma omp parallel for
		for (int k = 0; k < i; k++) {
			TEMP += A[i][k] * u[k];
		}
		u[i] = b[P[i]] - TEMP;
	}

	//後進代入
	for (int i = Matrix_size - 1; i >= 0; i--) {
		TEMP = 0;
		for (int k = i + 1; k < Matrix_size; k++) {
			TEMP += A[i][k] * u[k];
		}
		u[i] = (u[i] - TEMP) / A[i][i];
	}

	return 0;

}

int Matrix::SOR(vector<vector<double>>& A, vector<double>& u, vector<double>& b) {

	cout << "Pivoting";
	int ABC = Matrix_size / 20;

	//Pivot
	for (int k = 0; k < Matrix_size; k++) {
		int Pivot = k;
		for (int i = 0; i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
			if (abs(A[i][k]) > abs(A[k][k])) {
				Pivot = i;
			}
		}

		if (A[Pivot][k] == 0) {
			cout << "singular matrix";
			return 1;
		}

		if (Pivot != k) {
			double tempd;
			tempd = b[k];
			b[k] = b[Pivot];
			b[Pivot] = tempd;
			for (int j = 0; j < Matrix_size; j++) {
				tempd = A[k][j];
				A[k][j] = A[Pivot][j];
				A[Pivot][j] = tempd;
			}
		}

		double m = 1 / A[k][k];
		for (int i = MAX(k - Matrix_band, 0); i < MIN(k + Matrix_band + 1, Matrix_size); i++) {
			A[k][i] = A[k][i] * m;
		}
		b[k] = b[k] * m;

		if (k % ABC == 0) {
			cout << ".";
		}

	}

	cout << "done!\n";

	//SOR
	double tmp;
	double  SOR_Error = 0;
	vector<double> u_1(Matrix_size, 0);

	for (int i = 0; i < MAX_try; i++) {
		for (int r = 0; r < Matrix_size; r++) {
			tmp = 0.;
			for (int j = MAX(r - Matrix_band, 0); j < MIN(r + Matrix_band + 1, Matrix_size); j++) {
				tmp += A[r][j] * u[j];
			}
			tmp -= b[r];

			u_1[r] = u[r];
			u[r] = u[r] - Omega * tmp;
		}
		cout << "\rmax iteration:current step " << MAX_try << " : " << i << "   error:current error "
			 << Error << "  " << SOR_Error;

		//Check if the norm is less than epsilon 
		//if (i % 1 == 0) {
			//SOR_Error = Error_check(A, u, b);
			SOR_Error = Error_check_C(u, u_1);
			if (SOR_Error < Error) {
				cout << "\n" << "Number of trials was " << i << " times\n";
				return 0;
			}
		//}
	}

	cout << "\nThe computation ended because the number of trials hit the maximum tirials number of " << MAX_try << " times\n";
	return 1;

}

//int Matrix::CG(vector<vector<double>>& A, vector<double>& u, vector<double>& b) {
//
//	vector<double> r = b; //残差ベクトル
//	vector<double> d = b; //勾配ベクトル
//	vector<double> v(Matrix_size, 0); //計算に使用
//	double r_norm, r_old_norm, beta, Omega_CG, ERROR_;
//	double norm = 0;
//
//	r_norm = Norm_2(r);
//
//	if (pow(r_norm, 0.5) < Error) {
//		return 0;
//	}
//
//	for (int k = 0; k < MAX_try; k++) {
//		//step 1
//		double product = 0;
//		v = vector<double>(Matrix_size, 0);
//
//		#pragma omp parallel
//		{
//			#pragma omp for
//			for (int i = 0; i < Matrix_size; i++) {
//				for (int j = MAX(i - Matrix_band, 0); j < MIN(i + Matrix_band + 1, Matrix_size); j++) {
//					v[i] += A[i][j] * d[j];
//				}
//			}
//
//			#pragma omp for reduction(+: product)
//			for (int i = 0; i < Matrix_size; i++) {
//				product += d[i] * v[i];
//			}
//
//			#pragma omp single
//			{
//				Omega_CG = r_norm / product;
//			}
//
//			//step 2
//			#pragma omp for nowait
//			for (int i = 0; i < Matrix_size; i++) {
//				u[i] = u[i] + Omega_CG * d[i];
//			}
//
//			#pragma omp for
//			for (int i = 0; i < Matrix_size; i++) {
//				r[i] = r[i] - Omega_CG * v[i];
//			}
//		}
//
//		ERROR_ = pow(r_norm, 0.5);
//		if (ERROR_ < Error) {
//			return 0;
//		}
//
//		//step 3
//		r_old_norm = r_norm;
//		norm = 0;
//
//		#pragma omp parallel
//		{
//			#pragma omp for reduction(+:norm)
//			for (int i = 0; i < r.size(); i++) {
//				norm += abs(r[i]) * abs(r[i]);
//			}
//
//			#pragma omp single
//			{
//				r_norm = norm;
//				beta = r_norm / r_old_norm;
//			}
//
//			//step 4
//			#pragma omp for
//			for (int i = 0; i < Matrix_size; i++) {
//				d[i] = r[i] + beta * d[i];
//			}
//		}
//		
//		cout << "\r" << Error << ":" << ERROR_;
//
//	}
//
//}

int Matrix::CG(vector<vector<double>>& A, vector<double>& u, vector<double>& b) {

	vector<double> r = b; //残差ベクトル
	vector<double> d = b; //勾配ベクトル
	vector<double> v(Matrix_size, 0); //計算に使用
	double beta, Omega_CG;
	double product_dAd = 0;
	double inner_product_rr = 0;
	double inner_product_rr_old = 0;

	for (int i = 0; i < Matrix_size; i++) {
		inner_product_rr += r[i] * r[i];
	}

	if (inner_product_rr < Error * Error) {
		return 0;
	}

	//番号はCG法実用版を参照
	for (int k = 0; k < MAX_try; k++) {

		product = 0;
		v = vector<double>(Matrix_size, 0);

		//1, 2を計算
		#pragma omp parallel
		{
			#pragma omp for reduction(+: product_dAd)
			for (int i = 0; i < Matrix_size; i++) {
				for (int j = MAX(i - Matrix_band, 0); j < MIN(i + Matrix_band + 1, Matrix_size); j++) {
					v[i] += A[i][j] * d[j];
				}
				product_dAd += d[i] * v[i];
			}

			#pragma omp single
			{
				Omega_CG = inner_product_rr / product_dAd;
				inner_product_rr_old = inner_product_rr;
			}

			//3
			#pragma omp for reduction(+: inner_product_rr)
			for (int i = 0; i < Matrix_size; i++) {
				u[i] = u[i] + Omega_CG * d[i];
				r[i] = r[i] - Omega_CG * v[i];
				inner_product_rr += r[i] * r[i];
			}

			#pragma omp single
			{
				if (inner_product_rr < Error * Error) {
					return 0;
				}
			}

			//4
			#pragma omp single
			{
				beta = inner_product_rr / inner_product_rr_old;
			}

			//5
			#pragma omp for
			for (int i = 0; i < Matrix_size; i++) {
				d[i] = r[i] + beta * d[i];
			}
		}

		cout << "\r" << Error << ":" << inner_product_rr;

	}

}