#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Matrix {

private:
	int Matrix_size; // = 0
	int Matrix_band; // = 0
	int MAX_try; // = 1000
	double Omega; // = 1
	double Error; //0.01
	double Norm_of_b; // = 0

public:
	Matrix(vector<vector<double>>&);
	Matrix(vector<vector<double>>&, int, double, double);
	int set_size(vector<vector<double>>&);	
	void Matrix_Band(vector<vector<double>>&);
	double Error_check(const vector<vector<double>>&, const vector<double>&, const vector<double>&);
	double Error_check_C(const vector<double>&, const vector<double>&);
	double Norm_2(const vector<double>&);

	int LU(vector<vector<double>>&, vector<double>&, vector<double>&);
	int SOR(vector<vector<double>>&, vector<double>&, vector<double>&);
	int CG(vector<vector<double>>&, vector<double>&, vector<double>&);

	vector<int> LU_D(vector<vector<double>>&);
	int LU_SUB(const vector<vector<double>>&, vector<double>&, vector<double>&, const vector<int>&);

};
