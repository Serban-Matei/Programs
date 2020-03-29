#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

const int N = 1; //Nは微分方程式のベクトル du/dt = f(t, u) の行数　
				 //ここを変更した場合class SolverのFuncも変更する必要がある

//Nが変更された場合この中の一部分も変更する必要がある
class Solvers {
private:

	int n_ini;												//微分方程式の離散個数
	double t_start, t_end;									//微分方程式の解く範囲　初期値は0と1
	double h_start;											//微分方程式のステップサイズ　初期値は0.1
	double delta;											//Adaptive stepsize における誤差の許容範囲 初期値は0.001
	double h_min, h_max;									//Adaptive stepsize における最小と最大のh 初期値は1と0.00000001
	vector<double> u_ini;									//微分方程式の初期値　初期値はすべて1
	int Flag_cout, Flag_write_xt, Flag_write_parameter_t;	//その他のフラグを定義している

public:

	//解きたい微分方程式たちの関数　N個必要 Nを変えた場合ここにも変更を加える必要がある
	//double Func1(const double t, const vector<double> u) {
	//	return -10 * u[0] + 10 * u[1];
	//}
	//double Func2(const double t, const vector<double> u) {
	//	return -u[0] * u[2] + 28 * u[0] - u[1];
	//}
	//double Func3(const double t, const vector<double> u) {
	//	return u[0] * u[1] - 2.666 * u[2];
	//}
	double Func1(const double t, const vector<double> u) {
		return -1 * u[0] + 2 * cos(t) * exp(-t);
	}

	//ここで初期化
	Solvers();

	//以下微分方程式を解くときの初期値などを変更する関数を宣言
	void set_t(double, double);
	void set_h(double);
	void set_n();
	void set_cout();
	void unset_cout();
	void set_delta(double);
	void set_h_min(double);
	void set_h_max(double);
	void set_write_xt();
	void unset_write_xt();
	void set_write_parameter();
	void unset_write_parameter();

	//よく使う機能を関数としてプロトタイプ宣言
	void Output_xt(const vector<double> &, const vector<vector<double>> &, string); //u[] - t関数としてのデータファイルをN個書き出す
	void Output_write_parameter_t(const vector<vector<double>> &, string); //ｔを曲線uのパラメータとして一つのファイルに出力　
	void Output_h(const vector<double> &, const vector<double> &, string s);
	double L_p_norm(const vector<double> &, int p); //Lpノルムを計算 p = 0　で　L∞

	//以下実際の数値計算の関数をプロトタイプ宣言している
	void Euler();
	void ERK4();
	void ERK23();
	void ERK45();
	void ERK23_O();
	void ERK45_O();

};
