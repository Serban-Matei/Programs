#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "ode.h"

//ここで初期化
Solvers::Solvers() {
	t_start = 0.;
	t_end = 1.;
	h_start = 0.1;
	u_ini = vector<double>(N, 1.);
	n_ini = 10;
	delta = 0.001;
	h_min = 0.00000001;
	h_max = 1;
	Flag_cout = 0;
	Flag_write_xt = 0;
	Flag_write_parameter_t = 0;
}

//以下微分方程式を解くときの初期値などを変更する関数を定義
void Solvers::set_t(double t1, double t2) {
	t_start = t1;
	t_end = t2;
	set_n();
}
void Solvers::set_h(double h1) {
	h_start = h1;
	set_n();
}
void Solvers::set_n() {
	double Temp = t_end - t_start;
	n_ini = Temp / h_start;
}
void Solvers::set_cout() {
	Flag_cout = 1;
}
void Solvers::unset_cout() {
	Flag_cout = 0;
}
void Solvers::set_delta(double a) {
	delta = a;
}
void Solvers::set_h_min(double h) {
	h_min = h;
}
void Solvers::set_h_max(double h) {
	h_max = h;
}
void Solvers::set_write_xt() {
	Flag_write_xt = 1;
}
void Solvers::unset_write_xt() {
	Flag_write_xt = 0;
}
void Solvers::set_write_parameter() {
	Flag_write_parameter_t = 1;
}
void Solvers::unset_write_parameter() {
	Flag_write_parameter_t = 0;
}

//関数を成分に持つ配列を定義している Nを変えた場合ここに変更が必要になる
double (Solvers::*Func[])(double, vector<double>) {

	&Solvers::Func1,

};

//便利な関数を定義
void Solvers::Output_xt(const vector<double> &t, const vector<vector<double>> &u, string s) {

	string dat(".dat");
	string Int;

	//書き出し　
	char Enter = '\n';
	char Space = ' ';

	for (int j = 0; j < N; j++) {

		string name;
		Int = to_string(j);
		name += s + Int + dat;

		ofstream outputfile(name);

		for (int i = 0; i < n_ini; i++) {

			outputfile << t[i];
			outputfile << Space;
			outputfile << u[i][j];
			outputfile << Space;
			outputfile << Enter;

		}
		outputfile.close();
		outputfile.clear();
	}

}

void Solvers::Output_write_parameter_t(const vector<vector<double>> &u, string s) {

	string dat(".dat");

	//書き出し　
	char Enter = '\n';
	char Space = ' ';

	string name;
	name += s + dat;

	ofstream outputfile(name);

	for (int i = 0; i < n_ini; i++) {
		for (int j = 0; j < N; j++) {

			outputfile << u[i][j];
			outputfile << Space;

		}
		outputfile << Enter;
	}

	outputfile.close();
	outputfile.clear();

}

void Solvers::Output_h(const vector<double> &t, const vector<double> &h, string s) {

	string dat(".dat");

	//書き出し　
	char Enter = '\n';
	char Space = ' ';

	string name;
	name += s + dat;

	ofstream outputfile(name);

	for (int i = 0; i < n_ini; i++) {

		outputfile << t[i];
		outputfile << Space;
		outputfile << 100 * h[i]; //グラフにしたときわかりやすいように100倍にしてある
		outputfile << Enter;

	}

	outputfile.close();
	outputfile.clear();

}

double Solvers::L_p_norm(const vector<double> &u, int p) {

	double Norm = 0;
	double Temp = 0;

	if (p > 0) {
		for (int i = 0; i < N; i++) {

			Temp = abs(u[i]);

			for (int j = 1; j < p; j++) {
				Temp *= abs(u[i]);
			}

			Norm += Temp;
		}
		Norm = pow(Norm, 1. / (double)p);
	}

	if (p == 0) {

		Temp = u[0];

		for (int i = 0; i < N; i++) {

			if (u[i] > Temp) {

				Temp = u[i];

			}
		}
		Norm = Temp;
	}

	return Norm;

}


//以下実際の数値計算の関数を定義している

//Euler法を定義する
void Solvers::Euler() {

	vector<vector<double>> u; //各ステップで関数uを計算、保存する
	vector<double> t; //パラメータtを保存

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	t = vector<double>(n_ini, 0);

	t[0] = t_start;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//以下オイラー法の計算
	for (int i = 1; i < n_ini; i++) {
		for (int j = 0; j < N; j++) {

			u[i][j] = u[i - 1][j] + h_start * (this->*Func[j])(t[i - 1], u[i - 1]);

			if (Flag_cout == 1) {
				cout << u[i][j] << "  ";
			}
		}

		t[i] = t[i - 1] + h_start;
		if (Flag_cout == 1) {
			cout << "\n";
		}

	}

	//Flag_write が1ならば結果を書き出す
	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_Euler");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_Euler_parameter");
	}

}

//4次のERKを定義する
void Solvers::ERK4() {

	vector<vector<double>> u, A; //AはRunge-Kutta matrix
	vector<double> t, k2, k3, k4, b, c; //b, c はそれぞれRunge-Kutta wight, Runge-Kutta node

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	A = vector<vector<double>>(4, vector<double>(4, 0));
	t = vector<double>(n_ini, 0);
	k2 = vector<double>(N, 0);
	k3 = vector<double>(N, 0);
	k4 = vector<double>(N, 0);
	b = vector<double>(4, 0);
	c = vector<double>(4, 0);

	t[0] = t_start;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//Rk matrix, wight, node を定義
	A[1][0] = 1. / 2.;
	A[2][1] = 1. / 2.;
	A[3][2] = 1.;

	c[1] = 1. / 2.;
	c[2] = 1. / 2.;
	c[3] = 1.;

	b[0] = 1. / 6.;
	b[1] = 1. / 3.;
	b[2] = 1. / 3.;
	b[3] = 1. / 6.;

	//以下ERK4法の計算
	for (int i = 1; i < n_ini; i++) {

		//kの計算
		for (int j = 0; j < N; j++) {

			k2[j] = u[i - 1][j] + h_start * A[1][0] * (this->*Func[j])(t[i - 1], u[i - 1]);

		}
		for (int j = 0; j < N; j++) {

			k3[j] = u[i - 1][j] + h_start * A[2][0] * (this->*Func[j])(t[i - 1], u[i - 1])
				+ h_start * A[2][1] * (this->*Func[j])(t[i - 1] + c[1] * h_start, k2);

		}
		for (int j = 0; j < N; j++) {

			k4[j] = u[i - 1][j] + h_start * A[3][0] * (this->*Func[j])(t[i - 1], u[i - 1])
				+ h_start * A[3][1] * (this->*Func[j])(t[i - 1] + c[1] * h_start, k2)
				+ h_start * A[3][2] * (this->*Func[j])(t[i - 1] + c[2] * h_start, k3);

		}

		//u[i]の計算
		for (int j = 0; j < N; j++) {

			u[i][j] = u[i - 1][j] + h_start * b[0] * (this->*Func[j])(t[i - 1], u[i - 1])
				+ h_start * b[1] * (this->*Func[j])(t[i - 1] + c[1] * h_start, k2)
				+ h_start * b[2] * (this->*Func[j])(t[i - 1] + c[2] * h_start, k3)
				+ h_start * b[3] * (this->*Func[j])(t[i - 1] + c[3] * h_start, k4);

			if (Flag_cout == 1) {
				cout << u[i][j] << "  ";
			}

		}

		t[i] = t[i - 1] + h_start;
		if (Flag_cout == 1) {
			cout << "\n";
		}

	}

	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_ERK4");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_ERK4_parameter");
	}

}

//Adaptive stepsize RKF23
void Solvers::ERK23() {

	vector<vector<double>> u, A; // AはRunge-Kutta matrix
	vector<double> t, x, k2, k3, b1, b2, c; //xはError control用のベクトル b, c はそれぞれRunge-Kutta wight, Runge-Kutta node
	vector<double> h; //adaptive stepsize h_adaptiveの記録用
	double K, h_adaptive = h_start; //可変h

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	A = vector<vector<double>>(3, vector<double>(3, 0));
	t = vector<double>(n_ini, 0);
	h = vector<double>(n_ini, 0);
	x = vector<double>(N, 0);
	k2 = vector<double>(N, 0);
	k3 = vector<double>(N, 0);
	b1 = vector<double>(2, 0);
	b2 = vector<double>(3, 0);
	c = vector<double>(3, 0);

	t[0] = t_start;
	h[0] = 0;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//Rk matrix, wight, node を定義
	A[1][0] = 2. / 3.;
	A[2][1] = 2. / 3.;

	c[1] = 2. / 3.;
	c[2] = 2. / 3.;

	b1[0] = 1. / 4.;
	b1[1] = 3. / 4.;

	b2[0] = 1. / 4.;
	b2[1] = 3. / 8.;
	b2[2] = 3. / 8.;

	cout << "ERK23\n\n";

	//以下ERK23法の計算
	for (int i = 1; i < n_ini; i++) {

		for (;;) {
			//kの計算
			for (int j = 0; j < N; j++) {

				k2[j] = u[i - 1][j] + h_adaptive * A[1][0] * (this->*Func[j])(t[i - 1], u[i - 1]);

			}
			for (int j = 0; j < N; j++) {

				k3[j] = u[i - 1][j] + h_adaptive * A[2][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[2][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

			}

			//u[i], x の計算
			for (int j = 0; j < N; j++) {

				u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

				x[j] = u[i - 1][j] + h_adaptive * b2[0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * b2[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * b2[2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3);
			}

			vector<double> Temp(N, 0);

			for (int j = 0; j < N; j++) {

				Temp[j] = u[i][j] - x[j];

			}

			//誤差ノルムを計算
			K = L_p_norm(Temp, 2);

			if ((K <= delta) && (h_adaptive < h_max)) {

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;

				if (h_adaptive * pow(abs(delta / K), 1. / 3.) < h_max) {
					h_adaptive = h_adaptive * pow(abs(delta / K), 1. / 3.);
				}

				break;

			}if (h_adaptive <= h_min) {

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				break;

			}

			//精度が満たされなかった場合stepsize hをきつくする
			h_adaptive = h_adaptive * pow(abs(delta / K), 1. / 3.);
		}

		if (t[i] >= t_end) {
			//tがt_endを超えた時点でループから抜ける
			for (int j = i; j < n_ini; j++) {
				t[j] = t[i];
			}
			break;
		}

		if (Flag_cout == 1) {
			for (int l = 0; l < N; l++) {
				cout << u[i][l];
				cout << " ";
			}
			cout << h_adaptive;
			cout << "\n";
		}
	}

	//ファイルに書き出し
	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_ERK23");
		Output_h(t, h, "Result_ERK23_stepsize");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_ERK23_parameter");
		Output_h(t, h, "Result_ERK23_stepsize");
	}

}

//Adaptive stepsize RKF45
void Solvers::ERK45() {

	vector<vector<double>> u, A; // AはRunge-Kutta matrix
	vector<double> t, x, k2, k3, k4, k5, k6, b1, b2, c; //xはError control用のベクトル b, c はそれぞれRunge-Kutta wight, Runge-Kutta node
	vector<double> h; //adaptive stepsize h_adaptiveの記録用
	double K, h_adaptive = h_start; //可変h

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	A = vector<vector<double>>(6, vector<double>(6, 0));
	t = vector<double>(n_ini, 0);
	h = vector<double>(n_ini, 0);
	x = vector<double>(N, 0);
	k2 = vector<double>(N, 0);
	k3 = vector<double>(N, 0);
	k4 = vector<double>(N, 0);
	k5 = vector<double>(N, 0);
	k6 = vector<double>(N, 0);
	b1 = vector<double>(5, 0);
	b2 = vector<double>(6, 0);
	c = vector<double>(6, 0);

	t[0] = t_start;
	h[0] = 0;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//Rk matrix, wight, node を定義
	A[1][0] = 1. / 4.;
	A[2][0] = 3. / 32.;
	A[3][0] = 1932. / 2197.;
	A[4][0] = 439. / 216.;
	A[5][0] = -8. / 27.;
	A[2][1] = 9. / 32.;
	A[3][1] = -7200. / 2197.;
	A[4][1] = -8.;
	A[5][1] = 2.;
	A[3][2] = 7296. / 2197.;
	A[4][2] = 3680. / 513.;
	A[5][2] = -3544. / 2565;
	A[4][3] = -845. / 4104;
	A[5][3] = 1859. / 4104.;
	A[5][4] = -11. / 40.;

	c[1] = 1. / 4.;
	c[2] = 3. / 8.;
	c[3] = 12. / 13.;
	c[4] = 1.;
	c[5] = 1. / 2.;

	b1[0] = 25. / 216.;
	b1[2] = 1408. / 2565.;
	b1[3] = 2197. / 4104.;
	b1[4] = -1. / 5.;

	b2[0] = 16. / 135.;
	b2[2] = 6656. / 12825.;
	b2[3] = 28561. / 56430.;
	b2[4] = -9. / 50;
	b2[5] = 2. / 55.;

	cout << "ERK45\n\n";

	//以下ERK45法の計算
	for (int i = 1; i < n_ini; i++) {

		for (;;) {
			//kの計算
			for (int j = 0; j < N; j++) {

				k2[j] = u[i - 1][j] + h_adaptive * A[1][0] * (this->*Func[j])(t[i - 1], u[i - 1]);

			}
			for (int j = 0; j < N; j++) {

				k3[j] = u[i - 1][j] + h_adaptive * A[2][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[2][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

			}
			for (int j = 0; j < N; j++) {

				k4[j] = u[i - 1][j] + h_adaptive * A[3][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[3][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[3][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3);

			}
			for (int j = 0; j < N; j++) {

				k5[j] = u[i - 1][j] + h_adaptive * A[4][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[4][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[4][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * A[4][3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4);

			}
			for (int j = 0; j < N; j++) {

				k6[j] = u[i - 1][j] + h_adaptive * A[5][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[5][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[5][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * A[5][3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
					+ h_adaptive * A[5][4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5);

			}

			//u[i], x の計算
			for (int j = 0; j < N; j++) {

				u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * b1[2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * b1[3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
					+ h_adaptive * b1[4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5);

				x[j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * b2[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * b2[2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * b2[3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
					+ h_adaptive * b2[4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5)
					+ h_adaptive * b2[5] * (this->*Func[j])(t[i - 1] + c[5] * h_adaptive, k6);

			}

			vector<double> Temp(N, 0);

			for (int j = 0; j < N; j++) {

				Temp[j] = u[i][j] - x[j];

			}

			//誤差ノルムを計算
			K = L_p_norm(Temp, 2);

			if ((K <= delta) && (h_adaptive < h_max)) {

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;

				if (h_adaptive * pow(abs(delta / K), 0.2) < h_max) {
					h_adaptive = h_adaptive * pow(abs(delta / K), 0.2);
				}

				break;

			}
			if (h_adaptive <= h_min) {

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				break;

			}

			//精度が満たされなかった場合stepsize hをきつくする
			h_adaptive = h_adaptive * pow(abs(delta / K), 0.2);
		}

		if (t[i] >= t_end) {
			//tがt_endを超えた時点でループから抜ける
			for (int j = i; j < n_ini; j++) {
				t[j] = t[i];
			}
			break;
		}

		if (Flag_cout == 1) {
			for (int l = 0; l < N; l++) {
				cout << u[i][l];
				cout << " ";
			}
			cout << h_adaptive;
			cout << "\n";
		}
	}

	//ファイルに書き出し
	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_ERK45");
		Output_h(t, h, "Result_ERK45_stepsize");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_ERK45_parameter");
		Output_h(t, h, "Result_ERK45_stepsize");
	}

}

//Adaptive stepsize RKF23 Optimazed
void Solvers::ERK23_O() {

	vector<vector<double>> u, A; // AはRunge-Kutta matrix
	vector<double> t, x, k2, k3, b1, b2, c; //xはError control用のベクトル b, c はそれぞれRunge-Kutta wight, Runge-Kutta node
	vector<double> h; //adaptive stepsize h_adaptiveの記録用
	double K, h_adaptive = h_start; //可変h

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	A = vector<vector<double>>(3, vector<double>(3, 0));
	t = vector<double>(n_ini, 0);
	h = vector<double>(n_ini, 0);
	x = vector<double>(N, 0);
	k2 = vector<double>(N, 0);
	k3 = vector<double>(N, 0);
	b1 = vector<double>(2, 0);
	b2 = vector<double>(3, 0);
	c = vector<double>(3, 0);

	t[0] = t_start;
	h[0] = 0;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//Rk matrix, wight, node を定義
	A[1][0] = 2. / 3.;
	A[2][1] = 2. / 3.;

	c[1] = 2. / 3.;
	c[2] = 2. / 3.;

	b1[0] = 1. / 4.;
	b1[1] = 3. / 4.;

	b2[0] = 1. / 4.;
	b2[1] = 3. / 8.;
	b2[2] = 3. / 8.;

	cout << "ERK23_O\n\n";

	//以下ERK23法の計算
	for (int i = 1; i < n_ini; i++) {

		for (;;) {
			//kの計算
			for (int j = 0; j < N; j++) {

				k2[j] = u[i - 1][j] + h_adaptive * A[1][0] * (this->*Func[j])(t[i - 1], u[i - 1]);

			}
			for (int j = 0; j < N; j++) {

				k3[j] = u[i - 1][j] + h_adaptive * A[2][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[2][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

			}

			//誤差の計算
			vector<double> Temp(N, 0);
			double q1, q2;
			q1 = 2. / 3.;
			q2 = 3. / 8.;

			for (int j = 0; j < N; j++) {

				Temp[j] = (this->*Func[j])(t[i - 1] + q1 * h_adaptive, k3) - (this->*Func[j])(t[i - 1] + q1 * h_adaptive, k2);
				Temp[j] *= q2 * h_adaptive;

			}

			//誤差ノルムを計算
			K = L_p_norm(Temp, 2);

			if ((K <= delta) && (h_adaptive < h_max)) {

				//u[i] の計算
				for (int j = 0; j < N; j++) {

					u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
						+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

				}

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				double Temp2 = h_adaptive * pow(abs(delta / K), 1. / 3.);

				if (Temp2 < h_max) {
					h_adaptive = Temp2;
				}

				break;

			}if (h_adaptive <= h_min) {

				//u[i] の計算
				for (int j = 0; j < N; j++) {

					u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
						+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

				}

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				break;

			}

			//精度が満たされなかった場合stepsize hをきつくする
			h_adaptive = h_adaptive * pow(abs(delta / K), 1. / 3.);
		}

		if (t[i] >= t_end) {
			//tがt_endを超えた時点でループから抜ける
			for (int j = i; j < n_ini; j++) {
				t[j] = t[i];
			}
			break;
		}

		if (Flag_cout == 1) {
			for (int l = 0; l < N; l++) {
				cout << u[i][l];
				cout << " ";
			}
			cout << h_adaptive;
			cout << "\n";
		}
	}

	//ファイルに書き出し
	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_ERK23_O");
		Output_h(t, h, "Result_ERK23_stepsize_O");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_ERK23_O_parameter");
		Output_h(t, h, "Result_ERK23_O_stepsize");
	}

}

//Adaptive stepsize RKF45 Optimazed
void Solvers::ERK45_O() {

	vector<vector<double>> u, A; // AはRunge-Kutta matrix
	vector<double> t, x, k2, k3, k4, k5, k6, b1, b2, c; //xはError control用のベクトル b, c はそれぞれRunge-Kutta wight, Runge-Kutta node
	vector<double> h; //adaptive stepsize h_adaptiveの記録用
	double K, h_adaptive = h_start; //可変h

	u = vector<vector<double>>(n_ini, vector<double>(N, 0));
	A = vector<vector<double>>(6, vector<double>(6, 0));
	t = vector<double>(n_ini, 0);
	h = vector<double>(n_ini, 0);
	x = vector<double>(N, 0);
	k2 = vector<double>(N, 0);
	k3 = vector<double>(N, 0);
	k4 = vector<double>(N, 0);
	k5 = vector<double>(N, 0);
	k6 = vector<double>(N, 0);
	b1 = vector<double>(5, 0);
	b2 = vector<double>(6, 0);
	c = vector<double>(6, 0);

	t[0] = t_start;
	h[0] = 0;

	//以下u_iniを使ってuを初期化
	for (int i = 0; i < N; i++) {
		u[0][i] = u_ini[i];
	}

	//Rk matrix, wight, node を定義
	A[1][0] = 1. / 4.;
	A[2][0] = 3. / 32.;
	A[3][0] = 1932. / 2197.;
	A[4][0] = 439. / 216.;
	A[5][0] = -8. / 27.;
	A[2][1] = 9. / 32.;
	A[3][1] = -7200. / 2197.;
	A[4][1] = -8.;
	A[5][1] = 2.;
	A[3][2] = 7296. / 2197.;
	A[4][2] = 3680. / 513.;
	A[5][2] = -3544. / 2565;
	A[4][3] = -845. / 4104;
	A[5][3] = 1859. / 4104.;
	A[5][4] = -11. / 40.;

	c[1] = 1. / 4.;
	c[2] = 3. / 8.;
	c[3] = 12. / 13.;
	c[4] = 1.;
	c[5] = 1. / 2.;

	b1[0] = 25. / 216.;
	b1[2] = 1408. / 2565.;
	b1[3] = 2197. / 4104.;
	b1[4] = -1. / 5.;

	b2[0] = 16. / 135.;
	b2[2] = 6656. / 12825.;
	b2[3] = 28561. / 56430.;
	b2[4] = -9. / 50;
	b2[5] = 2. / 55.;

	cout << "ERK45_O\n\n";

	//以下ERK45法の計算
	for (int i = 1; i < n_ini; i++) {

		int count = 0;

		for (;;) {
			//kの計算
			for (int j = 0; j < N; j++) {

				k2[j] = u[i - 1][j] + h_adaptive * A[1][0] * (this->*Func[j])(t[i - 1], u[i - 1]);

			}
			for (int j = 0; j < N; j++) {

				k3[j] = u[i - 1][j] + h_adaptive * A[2][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[2][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2);

			}
			for (int j = 0; j < N; j++) {

				k4[j] = u[i - 1][j] + h_adaptive * A[3][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[3][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[3][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3);

			}
			for (int j = 0; j < N; j++) {

				k5[j] = u[i - 1][j] + h_adaptive * A[4][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[4][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[4][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * A[4][3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4);

			}
			for (int j = 0; j < N; j++) {

				k6[j] = u[i - 1][j] + h_adaptive * A[5][0] * (this->*Func[j])(t[i - 1], u[i - 1])
					+ h_adaptive * A[5][1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
					+ h_adaptive * A[5][2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
					+ h_adaptive * A[5][3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
					+ h_adaptive * A[5][4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5);

			}

			//誤差ベクトルを計算
			vector<double> Temp(N, 0);

			for (int j = 0; j < N; j++) {

				Temp[j] = -1 * (this->*Func[j])(t[i - 1], u[i - 1]) / 360
					+ 128 * (this->*Func[j])(t[i - 1] + 3 * h_adaptive / 8, k3) / 4275
					+ 2197 * (this->*Func[j])(t[i - 1] + 12 * h_adaptive / 13, k4) / 75240
					- 1 * (this->*Func[j])(t[i - 1] + h_adaptive, k5) / 50
					- 2 * (this->*Func[j])(t[i - 1] + h_adaptive / 2, k6) / 55;
				Temp[j] *= h_adaptive;

			}

			//誤差ノルムを計算
			K = L_p_norm(Temp, 2);

			if ((K <= delta) && (h_adaptive < h_max)) {

				//u[i], x の計算
				for (int j = 0; j < N; j++) {

					u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
						+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
						+ h_adaptive * b1[2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
						+ h_adaptive * b1[3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
						+ h_adaptive * b1[4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5);

				}

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				double Temp2 = h_adaptive * pow(abs(delta / K), 0.2);

				if (Temp2 < h_max) {
					h_adaptive = Temp2;
				}

				break;

			}
			count++;
			if ((h_adaptive <= h_min) || (count == 100)) {

				for (int j = 0; j < N; j++) {

					u[i][j] = u[i - 1][j] + h_adaptive * b1[0] * (this->*Func[j])(t[i - 1], u[i - 1])
						+ h_adaptive * b1[1] * (this->*Func[j])(t[i - 1] + c[1] * h_adaptive, k2)
						+ h_adaptive * b1[2] * (this->*Func[j])(t[i - 1] + c[2] * h_adaptive, k3)
						+ h_adaptive * b1[3] * (this->*Func[j])(t[i - 1] + c[3] * h_adaptive, k4)
						+ h_adaptive * b1[4] * (this->*Func[j])(t[i - 1] + c[4] * h_adaptive, k5);

				}

				t[i] = t[i - 1] + h_adaptive;
				h[i] = h_adaptive;
				break;

			}

			//精度が満たされなかった場合stepsize hをきつくする
			h_adaptive = h_adaptive * pow(abs(delta / K), 0.2);
		}

		if (t[i] >= t_end) {
			//tがt_endを超えた時点でループから抜ける
			for (int j = i; j < n_ini; j++) {
				t[j] = t[i];
			}
			break;
		}

		if (Flag_cout == 1) {
			for (int l = 0; l < N; l++) {
				cout << u[i][l];
				cout << " ";
			}
			cout << h_adaptive;
			cout << "\n";
		}
	}

	//ファイルに書き出し
	if (Flag_write_xt == 1) {
		Output_xt(t, u, "Result_ERK45_O");
		Output_h(t, h, "Result_ERK45_stepsize_O");
	}
	if (Flag_write_parameter_t == 1) {
		Output_write_parameter_t(u, "Result_ERK45_parameter_O");
		Output_h(t, h, "Result_ERK45_stepsize_O");
	}

}