#pragma once
#ifndef GNUPLOT_H_
#define GNUPLOT_H_
#endif
#include <string>
#include <iostream>

using namespace std;

class Gnuplot {

public:
	Gnuplot();
	~Gnuplot();
	void operator () (const string & command);

private:
	FILE* gnuplotpipe;

};

Gnuplot::Gnuplot() {

	gnuplotpipe = _popen("gnuplot -persist", "w");

	if (!gnuplotpipe) {
		cerr << ("Gnuplot not found !");
	}
}

Gnuplot::~Gnuplot() {

	fprintf(gnuplotpipe, "exit\n");
	_pclose(gnuplotpipe);

}

void Gnuplot::operator()(const string & command) {

	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);

}