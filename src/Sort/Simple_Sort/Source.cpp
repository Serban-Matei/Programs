#include <iostream>
#include <vector>

using namespace std;

template<class T>
void swap_sort(T& a, T& b) {

	T temp;
	temp = a;
	a = b;
	b = temp;

}

template<class T>
void sort(T& v, int Begin, int End) {

	for (int i = Begin; i < End + 1; i++) {
		for (int j = i; j < End + 1; j++) {
			if (v[i] > v[j]) {
				swap_sort(v[i], v[j]);
			}
		}
	}

}

int main() {

	int N = 10000;
	vector<int> v(N);

	for (int i = 0; i < N; i++) {
		v[i] = rand();
	}

	sort(v, 0, N - 1);

	cout.precision(17);

	for (int i = 0; i < N; i++) {
		cout << v[i] << "\n";
	}

}