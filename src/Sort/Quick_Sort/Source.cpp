#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

template<class T>
inline void swap_sort(T& a, T& b) {

	T temp;
	temp = a;
	a = b;
	b = temp;

}

template<class T>
inline T select_middle(T a, T b, T c) {

	T v[3] = { a, b, c };
	if (v[0] > v[1]) { swap_sort(v[0], v[1]); }
	if (v[0] > v[2]) { swap_sort(v[0], v[2]); }
	if (v[1] > v[2]) { swap_sort(v[1], v[2]); }
	return v[1];

}


template<class T>
void quick_sort(T& v, int Begin, int End) {

	int i = Begin - 1;
	int j = End + 1;
	int pivot = v[Begin];
	//int pivot = v[(Begin + End) / 2];
	//int pivot = select_middle(v[Begin], v[(Begin + End) / 2], v[End]);

	if (Begin == End) {
		return;
	}

	while (1) {
		do { i++; } while (v[i] < pivot);
		do { j--; } while (pivot < v[j]);

		if (i < j) { 
			swap_sort(v[i], v[j]);
		}
		else { break; }
	}

	quick_sort(v, Begin, j);
	quick_sort(v, j + 1, End);

}

 int main() {

	int N = 20000000;
	vector<int> v(N, 0);

	for (int i = 0; i < N; i++) {
		v[i] = ((2123 * i + 1754) * i + i) % 1286;
	}

	quick_sort(v, 0, N - 1);

}