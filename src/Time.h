#include <iostream>
#include <chrono>

class Time {

private:
	std::chrono::system_clock::time_point startA, endA;

public:
	void start() { 
		startA = std::chrono::system_clock::now(); 
	};
	void end() {
		endA = std::chrono::system_clock::now();       // 計測終了時刻を保存
		auto dur = endA - startA;        // 要した時間を計算
		auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
		// 要した時間をミリ秒（1/1000秒）に変換して表示
		std::cout << msec << " milli sec \n";
	};

};
