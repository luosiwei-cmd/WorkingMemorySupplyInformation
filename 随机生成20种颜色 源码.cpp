#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdlib>
using namespace std;

int Rand(int l, int r) {
	return (rand() % (r - l + 1) + l);
} 

int main() {
	freopen("rgb.out", "w", stdout);
	srand(time(0));
	for (int i = 1; i <= 30; i++) {
		cout << Rand(0, 255) << " " << Rand(0, 255) << " " << Rand(0, 255) << endl;
	} 
	fclose(stdout);
	return 0;
}
