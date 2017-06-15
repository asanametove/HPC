#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main(int argc, char **argv) {
	int N = 100;
	if(argc > 1)
		N = atoi(argv[1]);
	
	ofstream fout;
	fout.open("input.txt");
	fout << N << endl;

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i == j) 
				fout << rand() % 100 + 10000 << " ";
			else 
				fout << rand() % 100 << " ";
		}
		fout << endl;
	}
	fout.close();

	return 1;
}
