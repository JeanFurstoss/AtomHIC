#include "compare.h"

#include <iostream>
#include <string>
#include <string.h>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {
    map<string, double> toleranceFields;
    double defaultTolerance = 0.001;
    string type = "exact";
    int status = 0;

    if (argc < 3) {
        cout << "Usage : Compare FILE1 FILE2 [DefaultTolerance]" << endl;
        cout << "        DefaultTolerance : default 0.001" << endl;
        return 1;
    }

    // Arguments

    // argv[1] = first file (first file)
    string file1 = argv[1];

    ifstream f1(file1.c_str());
    if (!f1.good()) {
	cout << "FILE1 (" << file1 << ") does not exist" << endl;
	return 1;	
    }

    string extension1 = strrchr(argv[1],'.');

    // argv[2] = second file (to compare the first one)
    string file2 = argv[2];

    ifstream f2(file2.c_str());
    if (!f2.good()) {
	cout << "FILE2 (" << file2 << ") does not exist" << endl;
	return 1;	
    }

    string extension2 = strrchr(argv[2],'.');

    // argv[3] = default tolerance
    if (argc > 3) {
	defaultTolerance = strtod(argv[3], 0);
    }

    if (type == "sum") {
	status = compareSum(file1, file2, defaultTolerance);
    } else {
	status = compareExact(file1, file2, defaultTolerance);
    }

    return status;
}
