#include <iostream>
#include <list>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

int compareSum(const string& file1, const string& file2, const double& tolerance)
{
    int compright = 0, compleft = 0;
    double ecart = 0;

    list<double> valleft, valright;
    list<double>::iterator itleft, itright;

    ifstream ifs;

    // Distinction entre fichiers binaire (.bin) et ASCII (non .bin)

    string binary = ".bin";

    // Extension

    int rank_file1 = (int) file1.find_last_of(".");
    string extend_file1 = string(file1, rank_file1, file1.size());

    int rank_file2 = (int) file2.find_last_of(".");
    string extend_file2 = string(file2, rank_file2, file2.size());

    // cout << "rank_file1 = " << rank_file1 << " " << extend_file1 << endl;
    // cout << "rank_file2 = " << rank_file2 << " " << extend_file2 << endl;

    if (extend_file1.compare(extend_file2)) {
        cout << "The two files are probably not of the same type !" << endl;
        cout << "Extension of file " << file1 << " is " << extend_file1 << endl;
        cout << "Extension of file " << file2 << " is " << extend_file2 << endl;
        return (1);
    }

    // Lecture du fichier gauche

    if (!extend_file1.compare(string(binary))) {
        // Ouverture en binaire
        ifs.open(file1.c_str(), ios::in | ios::binary);
    } else {
        // Ouverture en ASCII
        ifs.open(file1.c_str(), ios::in);
    }

    if (!extend_file1.compare(string(binary))) {
        unsigned long long val_long;
        int val_int;
        float val_float;

        while (ifs) {
            ifs.read(reinterpret_cast<char *>(&val_long), sizeof(val_long));
            valleft.push_back(static_cast<double>(val_long));

            ifs.read(reinterpret_cast<char *>(&val_int), sizeof(val_int));
            valleft.push_back(static_cast<double>(val_int));

            ifs.read(reinterpret_cast<char *>(&val_float), sizeof(val_float));
            valleft.push_back(static_cast<double>(val_float));

            // cout << val_long << " " << val_int << " " << val_float << endl;

            compleft++;
        }
    } else {
        string valtmp;

        while (!ifs.eof()) {
            ifs >> valtmp;

            if (valtmp.find('0') || valtmp.find('1') || valtmp.find('2') ||
                valtmp.find('3') || valtmp.find('4') || valtmp.find('5') ||
                valtmp.find('6') || valtmp.find('7') || valtmp.find('8') ||
                valtmp.find('9'))
                valleft.push_back(atof(valtmp.c_str()));

            compleft++;
        }
    }

    ifs.close();

    if (compleft == 0) {
        cout << "File " << file1 << " is empty!" << endl;
        return 1;
    }

    valleft.sort();

    // Lecture du fichier droit

    if (!extend_file2.compare(string(binary))) {
        // Ouverture en binaire
        ifs.open(file2.c_str(), ios::in | ios::binary);
    } else {
        // Ouverture en ASCII
        ifs.open(file2.c_str(), ios::in);
    }

    if (!extend_file2.compare(string(binary))) {
        unsigned long long val_long;
        int val_int;
        float val_float;

        while (ifs) {
            ifs.read(reinterpret_cast<char *>(&val_long), sizeof(val_long));
            valright.push_back(static_cast<double>(val_long));

            ifs.read(reinterpret_cast<char *>(&val_int), sizeof(val_int));
            valright.push_back(static_cast<double>(val_int));

            ifs.read(reinterpret_cast<char *>(&val_float), sizeof(val_float));
            valright.push_back(static_cast<double>(val_float));

            // cout << val_long << " " << val_int << " " << val_float << endl;

            compright++;
        }
    } else {
        string valtmp;

        while (!ifs.eof()) {
            ifs >> valtmp;

            if (valtmp.find('0') || valtmp.find('1') || valtmp.find('2') ||
                valtmp.find('3') || valtmp.find('4') || valtmp.find('5') ||
                valtmp.find('6') || valtmp.find('7') || valtmp.find('8') ||
                valtmp.find('9'))
                valright.push_back(atof(valtmp.c_str()));

            compright++;
        }
    }

    ifs.close();

    if (compright == 0) {
        cout << "File " << file2 << " is empty!" << endl;
        return 1;
    }

    valright.sort();

    // Comparaison

    for (itleft = valleft.begin(), itright = valright.begin();
         itleft != valleft.end(); itleft++, itright++) {
        // cout << (*itright) << " " << (*itleft) << " " << ecart << endl;

        ecart += fabs((*itright) - (*itleft));
    }

    // Bilan

    if (ecart < tolerance || ecart == 0) {
        return 0;
    } else {
        cout << "ecart = " << ecart << endl;

        if (compright != compleft) {
            cout << "Difference of number of lines : file1 = " << compleft
                 << " file2 = " << compright << endl;
        }

        cout << "Gap between files = " << ecart << endl;
        cout << "Test failed !" << endl;
        return 1;
    }
}
