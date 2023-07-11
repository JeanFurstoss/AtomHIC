#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

int compareExact(const string& file1, const string& file2, const double& tolerance)
{
    int cp_ok = 0, cp_tot = 0;
    double valeur;
    double ref_valeur;
    char mot[256];
    char ref_mot[256];
    char *fin;
    char *ref_fin;
    bool stop_after_first_error = true;

    ifstream in(file1.c_str());
    ifstream ref_in(file2.c_str());

    while (!in.eof() && !ref_in.eof()) {
        in >> mot;
        valeur = strtod(mot, &fin);

        ref_in >> ref_mot;
        ref_valeur = strtod(ref_mot, &ref_fin);

        if (!strcmp(fin, "") && !strcmp(ref_fin, "")) {
            if (!valeur && !ref_valeur) {
                cp_ok++;
            } else {
                if (fabs(valeur - ref_valeur) < tolerance || 2.0 * fabs(valeur - ref_valeur)/(fabs(valeur) + fabs(ref_valeur)) < tolerance) {
                    cp_ok++;
                } else {
                    cout << "num   : " << mot << " != " << ref_mot << endl;
                    if (stop_after_first_error) return 1;
                }
            }
        } else {
            if (!strcmp(mot, ref_mot)) {
                cp_ok++;
            } else {
                cout << "alpha : " << mot << " != " << ref_mot << endl;
                if (stop_after_first_error) return 1;
            }
        }
        cp_tot++;
    }

    while (!in.eof()) {
        in >> mot;
        valeur = strtod(mot, &fin);

        if (!strcmp(fin, "")) {
            cout << "num   : " << mot << " != " << endl;
            if (stop_after_first_error) return 1;
        } else {
            cout << "alpha : " << mot << " != " << endl;
            if (stop_after_first_error) return 1;
        }
        cp_tot++;
    }

    while (!ref_in.eof()) {
        ref_in >> ref_mot;
        ref_valeur = strtod(ref_mot, &ref_fin);

        if (!strcmp(ref_fin, "")) {
            cout << "num   : " << ref_mot << " != " << endl;
            if (stop_after_first_error) return 1;
        } else {
            cout << "alpha : " << ref_mot << " != " << endl;
            if (stop_after_first_error) return 1;
        }
        cp_tot++;
    }

    in.close();
    ref_in.close();

    return (cp_ok * 100) / (cp_tot);
}
